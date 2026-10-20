"""Rewrite per-residue CSVs into canonical UniProt residue numbering.

Reads 03_data/uniprot_mapping.json (canonical, evidence-documented) and
renumbers every per-residue file by ORDINAL POSITION within each chain
(rows sorted by the original PDB number preserved in the *_pdb columns),
which makes the mapping independent of each file's chain-numbering
convention (osr1 rmsf numbers its chains continuously 1..96/97..192/193..198/
199..203; other files restart per chain).

Rules
  protein chains  -> UniProt blocks (ordinal -> block offset)
  mapped ligands  -> position within the ligand (osr1: two copies of the
                     GRFQVT 6-mer -> WNK4 1015-1020, cycle mod 6)
  unmapped ligands-> keep the file's original PDB numbers

Files rewritten per system/replica in organized_by_system_bundle_20260327:
  csv/rmsf.csv, csv/sasa.csv, dssp/dssp_per_residue.csv
(03_data holds no loose per-replica files: everything lives in the bundle.)
rmsf residue_label is renumbered to NAME + UniProt position. Original PDB
numbering is preserved in the *_pdb columns.

Run:  py -3.13 02figure_scripts\\_apply_uniprot_numbering.py
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "03_data"
BUNDLE = DATA_DIR / "organized_by_system_bundle_20260327"
MAPPING = DATA_DIR / "uniprot_mapping.json"

PDB_COL = {"residue": "residue_pdb", "resSeq": "resseq_pdb", "resid": "resid_pdb"}


def protein_up(ordinal: int, blocks) -> int:
    for b in blocks:
        if b["pdb_start"] <= ordinal <= b["pdb_end"]:
            return ordinal + int(b["offset"])
    return ordinal  # outside documented blocks: keep as-is


def ligand_up(ordinal: int, lig) -> int | None:
    if lig.get("status") == "unmapped":
        return None
    if lig.get("mod_length"):  # repeated motif copies (OSR1 GRFQVT x2)
        return int(lig["uniprot_start"]) + (ordinal - 1) % int(lig["mod_length"])
    ranges = lig.get("ranges") or []
    for r in ranges:
        if r["pdb_start"] <= ordinal <= r["pdb_end"]:
            return ordinal + int(r["offset"])
    return None


def rewrite(path: Path, number_col: str, group_col: str, role_by_chain, blocks, lig,
            label_col: str | None = None, name_col: str = "residue_name") -> str:
    df = pd.read_csv(path)
    pdb_col = PDB_COL[number_col]
    if pdb_col not in df.columns:
        df[pdb_col] = df[number_col]

    notes = []
    new_vals = df[number_col].astype("Int64").copy()

    for chain, sub in df.groupby(group_col):
        role = role_by_chain(chain)
        order = sub.sort_values(pdb_col)
        ordinals = range(1, len(order) + 1)
        for ordinal, (row_idx, row) in zip(ordinals, order.iterrows()):
            if role == "protein":
                up = protein_up(ordinal, blocks)
            else:
                up = ligand_up(ordinal, lig)
            if up is None:
                # unmapped ligand (or row outside documented blocks):
                # restore the original PDB number
                up = int(row[pdb_col])
            new_vals.loc[row_idx] = up
        if role == "protein":
            notes.append(f"{chain}: {len(order)} protein rows")
        else:
            notes.append(f"{chain}: {len(order)} ligand rows"
                         + (" (unmapped, PDB kept)" if ligand_up(1, lig) is None else ""))

    changed = int((new_vals.astype("Int64") != df[number_col].astype("Int64")).sum())
    df[number_col] = new_vals

    if label_col and label_col in df.columns:
        for chain, sub in df.groupby(group_col):
            role = role_by_chain(chain)
            order = sub.sort_values(pdb_col)
            for ordinal, (row_idx, row) in zip(range(1, len(order) + 1), order.iterrows()):
                if pd.isna(df.at[row_idx, label_col]) or pd.isna(df.at[row_idx, number_col]):
                    continue
                df.at[row_idx, label_col] = (str(row[name_col])
                                             + str(int(df.at[row_idx, number_col])))

    cols = [c for c in df.columns if c != pdb_col]
    cols.insert(cols.index(number_col) + 1, pdb_col)
    df[cols].to_csv(path, index=False)
    return f"{changed} rows renumbered [{'; '.join(notes)}]"


def main() -> None:
    mapping = json.loads(MAPPING.read_text(encoding="utf-8"))
    prot_map = mapping["protein"]
    lig_map = mapping["ligands"]

    for slug in sorted(prot_map):
        blocks = prot_map[slug]["blocks"]
        lig = lig_map.get(slug, {"status": "unmapped"})
        for rep in ("R1", "R2"):
            sasa_path = BUNDLE / slug / rep / "csv" / "sasa.csv"
            sasa_df = pd.read_csv(sasa_path)
            role_by_label = dict(zip(
                sasa_df["chain_label"].astype(str),
                sasa_df["entity_role"].astype(str).str.lower()))

            def role_sasa(chain):
                return role_by_label.get(str(chain), "protein")

            def role_rmsf(chain):
                c = str(chain).lower()
                return "ligand" if ("pep" in c or "rfxv" in c) else "protein"

            def role_dssp(chain):
                return role_sasa(chain)

            jobs = [
                (BUNDLE / slug / rep / "csv" / "rmsf.csv", "residue", "chain", role_rmsf,
                 "residue_label", "residue_name"),
                (sasa_path, "resSeq", "chain_label", role_sasa, None, "residue_name"),
                (BUNDLE / slug / rep / "dssp" / "dssp_per_residue.csv", "resid", "chain",
                 role_dssp, None, "resname"),
            ]
            for path, col, group_col, role_fn, label_col, name_col in jobs:
                msg = rewrite(path, col, group_col, role_fn, blocks, lig, label_col, name_col)
                print(f"{slug} {rep} {path.name}: {msg}")

    print("[OK] Canonical UniProt numbering applied (protein blocks + mapped "
          "ligands); original PDB numbering preserved in *_pdb columns.")


if __name__ == "__main__":
    main()

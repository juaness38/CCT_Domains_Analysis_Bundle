#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pick the ensemble state whose contact distances are closest to the consensus.

Reads a multi-MODEL PDB and a table of reference interactions, measures every
listed interaction in every model, and selects the model with the smallest sum
of squared deviations from the reference distances. That model is the one that
best represents the ensemble's average interface geometry.

Distances are measured according to the interaction type:

  * PIPISTACK, between the centroids of the two aromatic rings;
  * IONIC, between the centroids of the two charged groups;
  * HBOND and anything else, as the shortest polar heavy-atom pair.

The CCT domain chain is renumbered to canonical numbering in the output, using
the full-length FASTA, so the selected structure can be read against the
literature. The peptide chain keeps its own numbering.

The reference table needs the columns interaction, source, target and
average_distance_mean. Residues are given as CHAIN/NUMBER/RESNAME, and the
number may be a range written START-END, in which case the closest valid pair
is used.
"""

import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

AA_POLAR = {
    "ALA": ["N", "O"], "ARG": ["N", "O", "NE", "NH1", "NH2"],
    "ASN": ["N", "O", "OD1", "ND2"], "ASP": ["N", "O", "OD1", "OD2"],
    "CYS": ["N", "O", "SG"], "GLN": ["N", "O", "OE1", "NE2"],
    "GLU": ["N", "O", "OE1", "OE2"], "GLY": ["N", "O"],
    "HIS": ["N", "O", "ND1", "NE2"], "ILE": ["N", "O"],
    "LEU": ["N", "O"], "LYS": ["N", "O", "NZ"], "MET": ["N", "O", "SD"],
    "PHE": ["N", "O"], "PRO": ["N", "O"], "SER": ["N", "O", "OG"],
    "THR": ["N", "O", "OG1"], "TRP": ["N", "O", "NE1"],
    "TYR": ["N", "O", "OH"], "VAL": ["N", "O"],
}

RING_ATOMS = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"],
}

CHARGED_GROUP = {
    "ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"],
    "ARG": ["CZ", "NH1", "NH2"], "LYS": ["NZ"], "HIS": ["ND1", "NE2"],
}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("pdb", type=Path, help="Multi-MODEL PDB of the ensemble")
    p.add_argument("table", type=Path,
                   help="Reference interactions, as .xlsx or .csv, with columns "
                        "interaction, source, target, average_distance_mean")
    p.add_argument("fasta", type=Path,
                   help="Full-length FASTA of the CCT domain protein")
    p.add_argument("--cct-domain-chain", default="A",
                   help="Chain holding the CCT domain, the one renumbered (default: A)")
    p.add_argument("--output-pdb", type=Path, default=Path("representative_state.pdb"),
                   help="Where to write the selected model (default: representative_state.pdb)")
    p.add_argument("--report", type=Path, default=Path("representative_state_report.csv"),
                   help="Per-model distances and scores (default: representative_state_report.csv)")
    p.add_argument("--renumbered-ensemble", type=Path, default=None,
                   help="Optionally also write the whole ensemble with the CCT "
                        "domain renumbered")
    return p.parse_args(argv)


def read_fasta(path):
    sequence = "".join(
        line.strip() for line in path.read_text().splitlines()
        if line.strip() and not line.startswith(">")
    ).upper()
    if not sequence:
        raise ValueError(f"No sequence found in FASTA: {path}")
    return sequence


def read_table(path):
    if path.suffix.lower() in (".csv", ".tsv", ".txt"):
        separator = "\t" if path.suffix.lower() == ".tsv" else ","
        return pd.read_csv(path, sep=separator)
    return pd.read_excel(path)


def parse_models(pdb_file):
    """Read every MODEL of a multi-MODEL PDB, keeping its raw lines."""
    models = []
    current = None
    raw_lines = None
    number = None

    with pdb_file.open() as handle:
        for line in handle:
            if line.startswith("MODEL"):
                current = {}
                raw_lines = [line]
                match = re.search(r"MODEL\s+(\d+)", line)
                number = int(match.group(1)) if match else len(models) + 1
            elif current is not None:
                raw_lines.append(line)
                if line.startswith(("ATOM  ", "HETATM")):
                    atom = line[12:16].strip()
                    resname = line[17:20].strip().upper()
                    chain = line[21].strip()
                    resi = int(line[22:26])
                    position = np.array([float(line[30:38]), float(line[38:46]),
                                         float(line[46:54])], dtype=float)
                    current.setdefault((chain, resi, resname), {})[atom] = position
                elif line.startswith("ENDMDL"):
                    models.append({"number": number, "atoms": current,
                                   "lines": raw_lines[:]})
                    current = None
                    raw_lines = None
                    number = None

    if not models:
        raise ValueError("No MODEL/ENDMDL records were found in the PDB. This "
                         "script expects the multi-MODEL ensemble, not a single "
                         "structure.")
    return models


def infer_domain_mapping(model_atoms, fasta_sequence, chain):
    """Map each PDB residue number of the domain chain to its canonical number."""
    residues = sorted(
        {(resi, resname) for (ch, resi, resname) in model_atoms if ch == chain},
        key=lambda item: item[0]
    )
    if not residues:
        raise ValueError(f"Chain {chain} was not found in the PDB. "
                         f"Check --cct-domain-chain.")
    try:
        pdb_sequence = "".join(AA3_TO_1[resname] for _, resname in residues)
    except KeyError as error:
        raise ValueError(f"Unsupported residue name in chain {chain}: "
                         f"{error.args[0]}") from error

    first = fasta_sequence.find(pdb_sequence)
    if first < 0:
        raise ValueError(f"The chain {chain} sequence was not found as an exact "
                         f"contiguous match in the FASTA. This script stops rather "
                         f"than guessing the numbering.")
    if fasta_sequence.find(pdb_sequence, first + 1) >= 0:
        raise ValueError(f"The chain {chain} sequence occurs more than once in the "
                         f"FASTA; canonical numbering is ambiguous.")

    canonical_start = first + 1
    mapping = {resi: canonical_start + i for i, (resi, _) in enumerate(residues)}
    return mapping, pdb_sequence, canonical_start, canonical_start + len(residues) - 1


def parse_residue_spec(spec, canonical_to_pdb, chain_domain):
    """CHAIN/NUMBER/RESNAME, where NUMBER may be START-END, to candidate PDB keys."""
    parts = str(spec).strip().split("/")
    if len(parts) != 3:
        raise ValueError(f"Invalid residue specification: {spec!r}; expected "
                         f"CHAIN/NUMBER/RESNAME")
    chain, number_text, resname = parts[0].strip(), parts[1], parts[2].strip().upper()

    if "-" in number_text:
        low, high = number_text.split("-", 1)
        numbers = list(range(int(low), int(high) + 1))
    else:
        numbers = [int(number_text)]

    keys = []
    for number in numbers:
        if chain == chain_domain:
            if number not in canonical_to_pdb:
                raise ValueError(f"Canonical residue {chain}/{number}/{resname} lies "
                                 f"outside the modelled CCT domain fragment.")
            number = canonical_to_pdb[number]
        keys.append((chain, number, resname))
    return keys


def centroid(atoms, names):
    points = [atoms[name] for name in names if name in atoms]
    if not points:
        raise KeyError(f"None of these atoms were found: {names}")
    return np.mean(np.vstack(points), axis=0)


def distance_between(model_atoms, interaction, key1, key2):
    for key in (key1, key2):
        if key not in model_atoms:
            raise KeyError(f"Missing residue in PDB: {key}")
    atoms1, atoms2 = model_atoms[key1], model_atoms[key2]
    interaction = interaction.upper()

    if interaction == "PIPISTACK":
        if key1[2] not in RING_ATOMS or key2[2] not in RING_ATOMS:
            raise ValueError(f"PIPISTACK needs two aromatic residues: {key1}, {key2}")
        return float(np.linalg.norm(centroid(atoms1, RING_ATOMS[key1[2]]) -
                                    centroid(atoms2, RING_ATOMS[key2[2]])))

    if interaction == "IONIC":
        if key1[2] not in CHARGED_GROUP or key2[2] not in CHARGED_GROUP:
            raise ValueError(f"IONIC needs two charged residues: {key1}, {key2}")
        return float(np.linalg.norm(centroid(atoms1, CHARGED_GROUP[key1[2]]) -
                                    centroid(atoms2, CHARGED_GROUP[key2[2]])))

    names1 = AA_POLAR.get(key1[2], ["N", "O"])
    names2 = AA_POLAR.get(key2[2], ["N", "O"])
    distances = [float(np.linalg.norm(atoms1[n1] - atoms2[n2]))
                 for n1 in names1 if n1 in atoms1
                 for n2 in names2 if n2 in atoms2]
    if not distances:
        raise ValueError(f"No polar atom pair found for {key1} / {key2}")
    return min(distances)


def distance_for_interaction(model_atoms, interaction, source, target,
                             canonical_to_pdb, chain_domain):
    sources = parse_residue_spec(source, canonical_to_pdb, chain_domain)
    targets = parse_residue_spec(target, canonical_to_pdb, chain_domain)
    values = []
    for key1 in sources:
        for key2 in targets:
            try:
                values.append(distance_between(model_atoms, interaction, key1, key2))
            except KeyError:
                # A range may list a residue number whose identity does not match
                # the label; skip that candidate but require one valid pair overall.
                continue
    if not values:
        raise KeyError(f"No valid residue pair found for {source} / {target}")
    return min(values)


def renumber_atom_line(line, mapping, chain):
    if not line.startswith(("ATOM  ", "HETATM", "TER   ")) or len(line) < 26:
        return line
    if line[21].strip() != chain:
        return line
    try:
        old = int(line[22:26])
    except ValueError:
        return line
    if old not in mapping:
        return line
    return line[:22] + f"{mapping[old]:4d}" + line[26:]


def write_single_model(model, path, mapping, chain):
    with path.open("w") as out:
        for line in model["lines"]:
            if line.startswith(("MODEL", "ENDMDL")):
                continue
            out.write(renumber_atom_line(line, mapping, chain))
        if not model["lines"][-1].startswith("END"):
            out.write("END\n")


def write_ensemble(models, path, mapping, chain):
    with path.open("w") as out:
        for model in models:
            for line in model["lines"]:
                out.write(renumber_atom_line(line, mapping, chain))
        out.write("END\n")


def main(argv=None):
    args = parse_args(argv)
    chain = args.cct_domain_chain.strip().upper()

    for path in (args.pdb, args.table, args.fasta):
        if not path.exists():
            sys.exit(f"File not found: {path}")

    table = read_table(args.table)
    required = {"interaction", "source", "target", "average_distance_mean"}
    missing = required - set(table.columns)
    if missing:
        sys.exit(f"The reference table is missing columns: {sorted(missing)}")
    if table.empty:
        sys.exit("The reference table holds no rows, so no model can be scored.")

    try:
        models = parse_models(args.pdb)
        fasta_sequence = read_fasta(args.fasta)
        mapping, pdb_sequence, canonical_start, canonical_end = infer_domain_mapping(
            models[0]["atoms"], fasta_sequence, chain)
    except ValueError as error:
        sys.exit(str(error))

    canonical_to_pdb = {canonical: pdb for pdb, canonical in mapping.items()}
    offset = canonical_start - min(mapping)

    print(f"Models read: {len(models)}")
    print(f"Chain {chain} exact FASTA match: {canonical_start}-{canonical_end} "
          f"({len(pdb_sequence)} residues), offset {offset:+d}")
    if offset == 0:
        print("  The domain already carries canonical numbering; renumbering "
              "changes no residue number.")

    try:
        first_atoms = models[0]["atoms"]
        for column in ("source", "target"):
            for spec in table[column].astype(str):
                for key in parse_residue_spec(spec, canonical_to_pdb, chain):
                    if key[0] == chain and key not in first_atoms:
                        raise ValueError(f"Residue identity mismatch between the "
                                         f"table and the PDB: {spec} -> {key}")

        rows = []
        records = list(table.itertuples(index=False))
        for model in models:
            distances, residuals = [], []
            for record in records:
                value = distance_for_interaction(
                    model["atoms"], record.interaction, record.source,
                    record.target, canonical_to_pdb, chain)
                distances.append(value)
                residuals.append(value - float(record.average_distance_mean))

            sse = sum(r * r for r in residuals)
            row = {"model": model["number"], "SSE": sse,
                   "RMSE": math.sqrt(sse / len(residuals))}
            for i, (record, value) in enumerate(zip(records, distances), start=1):
                row[f"d{i}_{record.interaction}_{record.source}_{record.target}"] = value
            rows.append(row)
    except (ValueError, KeyError) as error:
        sys.exit(str(error))

    report = pd.DataFrame(rows).sort_values(["SSE", "model"])
    best_number = int(report.iloc[0]["model"])
    best_model = next(m for m in models if m["number"] == best_number)

    report.to_csv(args.report, index=False)
    write_single_model(best_model, args.output_pdb, mapping, chain)
    if args.renumbered_ensemble is not None:
        write_ensemble(models, args.renumbered_ensemble, mapping, chain)

    print(f"\nSelected MODEL: {best_number}")
    print(f"SSE:  {report.iloc[0]['SSE']:.6f}")
    print(f"RMSE: {report.iloc[0]['RMSE']:.6f} A")
    print(f"Selected structure: {args.output_pdb}")
    print(f"Report: {args.report}")
    if args.renumbered_ensemble is not None:
        print(f"Renumbered ensemble: {args.renumbered_ensemble}")


if __name__ == "__main__":
    main()

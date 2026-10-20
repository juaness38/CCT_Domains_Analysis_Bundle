"""Per-system composite panel figures: RMSD | RMSF | SASA  (R1 top, R2 bottom).

Two versions per system:
  V1 — DSSP secondary structure annotation track coupled to RMSF & SASA x-axes
  V2 — Minimal / professional without DSSP

Panels show protein chain (green) + peptide/ligand chain (orange) on the same axis
for both RMSF and SASA, matching the original residue SASA visual style.

Outputs: SVG + PNG (dpi 600) + DOCX (one page per system).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.patches import FancyArrow, FancyBboxPatch
from matplotlib.colors import to_rgba, LinearSegmentedColormap
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from docx import Document
from docx.enum.section import WD_ORIENTATION
from docx.enum.table import WD_ALIGN_VERTICAL, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt

import _replot_md_controls_rmsd_rmsf as base

# ── Paths (canonical package layout — relative to this file) ─────────
BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / "03_data"
RMSD_MANIFEST = DATA_DIR / "replot_manifest.csv"
DSSP_ROOT = DATA_DIR / "organized_by_system_bundle_20260327"


def _sasa_csv(slug: str, replica_label: str) -> Path:
    return DSSP_ROOT / slug / replica_label / "csv" / "sasa.csv"

OUT_DIR = BASE_DIR / "01_figures" / "panels"
OUT_DOCX = OUT_DIR / "md_controls_per_system_panels.docx"

TIME_CUTOFF_NS = 200.0

# ── Style ─────────────────────────────────────────────────────────────
SYSTEM_COLORS: Dict[str, str] = {
    "OSR1 CCT": "#6DB7E8",
    "SPAK CCT": "#6FAA55",
    "WNK2 CCTL1": "#C68A4B",
    "NRBP1 CCT": "#8E7BB7",
    "FERRY CCTL": "#D07286",
    "WNK4 CCTL2": "#4EA39B",
}
DISPLAY_LABELS: Dict[str, str] = {"NRBP1 CCT": "NRBP1 CCTL"}

DSSP_PALETTE = {"helix": "#E07070", "strand": "#5090D0", "coil": "#D0D0D0"}
PROTEIN_COLOR = "#90D090"
PEPTIDE_COLOR = "#F28E2B"
GRID_COLOR = "#D9E5EE"

SPEC_BY_KEY = {(s.system_label, s.replica_label): s for s in base.REPLICA_SPECS}

# ── Unified slug mapping ─────────────────────────────────────────────
_SLUGS: Dict[str, str] = {
    "OSR1 CCT": "osr1_cct",
    "SPAK CCT": "spak_cct",
    "WNK2 CCTL1": "wnk2_cctl1",
    "NRBP1 CCT": "nrbp1_cct",
    "FERRY CCTL": "ferry_cctl",
    "WNK4 CCTL2": "wnk4_cctl2",
}

# ── UniProt numbering ────────────────────────────────────────────────
# All per-residue data is ALREADY in canonical UniProt numbering for protein
# chains (applied by _apply_uniprot_numbering.py from uniprot_mapping.json);
# labels/axes use the data numbers directly — no plot-time offset.
with open(DATA_DIR / "uniprot_mapping.json", encoding="utf-8") as _fh:
    _UNIPROT_MAP = json.load(_fh)
UNIPROT_ACCESSIONS: Dict[str, str] = {
    slug: str(info["accession"]) for slug, info in _UNIPROT_MAP["protein"].items()
}

# ── OSR1 dimer: only show physically paired chains (pair A↔D) ────────
_DIMER_SYSTEMS: Dict[str, Dict[str, List[str]]] = {
    "OSR1 CCT": {
        "protein_chains": ["A"],         # first monomer
        "peptide_chains": ["D"],         # physically bound RFxV ligand
    },
}

SYSTEM_ORDER: List[str] = [
    "OSR1 CCT", "SPAK CCT", "WNK2 CCTL1",
    "NRBP1 CCT", "FERRY CCTL", "WNK4 CCTL2",
]

# ── Per-system RMSF gradient cap (Å) — avoids IDR saturation ────────
RMSF_GRADIENT_CAP: Dict[str, float] = {
    "FERRY CCTL": 4.0,
    "OSR1 CCT": 4.0,
}


# =====================================================================
#  Data loaders
# =====================================================================
def _load_manifest() -> pd.DataFrame:
    frame = pd.read_csv(RMSD_MANIFEST)
    frame["time_end_ns"] = pd.to_numeric(frame["time_end_ns"], errors="coerce")
    return frame


def _load_rmsd(system_label: str, replica_label: str) -> pd.DataFrame:
    """Return protein + peptide RMSD time-series.

    Columns: time_ns, rmsd_A (protein), rmsd_peptide_A (peptide).
    """
    manifest = _load_manifest()
    row = manifest[(manifest["system_label"] == system_label) & (manifest["replica_label"] == replica_label)]
    if row.empty:
        raise FileNotFoundError(f"No manifest entry for {system_label} {replica_label}")
    rmsd_csv = Path(str(row.iloc[0]["rmsd_csv"]))
    spec = SPEC_BY_KEY[(system_label, replica_label)]
    if rmsd_csv.exists():
        rmsd_df = pd.read_csv(str(rmsd_csv))
        prepared, _, _ = base._prepare_rmsd_time_axis(spec, rmsd_df)
    else:
        # Fallback: organized bundle on C: (already has time_ns, same columns)
        slug = _SLUGS[system_label]
        bundle_csv = DSSP_ROOT / slug / replica_label / "csv" / "rmsd.csv"
        rmsd_df = pd.read_csv(str(bundle_csv))
        rmsd_df["_plot_time_ns"] = rmsd_df["time_ns"]
        prepared = rmsd_df
    col = "rmsd_protein_ca_A" if "rmsd_protein_ca_A" in prepared.columns else "rmsd_protein_backbone_A"
    out: Dict[str, Any] = {
        "time_ns": pd.to_numeric(prepared["_plot_time_ns"], errors="coerce"),
        "rmsd_A": pd.to_numeric(prepared[col], errors="coerce"),
    }
    # Peptide RMSD (backbone)
    pep_col = "rmsd_peptide_backbone_A"
    if pep_col in prepared.columns:
        out["rmsd_peptide_A"] = pd.to_numeric(prepared[pep_col], errors="coerce")
    series = pd.DataFrame(out).dropna(subset=["time_ns", "rmsd_A"])
    return series[series["time_ns"] <= TIME_CUTOFF_NS].copy()


def _load_rmsf(system_label: str, replica_label: str) -> pd.DataFrame:
    """Load RMSF with protein + peptide chains, sequential plot_x.

    Returns DataFrame with: residue, residue_name, residue_label, rmsf_A,
    entity_role, chain_label, plot_x, tick_label.
    """
    slug = _SLUGS[system_label]
    rep = replica_label.lower()
    csv_path = DSSP_ROOT / slug / replica_label / "csv" / "rmsf.csv"
    df = pd.read_csv(csv_path)

    chains_lower = df["chain"].astype(str).str.lower()
    protein_df = df[chains_lower.isin(["protein_ca", "protein"])].copy()
    peptide_df = df[chains_lower.isin(["peptide_ca", "peptide"])].copy()

    # Legacy per-chain layout (e.g. OSR1: cct_a_ca / rfxv_c_ca / rfxv_d_ca):
    # peptide = fragment chains (contiguous residue ranges concatenated),
    # protein = the first non-peptide chain.
    if protein_df.empty and peptide_df.empty:
        is_peptide = chains_lower.str.contains("pep|rfxv", regex=True)
        peptide_df = df[is_peptide].copy()
        protein_chains = df[~is_peptide]
        if not protein_chains.empty:
            first = protein_chains["chain"].iloc[0]
            protein_df = protein_chains[protein_chains["chain"] == first].copy()

    if protein_df.empty:
        first = df.iloc[0]["chain"]
        protein_df = df[df["chain"] == first].copy()

    protein_df["entity_role"] = "protein"
    peptide_df["entity_role"] = "peptide"

    for chunk in [protein_df, peptide_df]:
        chunk["residue"] = pd.to_numeric(chunk["residue"], errors="coerce")
        chunk["rmsf_A"] = pd.to_numeric(chunk["rmsf_A"], errors="coerce")
    protein_df = protein_df.dropna(subset=["residue", "rmsf_A"]).sort_values("residue").reset_index(drop=True)
    peptide_df = peptide_df.dropna(subset=["residue", "rmsf_A"]).sort_values("residue").reset_index(drop=True)

    # For dimer systems, filter peptide to the physically paired chain
    if system_label in _DIMER_SYSTEMS and not peptide_df.empty:
        pair = _DIMER_SYSTEMS[system_label]
        sasa_csv = _sasa_csv(slug, replica_label)
        if sasa_csv.exists():
            sasa_df = pd.read_csv(sasa_csv)
            paired_pep = sasa_df[
                (sasa_df["entity_role"].astype(str).str.lower() == "peptide") &
                (sasa_df["chain_label"].isin(pair["peptide_chains"]))
            ]
            if not paired_pep.empty:
                # compare on the ORIGINAL PDB numbering (*_pdb) so the paired-
                # chain selection is invariant to the UniProt renumbering
                pdb_col = "resseq_pdb" if "resseq_pdb" in paired_pep.columns else "resSeq"
                pep_pdb_col = "residue_pdb" if "residue_pdb" in peptide_df.columns else "residue"
                min_res = int(paired_pep[pdb_col].min())
                max_res = int(paired_pep[pdb_col].max())
                peptide_df = peptide_df[
                    (peptide_df[pep_pdb_col] >= min_res) &
                    (peptide_df[pep_pdb_col] <= max_res)
                ].reset_index(drop=True)

    # Assign chain_label
    if system_label in _DIMER_SYSTEMS:
        pair = _DIMER_SYSTEMS[system_label]
        protein_df["chain_label"] = pair["protein_chains"][0]
        if not peptide_df.empty:
            peptide_df["chain_label"] = pair["peptide_chains"][0]
    else:
        protein_df["chain_label"] = "A"
        if not peptide_df.empty:
            peptide_df["chain_label"] = "B"

    # Build sequential plot_x across protein → peptide
    # Residue numbers are already canonical UniProt (protein chain) in the data
    parts: List[pd.DataFrame] = []
    idx = 0
    for role, chunk in [("protein", protein_df), ("peptide", peptide_df)]:
        if chunk.empty:
            continue
        chunk = chunk.copy()
        chunk["plot_x"] = range(idx, idx + len(chunk))
        if role == "protein":
            chunk["tick_label"] = (chunk["residue_name"].astype(str)
                                   + chunk["residue"].astype(int).astype(str))
        elif "residue_label" in chunk.columns and not chunk["residue_label"].isna().all():
            chunk["tick_label"] = chunk["residue_label"].astype(str)
        else:
            chunk["tick_label"] = chunk["residue_name"].astype(str) + chunk["residue"].astype(int).astype(str)
        parts.append(chunk)
        idx += len(chunk)

    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def _load_sasa(system_label: str, replica_label: str) -> pd.DataFrame:
    """Load SASA with protein + peptide chains, sequential plot_x.

    Returns DataFrame with: resSeq, residue_name, mean_sasa_nm2, std_sasa_nm2,
    entity_role, chain_label, plot_x, tick_label.
    """
    slug = _SLUGS[system_label]
    rep = replica_label.lower()
    csv_path = _sasa_csv(slug, replica_label)
    df = pd.read_csv(csv_path)
    df["mean_sasa_nm2"] = pd.to_numeric(df["mean_sasa_nm2"], errors="coerce")
    df["std_sasa_nm2"] = pd.to_numeric(df["std_sasa_nm2"], errors="coerce").fillna(0)
    df["resSeq"] = pd.to_numeric(df["resSeq"], errors="coerce")
    df["entity_role"] = df["entity_role"].fillna("unknown").astype(str)
    df["chain_label"] = df["chain_label"].fillna("?").astype(str)
    df = df.dropna(subset=["mean_sasa_nm2", "resSeq"])

    # For dimer systems, filter to physically paired chains only
    if system_label in _DIMER_SYSTEMS:
        pair = _DIMER_SYSTEMS[system_label]
        keep_chains = pair["protein_chains"] + pair["peptide_chains"]
        df = df[df["chain_label"].isin(keep_chains)].copy()

    protein_df = df[df["entity_role"].str.lower() == "protein"].sort_values("resSeq").reset_index(drop=True)
    peptide_df = df[df["entity_role"].str.lower() == "peptide"].sort_values("resSeq").reset_index(drop=True)

    # Residue numbers are already canonical UniProt (protein chain) in the data
    parts: List[pd.DataFrame] = []
    idx = 0
    for role, chunk in [("protein", protein_df), ("peptide", peptide_df)]:
        if chunk.empty:
            continue
        chunk = chunk.copy()
        chunk["plot_x"] = range(idx, idx + len(chunk))
        if role == "protein":
            chunk["tick_label"] = (chunk["residue_name"].astype(str)
                                   + chunk["resSeq"].astype(int).astype(str))
        else:
            chunk["tick_label"] = chunk["residue_name"].astype(str) + chunk["resSeq"].astype(int).astype(str)
        parts.append(chunk)
        idx += len(chunk)

    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def _load_dssp(system_label: str, replica_label: str) -> pd.DataFrame:
    slug = _SLUGS[system_label]
    csv_path = DSSP_ROOT / slug / replica_label / "dssp" / "dssp_per_residue.csv"
    df = pd.read_csv(csv_path)
    first_chain = str(df.iloc[0]["chain"])
    df = df[df["chain"] == first_chain].copy()
    df["resid"] = pd.to_numeric(df["resid"], errors="coerce")
    return df.sort_values("resid").reset_index(drop=True)


# =====================================================================
#  Helpers
# =====================================================================
def _ma(values, window: int = 31) -> list:
    return base._moving_average(list(values), window)


def _style_axis(ax: plt.Axes) -> None:
    ax.grid(True, color=GRID_COLOR, linewidth=0.6, alpha=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_facecolor("white")
    ax.tick_params(labelsize=7)


def _display(system_label: str) -> str:
    return DISPLAY_LABELS.get(system_label, system_label)


def _color(system_label: str) -> str:
    dl = _display(system_label)
    return SYSTEM_COLORS.get(dl, SYSTEM_COLORS.get(system_label, "#888888"))


def _slugify(value: str) -> str:
    return "_".join(value.lower().replace("/", " ").split())


# =====================================================================
#  DSSP annotation track
# =====================================================================
def _draw_dssp_track(ax: plt.Axes, dssp_df: pd.DataFrame,
                     resid_to_plotx: Optional[Dict[int, int]] = None,
                     y_center: float = 0, height: float = 0.6) -> None:
    """Draw secondary structure elements as colored bars.

    Helices → red rounded rectangles, Strands → blue filled arrows,
    Coils → thin grey line.

    *resid_to_plotx*: optional mapping from DSSP resid → sequential plot_x
    for alignment with multi-chain panels. If None, raw resid is used.
    """
    if dssp_df.empty:
        return

    dssp_df = dssp_df.copy()
    if resid_to_plotx is not None:
        dssp_df["_x"] = dssp_df["resid"].map(resid_to_plotx)
        dssp_df = dssp_df.dropna(subset=["_x"])
        dssp_df["_x"] = dssp_df["_x"].astype(int)
    else:
        dssp_df["_x"] = dssp_df["resid"].astype(int)

    if dssp_df.empty:
        return

    # Group contiguous runs of the same ss_type
    runs: List[Tuple[str, int, int]] = []
    current_ss = str(dssp_df.iloc[0]["ss_type"])
    start_x = int(dssp_df.iloc[0]["_x"])
    end_x = start_x
    prev_x = start_x
    for _, row in dssp_df.iloc[1:].iterrows():
        ss = str(row["ss_type"])
        x = int(row["_x"])
        if ss == current_ss and x == prev_x + 1:
            end_x = x
        else:
            runs.append((current_ss, start_x, end_x))
            current_ss = ss
            start_x = x
            end_x = x
        prev_x = x
    runs.append((current_ss, start_x, end_x))

    half_h = height / 2.0
    for ss_type, x_start, x_end in runs:
        color = DSSP_PALETTE.get(ss_type, DSSP_PALETTE["coil"])
        span = x_end - x_start + 1
        if ss_type == "helix":
            patch = FancyBboxPatch(
                (x_start - 0.5, y_center - half_h),
                span, height,
                boxstyle="round,pad=0.1",
                facecolor=color, edgecolor="none", alpha=0.85,
                zorder=3,
            )
            ax.add_patch(patch)
        elif ss_type == "strand":
            # Filled arrow occupying the full zone
            head_len = min(1.5, span * 0.35)
            arrow = FancyArrow(
                x_start - 0.5, y_center,
                span, 0,
                width=height * 0.85,
                head_width=height * 1.2,
                head_length=head_len,
                fc=color, ec="none",
                alpha=0.85,
                length_includes_head=True,
                zorder=3,
            )
            ax.add_patch(arrow)
        else:
            ax.plot([x_start - 0.5, x_end + 0.5], [y_center, y_center],
                    color=color, linewidth=1.5, solid_capstyle="round", zorder=2)


# =====================================================================
#  Gradient fill helper
# =====================================================================
def _gradient_fill(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    base_color: str,
    y_max: float,
) -> None:
    """Fill area under *y* with a value-mapped single-colour heatmap.

    Each x-position's fill intensity is proportional to y(x)/y_max:
    peaks → dark base_color, valleys → nearly white.  Combined with a
    polygon clip so only the area under the curve is coloured.
    """
    from matplotlib.path import Path as MPath

    if len(x) < 2 or y_max <= 0:
        return

    rgba_base = np.array(to_rgba(base_color))

    # Polygon clip: (x0,0) → data curve → (xN,0)
    verts = np.column_stack([
        np.concatenate([[x[0]], x, [x[-1]]]),
        np.concatenate([[0.0], y, [0.0]]),
    ])

    # Value-mapped intensity: each column = y_value / y_max (uniform vertically)
    n_rows = 4  # thin; clip hides everything outside curve anyway
    y_norm = np.clip(y / y_max, 0.0, 1.0)
    intensity = np.tile(y_norm, (n_rows, 1))  # shape (n_rows, len(x))

    # Colormap: white → base_color (RGB blending, not just alpha)
    white = np.array([1.0, 1.0, 1.0, 0.45])
    dark = rgba_base.copy()
    dark[3] = 0.72
    cmap = LinearSegmentedColormap.from_list("heatfill", [white, dark], N=256)

    im = ax.imshow(
        intensity, aspect="auto", origin="lower",
        extent=[x.min() - 0.5, x.max() + 0.5, 0, y_max],
        cmap=cmap, vmin=0.0, vmax=1.0,
        zorder=0, interpolation="bilinear",
    )
    im.set_clip_path(plt.Polygon(verts, transform=ax.transData))


# =====================================================================
#  Multi-chain residue panel helper
# =====================================================================
def _plot_multi_chain_panel(
    ax: plt.Axes,
    df: pd.DataFrame,
    y_col: str,
    y_label: str,
    ma_window: int = 7,
    fill: bool = False,
    y_cap: Optional[float] = None,
) -> None:
    """Plot protein + peptide chains on one axis with role-based colors.

    Protein → PROTEIN_COLOR (green), Peptide → PEPTIDE_COLOR (orange),
    with gradient fill (darker = higher values), dashed boundary line,
    and chain labels at the top.

    *y_cap*: if given, caps the gradient-fill reference maximum so that
    IDR peaks don't wash out interior flexibility differences.
    """
    chain_boundaries: List[Dict[str, Any]] = []
    y_global_max = float(df[y_col].max()) * 1.15 if not df.empty else 1.0
    if y_cap is not None:
        y_global_max = min(y_global_max, y_cap)

    for role in ["protein", "peptide"]:
        chunk = df[df["entity_role"].str.lower() == role]
        if chunk.empty:
            continue
        color = PROTEIN_COLOR if role == "protein" else PEPTIDE_COLOR
        label_text = role.capitalize()

        x_arr = chunk["plot_x"].values.astype(float)
        y_arr = chunk[y_col].values.astype(float)

        # Gradient fill (replaces flat fill_between)
        _gradient_fill(ax, x_arr, y_arr, color, y_global_max)

        # Raw line (light) + smoothed line (bold)
        ax.plot(x_arr, y_arr, color=color, linewidth=0.6, alpha=0.22)
        ax.plot(x_arr, _ma(y_arr, ma_window),
                color=color, linewidth=2.2, label=label_text)

        chain_boundaries.append({
            "role": role,
            "start": int(chunk["plot_x"].iloc[0]),
            "end": int(chunk["plot_x"].iloc[-1]),
            "chain_label": str(chunk["chain_label"].iloc[0]),
        })

    # Dashed boundary lines between chains
    for boundary in chain_boundaries[:-1]:
        bx = boundary["end"] + 0.5
        ax.axvline(bx, color="#C7D6E5", linewidth=0.8, linestyle="--", alpha=0.8)

    ax.set_ylabel(y_label, fontsize=8)
    ax.set_xlim(-1, len(df))

    # Chain labels at top (uses blended transform: data-x, axes-y)
    blend = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
    for boundary in chain_boundaries:
        mid = (boundary["start"] + boundary["end"]) / 2.0
        ax.text(mid, 1.02, f"Chain {boundary['chain_label']}",
                transform=blend, ha="center", va="bottom",
                fontsize=7, color="#425466", clip_on=False)

    # Compact legend when multi-chain
    if len(chain_boundaries) > 1:
        ax.legend(loc="upper right", frameon=False, fontsize=7, ncol=2)


# =====================================================================
#  Panel plotter
# =====================================================================
def _plot_system_panels(
    system_label: str,
    version: str,            # "v1" (DSSP) or "v2" (minimal)
    out_svg: Path,
    out_png: Path,
) -> None:
    """Generate a 2-row × 3-col panel figure for one system."""
    color = _color(system_label)
    display = _display(system_label)
    with_dssp = (version == "v1")

    if with_dssp:
        fig = plt.figure(figsize=(20.0, 12.0), dpi=600)
        fig.patch.set_facecolor("white")
        outer = gridspec.GridSpec(2, 3, hspace=0.28, wspace=0.20,
                                  left=0.05, right=0.97, top=0.86, bottom=0.07)
    else:
        fig, axes_flat = plt.subplots(2, 3, figsize=(20.0, 10.5), dpi=600,
                                       constrained_layout=True)
        fig.patch.set_facecolor("white")

    # Panel letters only on top row; R1/R2 label only on RMSD column
    panel_letters_top = ["A", "B", "C"]
    replica_labels = ["R1", "R2"]

    for row_idx, rep_label in enumerate(replica_labels):
        rmsd = _load_rmsd(system_label, rep_label)
        rmsf = _load_rmsf(system_label, rep_label)
        sasa = _load_sasa(system_label, rep_label)

        if with_dssp:
            dssp = _load_dssp(system_label, rep_label)
            # Build resid → plot_x mapping for DSSP alignment
            rmsf_prot = rmsf[rmsf["entity_role"].str.lower() == "protein"]
            rmsf_resid_map: Dict[int, int] = dict(
                zip(rmsf_prot["residue"].astype(int), rmsf_prot["plot_x"].astype(int))
            )
            sasa_prot = sasa[sasa["entity_role"].str.lower() == "protein"]
            sasa_resid_map: Dict[int, int] = dict(
                zip(sasa_prot["resSeq"].astype(int), sasa_prot["plot_x"].astype(int))
            )

        for col_idx in range(3):
            if with_dssp:
                if col_idx == 0:
                    inner = gridspec.GridSpecFromSubplotSpec(
                        1, 1, subplot_spec=outer[row_idx, col_idx])
                    ax_main = fig.add_subplot(inner[0])
                    ax_dssp = None
                else:
                    inner = gridspec.GridSpecFromSubplotSpec(
                        2, 1, subplot_spec=outer[row_idx, col_idx],
                        height_ratios=[5, 1], hspace=0.05)
                    ax_main = fig.add_subplot(inner[0])
                    ax_dssp = fig.add_subplot(inner[1], sharex=ax_main)
            else:
                ax_main = axes_flat[row_idx, col_idx]
                ax_dssp = None

            _style_axis(ax_main)

            if col_idx == 0:
                # ── RMSD (protein green + peptide orange) ──
                ax_main.plot(rmsd["time_ns"], rmsd["rmsd_A"],
                             color=PROTEIN_COLOR, linewidth=0.6, alpha=0.22)
                ax_main.plot(rmsd["time_ns"], _ma(rmsd["rmsd_A"], 41),
                             color=PROTEIN_COLOR, linewidth=2.2, label="Protein")
                rmsd_max = float(rmsd["rmsd_A"].max())
                # Peptide RMSD in orange
                if "rmsd_peptide_A" in rmsd.columns:
                    pep = rmsd.dropna(subset=["rmsd_peptide_A"])
                    if not pep.empty:
                        ax_main.plot(pep["time_ns"], pep["rmsd_peptide_A"],
                                     color=PEPTIDE_COLOR, linewidth=0.6, alpha=0.22)
                        ax_main.plot(pep["time_ns"], _ma(pep["rmsd_peptide_A"], 41),
                                     color=PEPTIDE_COLOR, linewidth=2.2, label="Peptide")
                        rmsd_max = max(rmsd_max, float(pep["rmsd_peptide_A"].max()))
                        ax_main.legend(loc="upper right", frameon=False, fontsize=7, ncol=2)
                ax_main.set_xlabel("Time (ns)", fontsize=8)
                ax_main.set_ylabel("Cα RMSD (Å)", fontsize=8)
                ax_main.set_xlim(0, TIME_CUTOFF_NS)
                ax_main.set_ylim(0, rmsd_max * 1.15)

            elif col_idx == 1:
                # ── RMSF (multi-chain: protein green + peptide orange) ──
                rmsf_cap = RMSF_GRADIENT_CAP.get(system_label)
                _plot_multi_chain_panel(
                    ax_main, rmsf, y_col="rmsf_A",
                    y_label="RMSF (Å)", ma_window=7, fill=False,
                    y_cap=rmsf_cap,
                )
                rmsf_max = float(rmsf["rmsf_A"].max())
                ax_main.set_ylim(0, rmsf_max * 1.15)

                # Tick labels
                n_ticks = min(20, len(rmsf))
                step = max(1, len(rmsf) // n_ticks)
                tick_rows = rmsf.iloc[::step]
                tick_positions = tick_rows["plot_x"].values
                tick_labels = tick_rows["tick_label"].values

                if with_dssp:
                    plt.setp(ax_main.get_xticklabels(), visible=False)
                else:
                    ax_main.set_xticks(tick_positions)
                    ax_main.set_xticklabels(tick_labels, rotation=55, ha="right", fontsize=5.5)
                    ax_main.set_xlabel("Residue", fontsize=8)

                if with_dssp and ax_dssp is not None:
                    _draw_dssp_track(ax_dssp, dssp, resid_to_plotx=rmsf_resid_map,
                                     y_center=0, height=0.6)
                    ax_dssp.set_ylim(-0.6, 0.6)
                    ax_dssp.set_xticks(tick_positions)
                    ax_dssp.set_xticklabels(tick_labels, rotation=55, ha="right", fontsize=5.5)
                    ax_dssp.set_xlabel("Residue", fontsize=8)
                    ax_dssp.set_yticks([])
                    ax_dssp.spines["top"].set_visible(False)
                    ax_dssp.spines["right"].set_visible(False)
                    ax_dssp.spines["left"].set_visible(False)

            else:
                # ── SASA (multi-chain: protein green + peptide orange) ──
                _plot_multi_chain_panel(
                    ax_main, sasa, y_col="mean_sasa_nm2",
                    y_label="Mean SASA (nm²)", ma_window=5, fill=True,
                )
                sasa_max = float(sasa["mean_sasa_nm2"].max())
                ax_main.set_ylim(0, max(sasa_max * 1.15, 0.25))

                # Tick labels
                n_ticks = min(20, len(sasa))
                step = max(1, len(sasa) // n_ticks)
                tick_rows = sasa.iloc[::step]
                tick_positions = tick_rows["plot_x"].values
                tick_labels = tick_rows["tick_label"].values

                if with_dssp:
                    plt.setp(ax_main.get_xticklabels(), visible=False)
                else:
                    ax_main.set_xticks(tick_positions)
                    ax_main.set_xticklabels(tick_labels, rotation=55, ha="right", fontsize=5.5)
                    ax_main.set_xlabel("Residue", fontsize=8)

                if with_dssp and ax_dssp is not None:
                    _draw_dssp_track(ax_dssp, dssp, resid_to_plotx=sasa_resid_map,
                                     y_center=0, height=0.6)
                    ax_dssp.set_ylim(-0.6, 0.6)
                    ax_dssp.set_xticks(tick_positions)
                    ax_dssp.set_xticklabels(tick_labels, rotation=55, ha="right", fontsize=5.5)
                    ax_dssp.set_xlabel("Residue", fontsize=8)
                    ax_dssp.set_yticks([])
                    ax_dssp.spines["top"].set_visible(False)
                    ax_dssp.spines["right"].set_visible(False)
                    ax_dssp.spines["left"].set_visible(False)

            # Panel letter (top row only) + R1/R2 label (RMSD column only)
            if row_idx == 0:
                ax_main.text(
                    -0.06, 1.08, panel_letters_top[col_idx],
                    transform=ax_main.transAxes,
                    fontsize=14, fontweight="bold", va="top", ha="right",
                )
            if col_idx == 0:
                ax_main.text(
                    0.02, 0.96, rep_label,
                    transform=ax_main.transAxes,
                    fontsize=9, fontweight="bold", va="top", ha="left",
                    color="#555555",
                )

    # Suptitle
    version_tag = "DSSP" if with_dssp else "Minimal"
    fig.suptitle(f"{display}  —  RMSD · RMSF · SASA  ({version_tag})",
                 fontsize=16, fontweight="bold", y=0.99 if with_dssp else 0.99)

    # DSSP legend for V1
    if with_dssp:
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=DSSP_PALETTE["helix"], label="α-helix"),
            Patch(facecolor=DSSP_PALETTE["strand"], label="β-strand"),
            Patch(facecolor=DSSP_PALETTE["coil"], label="Coil"),
        ]
        fig.legend(handles=legend_elements, loc="upper right",
                   frameon=True, fontsize=8, ncol=3,
                   bbox_to_anchor=(0.97, 0.97))

    out_svg.parent.mkdir(parents=True, exist_ok=True)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=600, facecolor="white", bbox_inches="tight")
    fig.savefig(out_svg, format="svg", facecolor="white", bbox_inches="tight")
    plt.close(fig)


# =====================================================================
#  All-systems overlay panel (1 replica per system, protein chain only)
# =====================================================================
def _plot_all_systems_overlay(out_dir: Path) -> Path:
    """Overlay all systems on 3 panels: RMSD · RMSF · SASA.

    One line per system (R1, protein chain only).
    * RMSD  — x = actual time (0–200 ns)
    * RMSF  — x = UniProt residue number per system
    * SASA  — x = UniProt residue number per system
    Systems are colour-coded by SYSTEM_COLORS.  Each shows a thin raw
    line (α=0.20) + bold smoothed line, matching the per-system panel style.
    """
    fig, axes = plt.subplots(1, 3, figsize=(22, 5.5), dpi=600)
    fig.patch.set_facecolor("white")
    ax_rmsd, ax_rmsf, ax_sasa = axes

    for sys_label in SYSTEM_ORDER:
        color = _color(sys_label)
        display_name = _display(sys_label)
        rep = "R1"

        # ── RMSD ──────────────────────────────────────────────────
        try:
            rmsd = _load_rmsd(sys_label, rep)
            t = rmsd["time_ns"].values
            y_r = rmsd["rmsd_A"].values
            ax_rmsd.plot(t, y_r, color=color, linewidth=0.5, alpha=0.20)
            ax_rmsd.plot(t, _ma(y_r, 41), color=color, linewidth=2.2, label=display_name)
        except Exception as exc:
            print(f"  [overlay RMSD skip] {sys_label}: {exc}")

        # ── RMSF (protein, UniProt residue numbering) ─────────────
        try:
            rmsf = _load_rmsf(sys_label, rep)
            prot_rmsf = rmsf[rmsf["entity_role"].str.lower() == "protein"].copy()
            if not prot_rmsf.empty:
                x_rf = prot_rmsf["residue"].values.astype(float)
                y_rf = prot_rmsf["rmsf_A"].values.astype(float)
                ax_rmsf.plot(x_rf, y_rf, color=color, linewidth=0.5, alpha=0.20)
                ax_rmsf.plot(x_rf, _ma(y_rf, 7), color=color, linewidth=2.2, label=display_name)
        except Exception as exc:
            print(f"  [overlay RMSF skip] {sys_label}: {exc}")

        # ── SASA (protein, UniProt residue numbering) ─────────────
        try:
            sasa = _load_sasa(sys_label, rep)
            prot_sasa = sasa[sasa["entity_role"].str.lower() == "protein"].copy()
            if not prot_sasa.empty:
                x_sa = prot_sasa["resSeq"].values.astype(float)
                y_sa = prot_sasa["mean_sasa_nm2"].values.astype(float)
                ax_sasa.plot(x_sa, y_sa, color=color, linewidth=0.5, alpha=0.20)
                ax_sasa.plot(x_sa, _ma(y_sa, 5), color=color, linewidth=2.2, label=display_name)
        except Exception as exc:
            print(f"  [overlay SASA skip] {sys_label}: {exc}")

    # ── Style ──────────────────────────────────────────────────────
    panel_configs = [
        (ax_rmsd, "Cα RMSD (Å)",    "Time (ns)"),
        (ax_rmsf, "RMSF (Å)",        "Residue (UniProt)"),
        (ax_sasa, "Mean SASA (nm²)", "Residue (UniProt)"),
    ]
    for ax, ylabel, xlabel in panel_configs:
        _style_axis(ax)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.tick_params(labelsize=8)

    ax_rmsd.set_xlim(0, TIME_CUTOFF_NS)
    ax_rmsd.legend(loc="upper left", frameon=False, fontsize=8, ncol=1)

    fig.suptitle(
        "All Systems — Cα RMSD · RMSF · SASA  (R1, protein chain)",
        fontsize=14, fontweight="bold", y=1.02,
    )
    plt.tight_layout()

    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / "all_systems_overlay_v1.png"
    out_svg = out_dir / "all_systems_overlay_v1.svg"
    fig.savefig(out_png, dpi=600, facecolor="white", bbox_inches="tight")
    fig.savefig(out_svg, format="svg", facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ All-systems overlay → {out_png.name}")
    return out_png


# =====================================================================
#  DOCX builder
# =====================================================================
def _set_landscape(section) -> None:
    section.orientation = WD_ORIENTATION.LANDSCAPE
    section.page_width, section.page_height = section.page_height, section.page_width
    section.top_margin = Inches(0.3)
    section.bottom_margin = Inches(0.3)
    section.left_margin = Inches(0.35)
    section.right_margin = Inches(0.35)


def _build_docx(panel_files: List[Dict[str, Any]]) -> Path:
    doc = Document()
    _set_landscape(doc.sections[0])
    normal = doc.styles["Normal"]
    normal.font.name = "Arial"
    normal.font.size = Pt(9)

    # Title page
    title = doc.add_paragraph()
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = title.add_run("Control MD — Per-System Analysis Panels\nRMSD · RMSF · SASA")
    run.bold = True
    run.font.size = Pt(18)
    run.font.name = "Arial"

    for entry in panel_files:
        doc.add_page_break()
        system_label = entry["system_label"]
        display = _display(system_label)
        version_label = entry["version_label"]

        heading = doc.add_paragraph()
        heading.alignment = WD_ALIGN_PARAGRAPH.CENTER
        hrun = heading.add_run(f"{display}  —  {version_label}")
        hrun.bold = True
        hrun.font.size = Pt(14)
        hrun.font.name = "Arial"

        img_para = doc.add_paragraph()
        img_para.alignment = WD_ALIGN_PARAGRAPH.CENTER
        img_run = img_para.add_run()
        img_run.add_picture(str(entry["png_path"]), width=Inches(9.5))

        caption = doc.add_paragraph()
        caption.alignment = WD_ALIGN_PARAGRAPH.CENTER
        cap_run = caption.add_run(
            f"Figure: {display} — Panels A: RMSD, B: RMSF, C: SASA. "
            f"Top row: R1, Bottom row: R2. {version_label}."
        )
        cap_run.font.size = Pt(8)
        cap_run.italic = True

    OUT_DOCX.parent.mkdir(parents=True, exist_ok=True)
    doc.save(str(OUT_DOCX))
    return OUT_DOCX


# =====================================================================
#  Main
# =====================================================================
def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    panel_records: List[Dict[str, Any]] = []

    # 1. Generate per-system panel figures
    print("[1/2] Generating per-system panel figures …")
    for sys_label in SYSTEM_ORDER:
        slug = _slugify(sys_label)
        display = _display(sys_label)
        for version, vtag in [("v1", "DSSP")]:
            out_png = OUT_DIR / f"{slug}_panels_{version}.png"
            out_svg = OUT_DIR / f"{slug}_panels_{version}.svg"
            _plot_system_panels(sys_label, version, out_svg, out_png)
            panel_records.append({
                "system_label": sys_label,
                "display_label": display,
                "version": version,
                "version_label": vtag,
                "png_path": out_png,
                "svg_path": out_svg,
            })
            print(f"  ✓ {display} {version} → PNG + SVG")

    # 2b. All-systems overlay
    print("[2b] Generating all-systems overlay …")
    _plot_all_systems_overlay(OUT_DIR)

    # 2. Build DOCX
    print("[2/2] Building DOCX …")
    docx_path = _build_docx(panel_records)

    print(f"\n[OK] All panel figures in: {OUT_DIR}")
    print(f"[OK] DOCX report: {docx_path}")


if __name__ == "__main__":
    main()

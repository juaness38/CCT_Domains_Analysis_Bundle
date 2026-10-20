# -*- coding: utf-8 -*-
"""
CCT-OSR1 Binding Analysis - RMSD/RMSF Analysis
===============================================
Structural deviation analysis with interactive Plotly visualizations.
"""

import numpy as np
import mdtraj as md
import pandas as pd
from pathlib import Path
from typing import Optional, Tuple, Dict, Union, List, Any


def _pdb_has_non_numeric_resseq(pdb_path: Path, *, max_lines: int = 250_000) -> bool:
    """Detect non-numeric residue sequence fields in PDB ATOM/HETATM/TER records.

    Some OpenMM-written PDBs overflow the 4-char residue number field and emit
    base36-ish strings (e.g. 'F0C2'), which can crash strict parsers.
    """
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for idx, line in enumerate(handle, start=1):
            if idx > max_lines:
                break
            if not (line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER")):
                continue
            if len(line) < 26:
                continue
            raw = line[22:26]
            raw_stripped = raw.strip()
            if not raw_stripped:
                continue
            if not raw_stripped.isdigit():
                return True
    return False


def _sanitize_pdb_resseq(input_pdb: Path, output_pdb: Path) -> None:
    """Rewrite PDB with numeric residue sequence numbers while preserving atom order."""
    output_pdb.parent.mkdir(parents=True, exist_ok=True)

    current_key = None
    new_resseq = 0
    last_written_resseq = 1

    with open(input_pdb, "r", encoding="utf-8", errors="replace") as src, open(
        output_pdb, "w", encoding="utf-8", newline="\n"
    ) as dst:
        for line in src:
            # Avoid mdtraj triggering NumPy linalg det() on CRYST1/unitcell parsing
            # in some Windows environments.
            if line.startswith("CRYST1"):
                continue
            record = line[0:6]
            if record.startswith("ATOM") or record.startswith("HETATM"):
                if len(line) < 27:
                    dst.write(line)
                    continue

                chain_id = line[21:22]
                resname = line[17:20]
                resseq_raw = line[22:26]
                icode = line[26:27]
                key = (chain_id, resname, resseq_raw, icode)

                if key != current_key:
                    current_key = key
                    new_resseq += 1
                    if new_resseq > 9999:
                        new_resseq = 1
                    last_written_resseq = new_resseq

                fixed = f"{last_written_resseq:4d}"
                dst.write(line[:22] + fixed + line[26:])
                continue

            if record.startswith("TER") and len(line) >= 27:
                fixed = f"{last_written_resseq:4d}"
                dst.write(line[:22] + fixed + line[26:])
                continue

            dst.write(line)


def _prepare_topology_for_mdtraj(topology_path: Path, prepared_dir: Path) -> Path:
    """Return a mdtraj-safe topology PDB path.

    On some Windows builds, mdtraj's PDB parsing can trigger NumPy linalg
    operations (e.g. determinant of unitcell vectors) that hard-crash the
    interpreter. We defensively write a prepared PDB that:
    - strips CRYST1 records (unit cell)
    - sanitizes non-numeric resSeq fields

    The periodic box can still be supplied via OpenMM State XML (box_lengths_nm).
    """
    topology_path = Path(topology_path)
    if topology_path.suffix.lower() != ".pdb" or not topology_path.exists():
        return topology_path

    needs_sanitize = _pdb_has_non_numeric_resseq(topology_path)
    has_cryst1 = False
    try:
        with open(topology_path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("CRYST1"):
                    has_cryst1 = True
                    break
    except Exception:
        # If we can't inspect, don't rewrite.
        return topology_path

    if not needs_sanitize and not has_cryst1:
        return topology_path

    prepared_dir.mkdir(parents=True, exist_ok=True)
    suffix = []
    if needs_sanitize:
        suffix.append("sanres")
    if has_cryst1:
        suffix.append("nocryst")
    tag = "_".join(suffix) if suffix else "prepared"
    prepared_path = prepared_dir / f"{topology_path.stem}__{tag}.pdb"

    if needs_sanitize:
        _sanitize_pdb_resseq(topology_path, prepared_path)
        return prepared_path

    # Only strip CRYST1 (and leave residue numbering as-is)
    with open(topology_path, "r", encoding="utf-8", errors="replace") as src, open(
        prepared_path, "w", encoding="utf-8", newline="\n"
    ) as dst:
        for line in src:
            if line.startswith("CRYST1"):
                continue
            dst.write(line)
    return prepared_path

try:  # Package-first configuration
    from .config import (
        SAVE_FREQ_PS, PLOTLY_TEMPLATE, COLOR_SCHEME,
        FIG_WIDTH, FIG_HEIGHT, FIG_HEIGHT_TALL,
        HTML_DIR, FIGURES_DIR, DATA_OUTPUT_DIR, RFXV_MOTIF,
        TOPOLOGY_PDB, TRAJECTORIES, OUTPUT_DIR,
    )
except Exception:  # pragma: no cover
    try:  # Legacy configuration (optional)
        from config import (
            SAVE_FREQ_PS, PLOTLY_TEMPLATE, COLOR_SCHEME,
            FIG_WIDTH, FIG_HEIGHT, FIG_HEIGHT_TALL,
            HTML_DIR, FIGURES_DIR, DATA_OUTPUT_DIR, RFXV_MOTIF,
            TOPOLOGY_PDB, TRAJECTORIES, OUTPUT_DIR,
        )
    except ImportError:  # pragma: no cover - fallback defaults
        SAVE_FREQ_PS = None
        PLOTLY_TEMPLATE = 'plotly_white'
        COLOR_SCHEME = {
            'protein': '#1f77b4',
            'peptide': '#ff7f0e',
            'rfxv': '#d62728',
            'water': '#17becf',
            'contact': '#2ca02c',
            'highlight': '#9467bd',
        }
        FIG_WIDTH = 1200
        FIG_HEIGHT = 600
        FIG_HEIGHT_TALL = 800
        HTML_DIR = Path("outputs/html")
        FIGURES_DIR = Path("outputs/figures")
        DATA_OUTPUT_DIR = Path("outputs/data")
        OUTPUT_DIR = Path("outputs")
        TOPOLOGY_PDB = None
        TRAJECTORIES = []
    RFXV_MOTIF = {'residues': [], 'core': [], 'sequence': ''}

try:
    from .runtime_config import AnalysisConfig, build_cli_parser, resolve_config
    from .trajectory_loader import (
        load_and_concatenate_trajectories,
        extract_protein_peptide,
    )
except ImportError:  # pragma: no cover
    from runtime_config import AnalysisConfig, build_cli_parser, resolve_config
    from trajectory_loader import (
        load_and_concatenate_trajectories,
        extract_protein_peptide,
    )

RMSD_COLORS = {
    'protein_backbone': '#e63946',      # red - backbone (N,Cα,C,O)
    'protein_ca': '#a8dadc',            # light blue - Cα only
    'peptide_backbone': '#e63946',      # red - backbone
    'peptide_all': '#1d3557',           # dark blue - all atoms
    'rfxv': '#9d4edd',                  # purple - RFxV motif (all atoms)
    'rfxv_sidechain': '#ff6b00',        # ORANGE - RFxV sidechains only (DESTACADO)
    'diagnostics': '#457b9d',           # slate - diagnostics
    'motion': '#f4a261',                # amber - motion metrics
}


def _hex_to_rgba(hex_color: str, alpha: float) -> str:
    hex_color = hex_color.lstrip('#')
    if len(hex_color) != 6:
        raise ValueError("Hex color must be 6 characters")
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


def _resolve_time_axis(
    traj: md.Trajectory,
    config: Optional[AnalysisConfig] = None
) -> Tuple[np.ndarray, float]:
    """Return per-frame times (ps) and timestep, falling back to config/env."""

    def _try_time_end_ps_from_state_xml() -> Optional[float]:
        if config is None:
            return None
        xml_path = getattr(config, "box_source", None)
        if not xml_path:
            return None
        try:
            from xml.etree import ElementTree as ET

            root = ET.parse(str(xml_path)).getroot()
            val = root.attrib.get("time")
            if val is None:
                return None
            time_end = float(val)
            return time_end if (np.isfinite(time_end) and time_end > 0) else None
        except Exception:
            return None

    time_attr = getattr(traj, "time", None)
    time_ps = None
    frame_dt = None

    if time_attr is not None and len(time_attr) == traj.n_frames:
        time_ps = np.asarray(time_attr, dtype=float)
        diffs = np.diff(time_ps) if traj.n_frames > 1 else np.array([])
        if diffs.size and np.any(diffs > 0):
            frame_dt = float(np.median(diffs[diffs > 0]))

    # If we have an OpenMM State XML with an authoritative end time, prefer it
    # over guessed/default time axes (common on Windows when DCD/XTC lack time).
    time_end_ps = _try_time_end_ps_from_state_xml()
    if time_end_ps is not None and traj.n_frames > 1:
        use_xml = False
        if time_ps is None or len(time_ps) != traj.n_frames:
            use_xml = True
        elif not np.isfinite(time_ps[-1]) or time_ps[-1] <= 0:
            use_xml = True
        else:
            rel_err = abs(float(time_ps[-1]) - float(time_end_ps)) / float(time_end_ps)
            if rel_err > 0.10:
                use_xml = True

        if use_xml:
            time_ps = np.linspace(0.0, float(time_end_ps), num=traj.n_frames, dtype=float)
            traj.time = time_ps
            diffs = np.diff(time_ps)
            frame_dt = float(np.median(diffs[diffs > 0])) if diffs.size else 0.0

    if frame_dt is None:
        frame_dt = (config.time_step_ps if config and config.time_step_ps else SAVE_FREQ_PS)
    if frame_dt is None or frame_dt <= 0:
        frame_dt = 1.0  # conservative fallback

    if time_ps is None or len(time_ps) != traj.n_frames:
        time_ps = np.arange(traj.n_frames, dtype=float) * frame_dt
        traj.time = time_ps

    return time_ps, frame_dt


def _copy_trajectory(traj: md.Trajectory) -> md.Trajectory:
    """Return an explicit trajectory copy to avoid in-place mutation side effects."""
    return traj.slice(np.arange(traj.n_frames), copy=True)


def compute_rmsd_timeseries(
    traj: md.Trajectory,
    reference_frame: int = 0,
    atom_indices: Optional[np.ndarray] = None,
    superpose: bool = True
) -> np.ndarray:
    """
    Compute RMSD vs time relative to a reference frame.
    
    Parameters
    ----------
    traj : md.Trajectory
    reference_frame : int
        Frame to use as reference (default: first frame)
    atom_indices : array-like, optional
        Atoms to include in RMSD calculation
    superpose : bool
        Whether to superpose structures before RMSD calculation
        
    Returns
    -------
    np.ndarray
        RMSD values in nm for each frame
    """
    traj_work = _copy_trajectory(traj)

    if superpose and atom_indices is not None:
        traj_aligned = traj_work.superpose(traj_work, frame=reference_frame, atom_indices=atom_indices)
    elif superpose:
        traj_aligned = traj_work.superpose(traj_work, frame=reference_frame)
    else:
        traj_aligned = traj_work
    
    if atom_indices is not None:
        rmsd = md.rmsd(traj_aligned, traj_aligned, frame=reference_frame, atom_indices=atom_indices)
    else:
        rmsd = md.rmsd(traj_aligned, traj_aligned, frame=reference_frame)
    
    return rmsd


def compute_rmsf(
    traj: md.Trajectory,
    atom_indices: Optional[np.ndarray] = None,
    superpose_indices: Optional[np.ndarray] = None
) -> np.ndarray:
    """
    Compute Root Mean Square Fluctuation per atom.
    
    Parameters
    ----------
    traj : md.Trajectory
    atom_indices : array-like, optional
        Atoms for which to compute RMSF
    superpose_indices : array-like, optional
        Atoms to use for superposition (default: same as atom_indices)
        
    Returns
    -------
    np.ndarray
        RMSF values in nm for each atom
    """
    if superpose_indices is None:
        superpose_indices = atom_indices
    
    # Superpose
    traj_work = _copy_trajectory(traj)

    if superpose_indices is not None:
        traj_aligned = traj_work.superpose(traj_work, frame=0, atom_indices=superpose_indices)
    else:
        traj_aligned = traj_work.superpose(traj_work, frame=0)
    
    # Get positions
    if atom_indices is not None:
        xyz = traj_aligned.xyz[:, atom_indices, :]
    else:
        xyz = traj_aligned.xyz
    
    # Compute mean position
    mean_xyz = xyz.mean(axis=0)
    
    # RMSF = sqrt(mean((x - <x>)^2))
    deviations = xyz - mean_xyz
    rmsf = np.sqrt(np.mean(deviations**2, axis=0).sum(axis=1))
    
    return rmsf


def _moving_average(values: np.ndarray, window: int = 25) -> np.ndarray:
    """Centered moving average; falls back gracefully for short series."""
    if values.size == 0 or window <= 1:
        return values
    series = pd.Series(values)
    return series.rolling(window=window, center=True, min_periods=max(1, window // 2)).mean().to_numpy()


def _basic_stats(values: np.ndarray) -> Dict[str, float]:
    """Return mean, median, std, and p95 (Å units expected downstream)."""
    if values.size == 0:
        return {"mean": float("nan"), "median": float("nan"), "std": float("nan"), "p95": float("nan")}
    return {
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "std": float(np.std(values, ddof=1)) if values.size > 1 else 0.0,
        "p95": float(np.percentile(values, 95)),
    }


def _safe_write_image(fig: Any, output_png: Path, *, scale: int = 2) -> None:
    """Write static image, falling back gracefully when Kaleido is absent."""
    try:
        fig.write_image(str(output_png), scale=scale)
        print(f"Saved: {output_png}")
    except ValueError as exc:
        if "kaleido" in str(exc).lower():
            print(f"Skipping PNG export (install 'kaleido' to enable): {output_png}")
        else:
            raise


def _resolve_chain_roles(topology: md.Topology) -> Tuple[Optional[int], Optional[int]]:
    """Infer protein/peptide-like chains by atom count.

    Returns (protein_chainid, peptide_chainid). If only one chain exists,
    peptide_chainid is None.
    """
    chain_sizes: List[Tuple[int, int]] = []
    for chain in topology.chains:
        atom_count = sum(1 for _ in chain.atoms)
        chain_sizes.append((chain.index, atom_count))

    if not chain_sizes:
        return None, None

    chain_sizes.sort(key=lambda x: x[1], reverse=True)
    protein_chain = chain_sizes[0][0]
    peptide_chain = chain_sizes[1][0] if len(chain_sizes) > 1 else None
    return protein_chain, peptide_chain


def compute_all_rmsd(traj: md.Trajectory) -> Dict[str, np.ndarray]:
    """Compute RMSD for multiple selections."""
    top = traj.topology
    protein_chain, peptide_chain = _resolve_chain_roles(top)
    if protein_chain is None:
        return {}
    
    results = {}
    
    # Protein backbone
    protein_bb = top.select(f'chainid {protein_chain} and backbone')
    if len(protein_bb) > 0:
        results['protein_backbone'] = compute_rmsd_timeseries(traj, atom_indices=protein_bb)
    
    # Protein CA
    protein_ca = top.select(f'chainid {protein_chain} and name CA')
    if len(protein_ca) > 0:
        results['protein_ca'] = compute_rmsd_timeseries(traj, atom_indices=protein_ca)
    
    # Protein-aligned working copy for peptide-relative metrics
    traj_aligned_on_protein: Optional[md.Trajectory] = None
    if len(protein_ca) > 0:
        traj_aligned_on_protein = _copy_trajectory(traj)
        traj_aligned_on_protein.superpose(traj_aligned_on_protein, frame=0, atom_indices=protein_ca)

    # Peptide backbone (superposed on protein)
    peptide_bb = top.select(f'chainid {peptide_chain} and backbone') if peptide_chain is not None else np.array([])
    if len(peptide_bb) > 0 and traj_aligned_on_protein is not None:
        results['peptide_backbone'] = md.rmsd(
            traj_aligned_on_protein,
            traj_aligned_on_protein,
            frame=0,
            atom_indices=peptide_bb,
        )
    
    # Peptide all atoms
    peptide_all = top.select(f'chainid {peptide_chain}') if peptide_chain is not None else np.array([])
    if len(peptide_all) > 0 and traj_aligned_on_protein is not None:
        results['peptide_all'] = md.rmsd(
            traj_aligned_on_protein,
            traj_aligned_on_protein,
            frame=0,
            atom_indices=peptide_all,
        )
    
    return results


def compute_all_rmsf(traj: md.Trajectory) -> Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Compute per-residue RMSF for three atom sets: Cα, all heavy atoms, sidechains.

    Returns
    -------
    Dict mapping key → (resSeq_arr, resname_arr, rmsf_nm_arr).

    Keys produced:
      protein_ca, protein_all, protein_sidechain
      peptide_ca, peptide_all, peptide_sidechain

    Backwards-compat aliases (same object, no extra data):
      protein  → protein_ca
      peptide, peptide_relative → peptide_ca
    """
    from collections import OrderedDict

    top = traj.topology
    protein_chain, peptide_chain = _resolve_chain_roles(top)
    if protein_chain is None:
        return {}

    protein_ca_idx = top.select(f'chainid {protein_chain} and name CA')

    # Backbone heavy-atom names excluded from sidechain RMSF.
    _BB_NAMES: frozenset = frozenset({'N', 'CA', 'C', 'O', 'OXT'})

    def _per_residue(
        chain_id: int,
        ca_only: bool,
        sidechain_only: bool,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (resSeq_arr, resname_arr, rmsf_nm_arr) for one chain+atom-set."""
        if ca_only:
            sel = top.select(f'chainid {chain_id} and name CA')
        else:
            heavy = top.select(f'chainid {chain_id} and not element H')
            if sidechain_only:
                sel = np.array(
                    [i for i in heavy if top.atom(int(i)).name not in _BB_NAMES],
                    dtype=int,
                )
            else:
                sel = heavy

        if len(sel) == 0:
            return (np.array([], dtype=int),
                    np.array([], dtype='U8'),
                    np.array([], dtype=float))

        # Align trajectory on protein Cα
        sp_idx = protein_ca_idx if len(protein_ca_idx) > 0 else sel
        traj_work = _copy_trajectory(traj)
        traj_aligned = traj_work.superpose(traj_work, frame=0, atom_indices=sp_idx)

        xyz = traj_aligned.xyz[:, sel, :]            # (n_frames, n_atoms, 3)
        mean_xyz = xyz.mean(axis=0)                  # (n_atoms, 3)
        diff = xyz - mean_xyz
        # per-atom RMSF = sqrt( mean_t( ||r_i(t) - <r_i>||^2 ) )
        per_atom_rmsf = np.sqrt(
            np.mean(np.sum(diff ** 2, axis=2), axis=0)
        )  # (n_atoms,)

        # Aggregate by residue (preserve insertion order)
        residue_map: 'OrderedDict[int, list]' = OrderedDict()
        for i, atom_idx in enumerate(sel):
            res = top.atom(int(atom_idx)).residue
            key = res.index
            if key not in residue_map:
                residue_map[key] = [res.resSeq, res.name, []]
            residue_map[key][2].append(float(per_atom_rmsf[i]))

        resSeqs, resnames, rmsfs = [], [], []
        for _, (rsq, rnm, vals) in residue_map.items():
            resSeqs.append(rsq)
            resnames.append(rnm)
            rmsfs.append(float(np.mean(vals)))

        return (
            np.array(resSeqs, dtype=int),
            np.array(resnames, dtype='U8'),
            np.array(rmsfs, dtype=float),
        )

    results: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    # --- Protein chain ---
    if len(protein_ca_idx) > 0:
        results['protein_ca'] = _per_residue(protein_chain, ca_only=True, sidechain_only=False)
        results['protein_all'] = _per_residue(protein_chain, ca_only=False, sidechain_only=False)
        results['protein_sidechain'] = _per_residue(protein_chain, ca_only=False, sidechain_only=True)
        results['protein'] = results['protein_ca']  # backwards compat

    # --- Peptide chain ---
    if peptide_chain is not None:
        results['peptide_ca'] = _per_residue(peptide_chain, ca_only=True, sidechain_only=False)
        results['peptide_all'] = _per_residue(peptide_chain, ca_only=False, sidechain_only=False)
        results['peptide_sidechain'] = _per_residue(peptide_chain, ca_only=False, sidechain_only=True)
        results['peptide'] = results['peptide_ca']           # backwards compat
        results['peptide_relative'] = results['peptide_ca']  # backwards compat

    return results


def plot_rmsd_timeseries(
    rmsd_data: Dict[str, np.ndarray],
    time_ps: np.ndarray,
    *,
    label: str,
    output_html: Path = None,
    output_png: Path = None,
    smooth_window: int = 25,
) -> Any:
    """
    Create interactive RMSD vs time plot with 3 panels.
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    time_ns = time_ps / 1000
    
    # Calculate stats for title
    protein_stats = _basic_stats(rmsd_data.get('protein_ca', np.array([0])) * 10)
    peptide_stats = _basic_stats(rmsd_data.get('peptide_all', np.array([0])) * 10)

    # Always 2 panels: protein stability + peptide global structure
    has_core = False
    n_rows = 2

    subplot_titles = [
        "Protein Stability (Cα + Backbone)",
        "Peptide Global Structure (Backbone vs All Atoms)",
    ]
    
    fig = make_subplots(
        rows=n_rows,
        cols=1,
        subplot_titles=tuple(subplot_titles),
        vertical_spacing=0.08,
        shared_xaxes=True,
    )

    # Series specification with panel assignments
    series_spec: List[Dict[str, Union[str, int]]] = [
        # Panel 1: Protein
        {'key': 'protein_backbone', 'row': 1, 'name': 'Backbone', 'color': RMSD_COLORS['protein_backbone'], 'width': 1.8},
        {'key': 'protein_ca', 'row': 1, 'name': 'Cα only', 'color': RMSD_COLORS['protein_ca'], 'width': 2.2},
        # Panel 2: Peptide global
        {'key': 'peptide_backbone', 'row': 2, 'name': 'Backbone', 'color': RMSD_COLORS['peptide_backbone'], 'width': 1.8},
        {'key': 'peptide_all', 'row': 2, 'name': 'All Atoms', 'color': RMSD_COLORS['peptide_all'], 'width': 2.4},
    ]
    
    # (no 3rd panel — core/RFxV series removed)

    legend_groups: Dict[str, str] = {}
    for spec in series_spec:
        key = spec['key']
        if key not in rmsd_data:
            continue
        values_A = rmsd_data[key] * 10
        color = spec['color']
        legend_group = f"panel{spec['row']}"
        legend_groups.setdefault(legend_group, legend_group)

        fig.add_trace(
            go.Scatter(
                x=time_ns,
                y=values_A,
                name=spec['name'],
                line=dict(color=color, width=spec.get('width', 2), dash=spec.get('dash')),
                legendgroup=legend_group,
                hovertemplate="t = %{x:.1f} ns<br>RMSD = %{y:.2f} Å<extra></extra>",
            ),
            row=spec['row'],
            col=1,
        )

        # Rolling average (subtle)
        smoothed = _moving_average(values_A, window=smooth_window)
        fig.add_trace(
            go.Scatter(
                x=time_ns,
                y=smoothed,
                name=f"{spec['name']} · avg",
                line=dict(color=color, width=3, dash='solid'),
                legendgroup=legend_group,
                showlegend=False,
                opacity=0.3,
                hoverinfo='skip',
            ),
            row=spec['row'],
            col=1,
        )

    # === STATISTICAL SHADING FOR EACH PANEL ===
    
    # Panel 1: Protein stats
    if 'protein_ca' in rmsd_data:
        mu, sigma = protein_stats['mean'], protein_stats['std']
        if not np.isnan(mu):
            fig.add_hrect(y0=mu - sigma, y1=mu + sigma, line_width=0,
                         fillcolor=_hex_to_rgba(RMSD_COLORS['protein_ca'], 0.12), row=1, col=1)
            fig.add_hline(y=mu, line_dash='dash', line_color=RMSD_COLORS['protein_ca'],
                         row=1, col=1, annotation_text=f"μ={mu:.2f}Å", annotation_position='top left')

    # Panel 2: Peptide stats
    if 'peptide_all' in rmsd_data:
        mu, sigma = peptide_stats['mean'], peptide_stats['std']
        if not np.isnan(mu):
            fig.add_hrect(y0=mu - sigma, y1=mu + sigma, line_width=0,
                         fillcolor=_hex_to_rgba(RMSD_COLORS['peptide_all'], 0.10), row=2, col=1)
            fig.add_hline(y=mu, line_dash='dash', line_color=RMSD_COLORS['peptide_all'],
                         row=2, col=1, annotation_text=f"μ={mu:.2f}Å", annotation_position='top left')

    # === PEAK ANNOTATIONS ===
    # Peptide peak detection
    if 'peptide_all' in rmsd_data:
        pep_vals = rmsd_data['peptide_all'] * 10
        if np.max(pep_vals) > 3.0:
            peak_idx = np.argmax(pep_vals)
            peak_time = time_ns[peak_idx]
            peak_val = pep_vals[peak_idx]
            fig.add_annotation(
                x=peak_time, y=peak_val,
                text=f"Peak: {peak_val:.1f} Å @ {peak_time:.0f}ns",
                showarrow=True, arrowhead=2, arrowcolor='#E63946',
                font=dict(size=10, color='#E63946'),
                row=2, col=1
            )
    
    # === AXIS FORMATTING ===
    # Only show x-axis label on bottom panel
    fig.update_xaxes(title_text="Time (ns)", row=n_rows, col=1)
    
    # Y-axis labels and gridlines for all panels
    for row in range(1, n_rows + 1):
        fig.update_yaxes(
            title_text="RMSD (Å)", 
            row=row, col=1,
            gridcolor='rgba(128, 128, 128, 0.15)',
            gridwidth=1,
            zeroline=False,
        )
        fig.update_xaxes(
            row=row, col=1,
            gridcolor='rgba(128, 128, 128, 0.15)',
            gridwidth=1,
        )

    # === TITLE WITH STATS ===
    title_text = (
        f"<b>{label.upper()}</b> · RMSD Analysis<br>"
        f"<sub>Simulation: {time_ns[-1]:.0f} ns | "
        f"Protein μ={protein_stats['mean']:.2f}±{protein_stats['std']:.2f} Å | "
        f"Peptide μ={peptide_stats['mean']:.2f}±{peptide_stats['std']:.2f} Å</sub>"
    )

    # === LAYOUT ===
    fig.update_layout(
        title=dict(text=title_text, x=0.5, xanchor='center'),
        template=PLOTLY_TEMPLATE,
        width=FIG_WIDTH,
        height=280 * n_rows + 100,  # Dynamic height
        legend=dict(
            orientation="v",
            yanchor="middle",
            y=0.5,
            xanchor="left",
            x=1.02,
            bgcolor="rgba(255,255,255,0.9)",
            bordercolor="rgba(0,0,0,0.2)",
            borderwidth=1,
            font=dict(size=11),
        ),
        margin=dict(t=90, l=60, r=120, b=60),
        plot_bgcolor='white',
    )

    # Save
    if output_html:
        output_html.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(str(output_html))
        print(f"Saved: {output_html}")
    
    if output_png:
        output_png.parent.mkdir(parents=True, exist_ok=True)
        _safe_write_image(fig, output_png, scale=2)
    
    return fig


def plot_rmsf(
    rmsf_data: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]],
    *,
    label: str,
    output_html: Path = None,
    output_png: Path = None,
    topology: Optional[md.Topology] = None,
) -> Any:
    """
    Interactive RMSF plot: 3 rows × 2 columns.

    Rows  : Cα RMSF  |  All-heavy-atoms RMSF  |  Sidechain RMSF
    Cols  : Protein chain  |  Peptide chain

    Each bar label on the x-axis shows '<RESNAME><resSeq>' (e.g. 'GLY47').
    Every ~20th residue is ticked to avoid crowding.
    Hover shows the full per-residue AA label and RMSF value.
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    PROT_COLOR: str = COLOR_SCHEME.get('protein', '#1f77b4')
    PEP_COLOR: str  = COLOR_SCHEME.get('peptide', '#ff7f0e')

    def _aa_labels(resSeqs: np.ndarray, resnames: np.ndarray) -> List[str]:
        """Build 'AAA123' labels from resSeq + 3-letter residue name."""
        out: List[str] = []
        for rsq, rnm in zip(resSeqs.astype(int), resnames):
            nm = str(rnm).strip()
            out.append(f"{nm}{rsq}" if nm else str(rsq))
        return out

    def _add_panel(
        fig: Any,
        row: int,
        col: int,
        key: str,
        trace_name: str,
        color: str,
    ) -> None:
        """Add a bar trace + mean band/line to one subplot panel."""
        if key not in rmsf_data:
            return
        data_tuple = rmsf_data[key]
        # Accept both 3-tuple (new) and 2-tuple (legacy)
        if len(data_tuple) == 3:
            resSeqs, resnames, rmsf_nm = data_tuple
        else:
            resSeqs, rmsf_nm = data_tuple
            resnames = np.array([''] * len(resSeqs), dtype='U8')

        if len(resSeqs) == 0:
            return

        rmsf_A = (np.asarray(rmsf_nm, dtype=float)) * 10.0
        xlabels = _aa_labels(np.asarray(resSeqs, dtype=int), resnames)
        stats = _basic_stats(rmsf_A)

        # Tick every ~20 residues (30 ticks max)
        n = len(resSeqs)
        step = max(1, n // 20)
        tick_idx = list(range(0, n, step))

        fig.add_trace(
            go.Bar(
                x=list(resSeqs.astype(int)),
                y=list(rmsf_A),
                name=trace_name,
                marker_color=color,
                opacity=0.82,
                customdata=[[xlabels[i], round(float(rmsf_A[i]), 3)] for i in range(n)],
                hovertemplate=(
                    '<b>%{customdata[0]}</b><br>'
                    'RMSF = %{customdata[1]:.2f} Å'
                    '<extra></extra>'
                ),
            ),
            row=row, col=col,
        )

        fig.update_xaxes(
            tickmode='array',
            tickvals=[int(resSeqs[i]) for i in tick_idx],
            ticktext=[xlabels[i] for i in tick_idx],
            tickangle=45,
            row=row, col=col,
        )

        mu, sigma = stats['mean'], stats['std']
        if np.isfinite(mu):
            fig.add_hrect(
                y0=mu - sigma, y1=mu + sigma,
                line_width=0,
                fillcolor=_hex_to_rgba(color, 0.12),
                row=row, col=col,
            )
            fig.add_hline(
                y=mu,
                line_dash='dash',
                line_color=color,
                annotation_text=f'μ={mu:.2f} Å',
                annotation_position='top left',
                row=row, col=col,
            )

    has_peptide = (
        'peptide_ca' in rmsf_data
        and len(rmsf_data['peptide_ca'][0]) > 0
    )
    n_cols = 2 if has_peptide else 1

    chain_labels = ['Protein']
    if has_peptide:
        chain_labels.append('Peptide')

    rmsf_row_labels = ['Cα RMSF', 'Total RMSF (all heavy atoms)', 'Sidechain RMSF']
    subplot_titles: List[str] = []
    for rtl in rmsf_row_labels:
        for cl in chain_labels:
            subplot_titles.append(f'{rtl} — {cl}')

    fig = make_subplots(
        rows=3,
        cols=n_cols,
        subplot_titles=subplot_titles,
        vertical_spacing=0.12,
        horizontal_spacing=0.08,
    )

    # Row 1 — Cα
    _add_panel(fig, 1, 1, 'protein_ca',        'Protein Cα',          PROT_COLOR)
    if has_peptide:
        _add_panel(fig, 1, 2, 'peptide_ca',    'Peptide Cα',          PEP_COLOR)

    # Row 2 — All heavy atoms
    _add_panel(fig, 2, 1, 'protein_all',       'Protein (all heavy)', PROT_COLOR)
    if has_peptide:
        _add_panel(fig, 2, 2, 'peptide_all',   'Peptide (all heavy)', PEP_COLOR)

    # Row 3 — Sidechains
    _add_panel(fig, 3, 1, 'protein_sidechain', 'Protein sidechains',  PROT_COLOR)
    if has_peptide:
        _add_panel(fig, 3, 2, 'peptide_sidechain', 'Peptide sidechains', PEP_COLOR)

    # Axis labels
    for row in range(1, 4):
        for col in range(1, n_cols + 1):
            fig.update_xaxes(title_text='Residue', row=row, col=col)
            fig.update_yaxes(
                title_text='RMSF (Å)',
                row=row, col=col,
                gridcolor='rgba(128,128,128,0.15)',
                gridwidth=1,
                zeroline=False,
            )

    # Summary annotation
    summary_bits: List[str] = []
    for key, lbl in [('protein_ca', 'Protein Cα'), ('peptide_ca', 'Peptide Cα')]:
        if key in rmsf_data:
            dt = rmsf_data[key]
            v = (np.asarray(dt[2] if len(dt) == 3 else dt[1], dtype=float)) * 10.0
            if v.size > 0:
                summary_bits.append(
                    f"{lbl} μ={np.mean(v):.2f} Å · p95={np.percentile(v, 95):.2f} Å"
                )

    fig.update_layout(
        title=dict(
            text=f'<b>{label.upper()}</b> · RMSF Profiles',
            x=0.5,
            xanchor='center',
        ),
        template=PLOTLY_TEMPLATE,
        width=FIG_WIDTH * n_cols,
        height=320 * 3 + 120,
        showlegend=False,
        margin=dict(t=90, l=60, r=30, b=90),
        plot_bgcolor='white',
    )

    if summary_bits:
        fig.add_annotation(
            text='  |  '.join(summary_bits),
            align='left',
            showarrow=False,
            xref='paper',
            yref='paper',
            x=0.01,
            y=1.055,
            font=dict(size=11),
        )

    if output_html:
        output_html.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(str(output_html))
        print(f'Saved: {output_html}')

    if output_png:
        output_png.parent.mkdir(parents=True, exist_ok=True)
        _safe_write_image(fig, output_png, scale=2)

    return fig


def save_rmsd_data(
    rmsd_data: Dict[str, np.ndarray],
    time_ps: np.ndarray,
    output_path: Path = None
):
    """Save RMSD data to CSV."""
    if output_path is None:
        output_path = DATA_OUTPUT_DIR / "rmsd_data.csv"
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    df = pd.DataFrame({'time_ps': time_ps, 'time_ns': time_ps / 1000})
    for name, values in rmsd_data.items():
        df[f'rmsd_{name}_nm'] = values
        df[f'rmsd_{name}_A'] = values * 10
    
    df.to_csv(output_path, index=False)
    print(f"Saved: {output_path}")


def save_rmsf_data(
    rmsf_data: Dict[str, Any],
    output_path: Path = None,
    topology: Optional[md.Topology] = None,
):
    """Save RMSF data to CSV.

    Accepts 3-tuple (resSeq, resname, rmsf) from compute_all_rmsf(),
    or legacy 2-tuple (resSeq, rmsf).
    Backwards-compat alias keys ('protein', 'peptide', 'peptide_relative')
    are skipped to avoid duplicate rows.
    """
    if output_path is None:
        output_path = DATA_OUTPUT_DIR / "rmsf_data.csv"

    output_path.parent.mkdir(parents=True, exist_ok=True)

    _SKIP_ALIASES = frozenset({'protein', 'peptide', 'peptide_relative'})

    rows = []
    for chain_name, data_tuple in rmsf_data.items():
        if chain_name in _SKIP_ALIASES:
            continue  # Skip backward-compat aliases (duplicates)
        if len(data_tuple) == 3:
            residues, resnames, rmsf = data_tuple
        else:
            residues, rmsf = data_tuple
            resnames = np.array([''] * len(residues), dtype='U8')

        for i, (res, val) in enumerate(zip(residues, rmsf)):
            res_int = int(res)
            rnm = str(resnames[i]).strip() if i < len(resnames) else ''
            rows.append({
                'chain': chain_name,
                'residue': res_int,
                'residue_name': rnm,
                'residue_label': f"{rnm}{res_int}" if rnm else str(res_int),
                'rmsf_nm': float(val),
                'rmsf_A': float(val) * 10.0,
            })

    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    print(f"Saved: {output_path}")


# ============================================================================
# MAIN
# ============================================================================
def run_rmsd_analysis(
    config_or_traj: Union[AnalysisConfig, md.Trajectory, None] = None,
    *,
    stride: int = 1,
    skip_rmsf: bool = False,
    skip_plots: bool = False,
    skip_png: bool = False,
):
    """Run RMSD/RMSF analysis using either a config or a pre-loaded trajectory."""

    if isinstance(config_or_traj, AnalysisConfig):
        return _run_rmsd_analysis_impl(
            config=config_or_traj,
            traj=None,
            stride=stride,
            skip_rmsf=skip_rmsf,
            skip_plots=skip_plots,
            skip_png=skip_png,
        )

    return _run_rmsd_analysis_impl(
        config=None,
        traj=config_or_traj,
        stride=stride,
        skip_rmsf=skip_rmsf,
        skip_plots=skip_plots,
        skip_png=skip_png,
    )


def _run_rmsd_analysis_impl(
    *,
    config: Optional[AnalysisConfig],
    traj: Optional[md.Trajectory],
    stride: int = 1,
    skip_rmsf: bool = False,
    skip_plots: bool = False,
    skip_png: bool = False,
):
    print("=" * 60)
    print("RMSD/RMSF ANALYSIS")
    print("=" * 60)

    label = config.label if config else "default"
    analysis_root = (
        config.output_root / label / "rmsd"
        if config
        else (FIGURES_DIR.parent if FIGURES_DIR else Path("outputs")) / "rmsd" / label
    )
    figures_dir = analysis_root / "figures"
    html_dir = analysis_root / "html"
    data_dir = analysis_root / "data"
    for directory in (figures_dir, html_dir, data_dir):
        directory.mkdir(parents=True, exist_ok=True)

    if traj is None:
        topology_path = config.topology if config else None
        if config is not None and topology_path is not None:
            prepared_dir = analysis_root / "_prepared_topology"
            prepared = _prepare_topology_for_mdtraj(Path(topology_path), prepared_dir)
            if prepared != Path(topology_path):
                print(f"\n[WARN] Using prepared mdtraj-safe topology: {prepared}")
            topology_path = prepared

        traj = load_and_concatenate_trajectories(
            traj_paths=config.trajectories if config else None,
            topology=topology_path,
            stride=config.stride if config else stride,
            time_step_ps=config.time_step_ps if config else SAVE_FREQ_PS,
            box_lengths_nm=config.box_lengths_nm if config else None,
            autoimage_selection=config.autoimage_selection if config else None,
            apply_autoimage=config.should_autoimage() if config else False,
        )

    if traj.n_frames < 2:
        eff_stride = config.stride if config else stride
        raise RuntimeError(
            f"Trajectory loaded with only {traj.n_frames} frame(s). "
            f"This makes RMSD/RMSF degenerate (often flat-zero). "
            f"Rerun with a smaller --stride (current: {eff_stride})."
        )

    unitcell_lengths = getattr(traj, 'unitcell_lengths', None)
    has_box = (
        unitcell_lengths is not None
        and np.size(unitcell_lengths) > 0
        and np.all(np.isfinite(unitcell_lengths))
        and np.any(np.asarray(unitcell_lengths) > 0)
    )
    if not has_box:
        print(
            "\nWARNING: No periodic box vectors detected on the trajectory. "
            "Protein-peptide relative RMSD/RMSF can be inflated by PBC wrapping/jumps. "
            "Consider rerunning with `--box-xml` (OpenMM State XML) and autoimage enabled."
        )
    elif config is not None and not config.should_autoimage():
        print(
            "\nWARNING: Box vectors present but autoimage is disabled. "
            "Relative peptide metrics may still include PBC wrapping artifacts."
        )

    print("\nExtracting protein and peptide atoms...")
    traj_pp = extract_protein_peptide(traj)
    print(f"Reduced system: {traj_pp.n_atoms} atoms")

    if traj_pp.n_frames < 2:
        eff_stride = config.stride if config else stride
        raise RuntimeError(
            f"Protein+peptide sliced trajectory has only {traj_pp.n_frames} frame(s). "
            f"Rerun with a smaller --stride (current: {eff_stride})."
        )

    frame_times_ps, frame_dt_ps = _resolve_time_axis(traj_pp, config)
    print(f"Resolved frame time step: {frame_dt_ps:.3f} ps")

    print("\nComputing RMSD...")
    rmsd_data = compute_all_rmsd(traj_pp)

    for name, values in rmsd_data.items():
        mean_val = np.mean(values) * 10
        max_val = np.max(values) * 10
        print(f"  {name}: mean={mean_val:.2f} Å, max={max_val:.2f} Å")

    rmsf_data = None
    if not skip_rmsf:
        print("\nComputing RMSF...")
        rmsf_data = compute_all_rmsf(traj_pp)
        if rmsf_data:
            _SKIP_ALIASES = frozenset({'protein', 'peptide', 'peptide_relative'})
            print("RMSF summary (Å):")
            for key, data_tuple in rmsf_data.items():
                if key in _SKIP_ALIASES:
                    continue
                # Accept 3-tuple (new) or 2-tuple (legacy)
                if len(data_tuple) == 3:
                    _residues, _resnames, rmsf_vals = data_tuple
                else:
                    _residues, rmsf_vals = data_tuple
                values_A = np.asarray(rmsf_vals, dtype=float) * 10.0
                if values_A.size == 0:
                    continue
                mean_val = float(np.mean(values_A))
                p95_val = float(np.percentile(values_A, 95))
                max_val = float(np.max(values_A))
                print(f"  {key}: mean={mean_val:.2f} Å, p95={p95_val:.2f} Å, max={max_val:.2f} Å")

    print("\nSaving data...")
    save_rmsd_data(
        rmsd_data,
        frame_times_ps,
        output_path=data_dir / f"{label}_rmsd.csv",
    )
    if rmsf_data is not None:
        save_rmsf_data(
            rmsf_data,
            output_path=data_dir / f"{label}_rmsf.csv",
            topology=traj_pp.topology,
        )

    if not skip_plots:
        print("\nCreating plots...")
        plot_rmsd_timeseries(
            rmsd_data,
            frame_times_ps,
            label=label,
            output_html=html_dir / f"{label}_rmsd_timeseries.html",
            output_png=None if skip_png else (figures_dir / f"{label}_rmsd_timeseries.png"),
        )

        if rmsf_data is not None:
            plot_rmsf(
                rmsf_data,
                label=label,
                output_html=html_dir / f"{label}_rmsf.html",
                output_png=None if skip_png else (figures_dir / f"{label}_rmsf.png"),
                topology=traj_pp.topology,
            )
    else:
        print("\n[OK] Plots skipped (--skip-plots)")

    print("\n[OK] RMSD/RMSF analysis complete!")

    return rmsd_data, rmsf_data


if __name__ == "__main__":
    parser = build_cli_parser(
        "Protein-peptide RMSD/RMSF analysis",
        env_prefix="CCT",
    )
    parser.add_argument(
        "--skip-rmsf",
        action="store_true",
        help="Skip RMSF computation/plots for a faster RMSD-only run",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Skip Plotly HTML/PNG generation (writes CSVs only; much faster)",
    )
    parser.add_argument(
        "--skip-png",
        action="store_true",
        help="Skip static PNG export (Kaleido); keep HTML plots and CSV outputs",
    )
    args = parser.parse_args()

    cfg = resolve_config(
        args,
        default_topology=TOPOLOGY_PDB,
        default_trajectories=TRAJECTORIES,
        default_output=OUTPUT_DIR,
        default_stride=1,
        default_time_step_ps=SAVE_FREQ_PS,
        default_label=args.label or "rmsd",
    )

    run_rmsd_analysis(
        cfg,
        skip_rmsf=bool(getattr(args, "skip_rmsf", False)),
        skip_plots=bool(getattr(args, "skip_plots", False)),
        skip_png=bool(getattr(args, "skip_png", False)),
    )

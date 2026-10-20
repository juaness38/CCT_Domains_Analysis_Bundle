"""
CCT-OSR1 Binding Analysis - Water Interface Analysis
=====================================================
Analyze water molecules at the protein-peptide interface.

Includes:
- Bridging water detection
- Water residence time calculation (autocorrelation method)
- Water exchange kinetics

Based on:
- Impey, Madden & McDonald, J. Phys. Chem. 87, 5071 (1983)
- Garcia & Stiller, J. Comput. Chem. 14, 1396 (1993)
"""

import numpy as np
import mdtraj as md
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
from pathlib import Path
from typing import Optional, Tuple, Dict, List, Any
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
import warnings
import logging
import json
import os

logger = logging.getLogger(__name__)

try:
    from .trajectory_loader import get_protein_peptide_indices, _resolve_frame_dt_ps
except ImportError:  # pragma: no cover
    from trajectory_loader import get_protein_peptide_indices, _resolve_frame_dt_ps

try:
    from .config import (
        SAVE_FREQ_PS, PLOTLY_TEMPLATE, COLOR_SCHEME,
        FIG_WIDTH, FIG_HEIGHT, FIG_HEIGHT_TALL,
        HTML_DIR, FIGURES_DIR, DATA_OUTPUT_DIR,
        RFXV_MOTIF,
        TRAJECTORIES, TOPOLOGY_PDB
    )
except ImportError:  # pragma: no cover
    from config import (
        SAVE_FREQ_PS, PLOTLY_TEMPLATE, COLOR_SCHEME,
        FIG_WIDTH, FIG_HEIGHT, FIG_HEIGHT_TALL,
        HTML_DIR, FIGURES_DIR, DATA_OUTPUT_DIR,
        RFXV_MOTIF,
        TRAJECTORIES, TOPOLOGY_PDB
    )

def _prepare_selection_indices(topology_path: Optional[Path] = None) -> Dict[str, np.ndarray]:
    """Load topology once and collect reusable atom selections."""
    top_path = Path(topology_path) if topology_path else Path(TOPOLOGY_PDB)
    reference = md.load(str(top_path))
    top = reference.topology
    try:
        protein_idx, peptide_idx = get_protein_peptide_indices(reference)
    except Exception:
        protein_idx = top.select('protein and chainid 0')
        peptide_idx = top.select('protein and chainid 1')

    heavy_atoms = top.select('not element H')
    protein_heavy = np.intersect1d(np.asarray(protein_idx, dtype=int), heavy_atoms, assume_unique=False)
    peptide_heavy = np.intersect1d(np.asarray(peptide_idx, dtype=int), heavy_atoms, assume_unique=False)

    selections: Dict[str, np.ndarray] = {
        'water_oxygens': top.select('water and name O'),
        'protein_heavy': protein_heavy,
        'peptide_heavy': peptide_heavy,
    }
    if len(RFXV_MOTIF['residues']) > 0:
        start = min(RFXV_MOTIF['residues'])
        end = max(RFXV_MOTIF['residues'])
        query = f"protein and resSeq {start} to {end} and not element H"
        rfxv_atoms = top.select(query)
        rfxv_atoms = np.intersect1d(rfxv_atoms, peptide_heavy, assume_unique=False)
        if rfxv_atoms.size == 0:
            warnings.warn(
                "RFxV resSeq selection returned no atoms; falling back to relative index",
                RuntimeWarning
            )
            peptide_residue_ids = sorted({top.atom(int(ai)).residue.index for ai in peptide_idx})
            if peptide_residue_ids:
                fallback_ids = set(peptide_residue_ids[:5])
                rfxv_atoms = np.array(
                    [
                        atom.index for atom in top.atoms
                        if atom.residue.index in fallback_ids
                        and (atom.element is None or atom.element.symbol != 'H')
                    ],
                    dtype=int,
                )
            else:
                rfxv_atoms = np.array([], dtype=int)
    else:
        rfxv_atoms = np.array([], dtype=int)
    selections['rfxv_atoms'] = rfxv_atoms
    return selections


def _resolve_output_directories(output_dir: Optional[Path] = None) -> Tuple[Path, Path, Path]:
    """Return HTML/figure/data directories for either default config or a custom root."""
    if output_dir is None:
        return HTML_DIR, FIGURES_DIR, DATA_OUTPUT_DIR

    root = Path(output_dir)
    html_dir = root / 'html'
    figures_dir = root / 'figures'
    data_dir = root / 'data'
    html_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    return html_dir, figures_dir, data_dir


def _safe_write_plot_outputs(
    fig: go.Figure,
    output_html: Optional[Path] = None,
    output_png: Optional[Path] = None,
) -> None:
    """Write HTML always and PNG opportunistically if the image backend exists."""
    if output_html:
        output_html.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(str(output_html))
        logger.info(f"Saved: {output_html}")

    if output_png and os.getenv('SMIC_ENABLE_STATIC_IMAGE_EXPORT', '').strip() == '1':
        output_png.parent.mkdir(parents=True, exist_ok=True)
        try:
            fig.write_image(str(output_png), scale=2)
            logger.info(f"Saved: {output_png}")
        except Exception as exc:  # pragma: no cover - depends on optional backends
            logger.warning(f"Skipping PNG export for {output_png}: {exc}")
    elif output_png:
        logger.info(f"Skipping PNG export for {output_png} (set SMIC_ENABLE_STATIC_IMAGE_EXPORT=1 to enable)")


def _coerce_chain_label(chain: Any) -> str:
    """Return a stable human-readable chain label across MDTraj topology variants."""
    for attr in ('chain_id', 'id', 'segment_id'):
        value = getattr(chain, attr, None)
        if value not in (None, ''):
            return str(value)
    return str(chain.index)


def _prepare_solvent_accessibility_context(
    topology_path: Optional[Path] = None,
    pocket_residues: Optional[List[int]] = None,
    pocket_cutoff_nm: float = 0.45,
) -> Dict[str, Any]:
    """Build reusable topology metadata for streamed solvent-accessibility analysis."""
    top_path = Path(topology_path) if topology_path else Path(TOPOLOGY_PDB)
    reference = md.load(str(top_path))
    top = reference.topology

    try:
        protein_idx, peptide_idx = get_protein_peptide_indices(reference)
    except Exception:
        protein_idx = top.select('protein and chainid 0')
        peptide_idx = top.select('protein and chainid 1')

    protein_idx = np.asarray(protein_idx, dtype=int)
    peptide_idx = np.asarray(peptide_idx, dtype=int)
    solute_indices = np.unique(np.concatenate([protein_idx, peptide_idx])).astype(int)
    if solute_indices.size == 0:
        raise ValueError('Could not resolve any protein/peptide atoms for solvent accessibility analysis')

    solute_reference = reference.atom_slice(solute_indices)
    solute_top = solute_reference.topology
    global_to_local = {int(global_idx): int(local_idx) for local_idx, global_idx in enumerate(solute_indices.tolist())}

    protein_local = np.array([global_to_local[int(idx)] for idx in protein_idx if int(idx) in global_to_local], dtype=int)
    peptide_local = np.array([global_to_local[int(idx)] for idx in peptide_idx if int(idx) in global_to_local], dtype=int)
    heavy_local = solute_top.select('not element H')
    protein_heavy_local = np.intersect1d(protein_local, heavy_local, assume_unique=False)
    peptide_heavy_local = np.intersect1d(peptide_local, heavy_local, assume_unique=False)

    residue_order = list(solute_top.residues)
    residue_col_by_index = {int(residue.index): int(col) for col, residue in enumerate(residue_order)}
    atom_to_residue_col = np.array(
        [residue_col_by_index[int(atom.residue.index)] for atom in solute_top.atoms],
        dtype=int,
    )
    atom_to_residue_matrix = np.zeros((solute_top.n_atoms, len(residue_order)), dtype=float)
    atom_to_residue_matrix[np.arange(solute_top.n_atoms), atom_to_residue_col] = 1.0

    protein_local_set = set(protein_local.tolist())
    peptide_local_set = set(peptide_local.tolist())
    residue_metadata: List[Dict[str, Any]] = []
    chain_groups: Dict[str, Dict[str, Any]] = {}
    protein_residue_cols: List[int] = []
    peptide_residue_cols: List[int] = []

    for col, residue in enumerate(residue_order):
        residue_atom_indices = np.array([atom.index for atom in residue.atoms], dtype=int)
        is_peptide = any(int(atom_idx) in peptide_local_set for atom_idx in residue_atom_indices.tolist())
        is_protein = any(int(atom_idx) in protein_local_set for atom_idx in residue_atom_indices.tolist())
        entity_role = 'peptide' if is_peptide and not is_protein else 'protein'
        chain_label = _coerce_chain_label(residue.chain)
        chain_key = f"{entity_role}:{chain_label}"
        chain_entry = chain_groups.setdefault(
            chain_key,
            {
                'chain_key': chain_key,
                'chain_label': chain_label,
                'chain_index': int(residue.chain.index),
                'entity_role': entity_role,
                'residue_columns': [],
            },
        )
        chain_entry['residue_columns'].append(col)

        if entity_role == 'protein':
            protein_residue_cols.append(col)
        else:
            peptide_residue_cols.append(col)

        residue_metadata.append(
            {
                'local_residue_column': col,
                'residue_index': int(residue.index),
                'resSeq': int(getattr(residue, 'resSeq', residue.index)),
                'residue_name': residue.name,
                'chain_index': int(residue.chain.index),
                'chain_label': chain_label,
                'chain_key': chain_key,
                'entity_role': entity_role,
            }
        )

    if pocket_residues:
        requested = {int(value) for value in pocket_residues}
        pocket_residue_cols = [
            meta['local_residue_column']
            for meta in residue_metadata
            if meta['entity_role'] == 'protein'
            and (int(meta['resSeq']) in requested or int(meta['residue_index']) in requested)
        ]
        pocket_source = 'explicit_residue_list'
    else:
        if protein_heavy_local.size == 0 or peptide_heavy_local.size == 0:
            pocket_residue_cols = []
        else:
            contact_atoms = md.compute_neighbors(
                solute_reference,
                pocket_cutoff_nm,
                peptide_heavy_local,
                haystack_indices=protein_heavy_local,
            )[0]
            pocket_residue_cols = sorted(
                {
                    residue_col_by_index[int(solute_top.atom(int(atom_idx)).residue.index)]
                    for atom_idx in contact_atoms.tolist()
                }
            )
        pocket_source = 'auto_detected_peptide_contact'

    if pocket_residues and not pocket_residue_cols:
        raise ValueError(f'Pocket residues {sorted(requested)} did not match any protein residues in topology')

    return {
        'topology_path': top_path,
        'solute_indices': solute_indices,
        'solute_reference': solute_reference,
        'solute_topology': solute_top,
        'atom_to_residue_matrix': atom_to_residue_matrix,
        'residue_metadata': residue_metadata,
        'chain_groups': chain_groups,
        'pocket_residue_cols': pocket_residue_cols,
        'pocket_source': pocket_source,
        'protein_residue_cols': protein_residue_cols,
        'peptide_residue_cols': peptide_residue_cols,
    }


def _append_residue_timeseries_rows(
    rows: List[Dict[str, Any]],
    residue_sasa: np.ndarray,
    residue_metadata: List[Dict[str, Any]],
    frame_indices: np.ndarray,
    time_ps: np.ndarray,
) -> None:
    """Append long-format residue SASA records for a chunk."""
    for local_frame_idx, frame_index in enumerate(frame_indices.tolist()):
        time_ps_value = float(time_ps[local_frame_idx])
        for residue_col, meta in enumerate(residue_metadata):
            rows.append(
                {
                    'frame_index': int(frame_index),
                    'time_ps': time_ps_value,
                    'time_ns': time_ps_value / 1000.0,
                    **meta,
                    'sasa_nm2': float(residue_sasa[local_frame_idx, residue_col]),
                }
            )


def plot_chain_solvent_accessibility(
    chain_timeseries: pd.DataFrame,
    output_html: Optional[Path] = None,
    output_png: Optional[Path] = None,
) -> go.Figure:
    """Plot chain-level solvent accessibility over time."""
    fig = go.Figure()
    for chain_key in chain_timeseries['chain_key'].drop_duplicates().tolist():
        subset = chain_timeseries[chain_timeseries['chain_key'] == chain_key]
        fig.add_trace(
            go.Scatter(
                x=subset['time_ns'],
                y=subset['chain_sasa_nm2'],
                mode='lines',
                name=chain_key,
            )
        )

    fig.update_layout(
        title='Chain Solvent Accessibility Over Time',
        xaxis_title='Time (ns)',
        yaxis_title='Chain SASA (nm²)',
        template=PLOTLY_TEMPLATE,
        width=FIG_WIDTH,
        height=FIG_HEIGHT,
    )
    _safe_write_plot_outputs(fig, output_html=output_html, output_png=output_png)
    return fig


def plot_pocket_solvent_accessibility(
    pocket_timeseries: pd.DataFrame,
    output_html: Optional[Path] = None,
    output_png: Optional[Path] = None,
) -> go.Figure:
    """Plot pocket total and mean per-residue SASA over time."""
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=pocket_timeseries['time_ns'],
            y=pocket_timeseries['pocket_total_sasa_nm2'],
            mode='lines',
            name='Pocket total SASA',
        )
    )
    fig.add_trace(
        go.Scatter(
            x=pocket_timeseries['time_ns'],
            y=pocket_timeseries['pocket_mean_residue_sasa_nm2'],
            mode='lines',
            name='Pocket mean residue SASA',
            yaxis='y2',
        )
    )

    fig.update_layout(
        title='Pocket Solvent Accessibility Over Time',
        xaxis_title='Time (ns)',
        yaxis=dict(title='Pocket total SASA (nm²)'),
        yaxis2=dict(title='Mean residue SASA (nm²)', overlaying='y', side='right'),
        template=PLOTLY_TEMPLATE,
        width=FIG_WIDTH,
        height=FIG_HEIGHT,
    )
    _safe_write_plot_outputs(fig, output_html=output_html, output_png=output_png)
    return fig


def plot_residue_solvent_accessibility_summary(
    residue_summary: pd.DataFrame,
    output_html: Optional[Path] = None,
    output_png: Optional[Path] = None,
) -> go.Figure:
    """Plot mean residue SASA as a categorical bar chart."""
    labels = [
        f"{row.chain_label}:{row.residue_name}{int(row.resSeq)}"
        for row in residue_summary.itertuples(index=False)
    ]
    colors = [COLOR_SCHEME['protein'] if role == 'protein' else COLOR_SCHEME['peptide'] for role in residue_summary['entity_role']]
    fig = go.Figure(
        data=[
            go.Bar(
                x=labels,
                y=residue_summary['mean_sasa_nm2'],
                marker_color=colors,
                customdata=np.stack(
                    [
                        residue_summary['std_sasa_nm2'].to_numpy(dtype=float),
                        residue_summary['entity_role'].to_numpy(dtype=object),
                    ],
                    axis=1,
                ),
                hovertemplate='Residue=%{x}<br>Mean SASA=%{y:.3f} nm²<br>Std=%{customdata[0]:.3f} nm²<br>Role=%{customdata[1]}<extra></extra>',
            )
        ]
    )
    fig.update_layout(
        title='Mean Residue Solvent Accessibility',
        xaxis_title='Residue',
        yaxis_title='Mean SASA (nm²)',
        template=PLOTLY_TEMPLATE,
        width=max(FIG_WIDTH, 14 * max(20, len(labels))),
        height=FIG_HEIGHT,
    )
    _safe_write_plot_outputs(fig, output_html=output_html, output_png=output_png)
    return fig


def run_solvent_accessibility_analysis(
    topology_path: Optional[Path] = None,
    trajectory_paths: Optional[List[Path]] = None,
    output_dir: Optional[Path] = None,
    stride: int = 10,
    chunk_size: int = 25,
    probe_radius_nm: float = 0.14,
    n_sphere_points: int = 960,
    ps_per_frame: Optional[float] = None,
    pocket_residues: Optional[List[int]] = None,
    pocket_cutoff_nm: float = 0.45,
    save_residue_timeseries: bool = True,
    label: str = 'solvent_accessibility',
) -> Dict[str, Any]:
    """Compute streamed solvent accessibility by residue, chain, and pocket."""
    print('=' * 60)
    print('SOLVENT ACCESSIBILITY ANALYSIS')
    print('=' * 60)

    top_path = Path(topology_path) if topology_path else Path(TOPOLOGY_PDB)
    traj_paths = [Path(p) for p in trajectory_paths] if trajectory_paths is not None else [Path(p) for p in TRAJECTORIES]
    html_dir, figures_dir, data_dir = _resolve_output_directories(output_dir)
    context = _prepare_solvent_accessibility_context(
        topology_path=top_path,
        pocket_residues=pocket_residues,
        pocket_cutoff_nm=pocket_cutoff_nm,
    )

    residue_metadata = context['residue_metadata']
    atom_to_residue_matrix = context['atom_to_residue_matrix']
    chain_groups = context['chain_groups']
    pocket_residue_cols = context['pocket_residue_cols']
    solute_indices = context['solute_indices']

    residue_sum = np.zeros(len(residue_metadata), dtype=float)
    residue_sum_sq = np.zeros(len(residue_metadata), dtype=float)
    residue_min = np.full(len(residue_metadata), np.inf, dtype=float)
    residue_max = np.zeros(len(residue_metadata), dtype=float)
    residue_rows: List[Dict[str, Any]] = []
    chain_rows: List[Dict[str, Any]] = []
    pocket_rows: List[Dict[str, Any]] = []

    loaded_frame_dt_ps: Optional[float] = None
    total_frames = 0
    for traj_path in traj_paths:
        if not traj_path.exists():
            raise FileNotFoundError(f'Trajectory not found: {traj_path}')

        print(f"  Streaming {traj_path.name} for SASA (stride={stride}, chunk={chunk_size})...")
        for chunk in md.iterload(str(traj_path), top=str(top_path), chunk=chunk_size, stride=stride):
            if loaded_frame_dt_ps is None:
                loaded_frame_dt_ps = _resolve_frame_dt_ps(
                    chunk,
                    stride=stride,
                    override_time_step_ps=ps_per_frame,
                )
            if loaded_frame_dt_ps is None or loaded_frame_dt_ps <= 0:
                raise ValueError('Could not resolve a valid loaded-frame dt for solvent accessibility analysis')

            solute_chunk = chunk.atom_slice(solute_indices)
            atom_sasa = md.shrake_rupley(
                solute_chunk,
                probe_radius=probe_radius_nm,
                n_sphere_points=n_sphere_points,
                mode='atom',
            )
            residue_sasa = atom_sasa @ atom_to_residue_matrix

            frame_indices = np.arange(total_frames, total_frames + residue_sasa.shape[0], dtype=int)
            time_ps = frame_indices.astype(float) * float(loaded_frame_dt_ps)

            residue_sum += residue_sasa.sum(axis=0)
            residue_sum_sq += np.square(residue_sasa).sum(axis=0)
            residue_min = np.minimum(residue_min, residue_sasa.min(axis=0))
            residue_max = np.maximum(residue_max, residue_sasa.max(axis=0))

            if save_residue_timeseries:
                _append_residue_timeseries_rows(residue_rows, residue_sasa, residue_metadata, frame_indices, time_ps)

            for chain_key, chain_entry in chain_groups.items():
                chain_total = residue_sasa[:, chain_entry['residue_columns']].sum(axis=1)
                for local_frame_idx, frame_index in enumerate(frame_indices.tolist()):
                    chain_rows.append(
                        {
                            'frame_index': int(frame_index),
                            'time_ps': float(time_ps[local_frame_idx]),
                            'time_ns': float(time_ps[local_frame_idx]) / 1000.0,
                            'chain_key': chain_key,
                            'chain_index': int(chain_entry['chain_index']),
                            'chain_label': chain_entry['chain_label'],
                            'entity_role': chain_entry['entity_role'],
                            'chain_sasa_nm2': float(chain_total[local_frame_idx]),
                        }
                    )

            if pocket_residue_cols:
                pocket_total = residue_sasa[:, pocket_residue_cols].sum(axis=1)
                pocket_mean = residue_sasa[:, pocket_residue_cols].mean(axis=1)
                for local_frame_idx, frame_index in enumerate(frame_indices.tolist()):
                    pocket_rows.append(
                        {
                            'frame_index': int(frame_index),
                            'time_ps': float(time_ps[local_frame_idx]),
                            'time_ns': float(time_ps[local_frame_idx]) / 1000.0,
                            'pocket_source': context['pocket_source'],
                            'pocket_total_sasa_nm2': float(pocket_total[local_frame_idx]),
                            'pocket_mean_residue_sasa_nm2': float(pocket_mean[local_frame_idx]),
                            'pocket_residue_count': int(len(pocket_residue_cols)),
                        }
                    )

            total_frames += residue_sasa.shape[0]

    if total_frames == 0:
        raise RuntimeError('No frames were processed during solvent accessibility analysis')

    residue_summary = pd.DataFrame(residue_metadata)
    residue_summary['mean_sasa_nm2'] = residue_sum / float(total_frames)
    residue_summary['std_sasa_nm2'] = np.sqrt(np.clip(residue_sum_sq / float(total_frames) - np.square(residue_summary['mean_sasa_nm2']), 0.0, None))
    residue_summary['min_sasa_nm2'] = residue_min
    residue_summary['max_sasa_nm2'] = residue_max
    residue_summary['in_pocket'] = residue_summary['local_residue_column'].isin(pocket_residue_cols)
    residue_summary = residue_summary.sort_values(['chain_index', 'resSeq', 'local_residue_column']).reset_index(drop=True)

    chain_timeseries = pd.DataFrame(chain_rows)
    chain_summary = (
        chain_timeseries
        .groupby(['chain_key', 'chain_index', 'chain_label', 'entity_role'], as_index=False)['chain_sasa_nm2']
        .agg(['mean', 'std', 'min', 'max'])
        .reset_index()
        .rename(columns={'mean': 'mean_sasa_nm2', 'std': 'std_sasa_nm2', 'min': 'min_sasa_nm2', 'max': 'max_sasa_nm2'})
    )

    pocket_timeseries = pd.DataFrame(pocket_rows)
    if not pocket_timeseries.empty:
        pocket_summary = pd.DataFrame(
            [
                {
                    'pocket_source': context['pocket_source'],
                    'pocket_residue_count': int(len(pocket_residue_cols)),
                    'mean_total_sasa_nm2': float(pocket_timeseries['pocket_total_sasa_nm2'].mean()),
                    'std_total_sasa_nm2': float(pocket_timeseries['pocket_total_sasa_nm2'].std(ddof=0)),
                    'mean_residue_sasa_nm2': float(pocket_timeseries['pocket_mean_residue_sasa_nm2'].mean()),
                }
            ]
        )
        pocket_residue_summary = residue_summary[residue_summary['in_pocket']].copy()
    else:
        pocket_summary = pd.DataFrame(
            [{'pocket_source': context['pocket_source'], 'pocket_residue_count': 0, 'mean_total_sasa_nm2': 0.0, 'std_total_sasa_nm2': 0.0, 'mean_residue_sasa_nm2': 0.0}]
        )
        pocket_residue_summary = residue_summary.iloc[0:0].copy()

    metadata = {
        'label': label,
        'topology_path': str(top_path),
        'trajectory_paths': [str(path) for path in traj_paths],
        'stride': int(stride),
        'chunk_size': int(chunk_size),
        'probe_radius_nm': float(probe_radius_nm),
        'n_sphere_points': int(n_sphere_points),
        'loaded_frame_dt_ps': float(loaded_frame_dt_ps),
        'total_frames': int(total_frames),
        'total_time_ns': float(max(0, total_frames - 1) * float(loaded_frame_dt_ps) / 1000.0),
        'pocket_source': context['pocket_source'],
        'pocket_residues_resSeq': [int(row.resSeq) for row in pocket_residue_summary.itertuples(index=False)],
    }

    residue_summary_path = data_dir / f'{label}_residue_summary.csv'
    chain_timeseries_path = data_dir / f'{label}_chain_timeseries.csv'
    chain_summary_path = data_dir / f'{label}_chain_summary.csv'
    pocket_timeseries_path = data_dir / f'{label}_pocket_timeseries.csv'
    pocket_summary_path = data_dir / f'{label}_pocket_summary.csv'
    pocket_residue_summary_path = data_dir / f'{label}_pocket_residue_summary.csv'
    metadata_path = data_dir / f'{label}_metadata.json'

    residue_summary.to_csv(residue_summary_path, index=False)
    chain_timeseries.to_csv(chain_timeseries_path, index=False)
    chain_summary.to_csv(chain_summary_path, index=False)
    pocket_timeseries.to_csv(pocket_timeseries_path, index=False)
    pocket_summary.to_csv(pocket_summary_path, index=False)
    pocket_residue_summary.to_csv(pocket_residue_summary_path, index=False)
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding='utf-8')

    if save_residue_timeseries:
        residue_timeseries = pd.DataFrame(residue_rows)
        residue_timeseries_path = data_dir / f'{label}_residue_timeseries.csv'
        residue_timeseries.to_csv(residue_timeseries_path, index=False)
    else:
        residue_timeseries = pd.DataFrame()
        residue_timeseries_path = None

    plot_chain_solvent_accessibility(
        chain_timeseries,
        output_html=html_dir / f'{label}_chain_sasa.html',
        output_png=figures_dir / f'{label}_chain_sasa.png',
    )
    if not pocket_timeseries.empty:
        plot_pocket_solvent_accessibility(
            pocket_timeseries,
            output_html=html_dir / f'{label}_pocket_sasa.html',
            output_png=figures_dir / f'{label}_pocket_sasa.png',
        )
    plot_residue_solvent_accessibility_summary(
        residue_summary,
        output_html=html_dir / f'{label}_residue_sasa_summary.html',
        output_png=figures_dir / f'{label}_residue_sasa_summary.png',
    )

    print(f"  Processed {total_frames} frames ({metadata['total_time_ns']:.2f} ns)")
    print(f"  Pocket mode: {context['pocket_source']} ({len(pocket_residue_cols)} residues)")
    print(f"  Saved residue summary: {residue_summary_path}")
    print(f"  Saved chain summary:   {chain_summary_path}")
    print(f"  Saved pocket summary:  {pocket_summary_path}")

    return {
        'residue_summary': residue_summary,
        'residue_timeseries': residue_timeseries,
        'chain_timeseries': chain_timeseries,
        'chain_summary': chain_summary,
        'pocket_timeseries': pocket_timeseries,
        'pocket_summary': pocket_summary,
        'pocket_residue_summary': pocket_residue_summary,
        'metadata': metadata,
        'paths': {
            'residue_summary': residue_summary_path,
            'residue_timeseries': residue_timeseries_path,
            'chain_timeseries': chain_timeseries_path,
            'chain_summary': chain_summary_path,
            'pocket_timeseries': pocket_timeseries_path,
            'pocket_summary': pocket_summary_path,
            'pocket_residue_summary': pocket_residue_summary_path,
            'metadata': metadata_path,
        },
    }


def _count_bridging_waters(protein_hits: np.ndarray, peptide_hits: np.ndarray) -> int:
    """Count water oxygens that appear in both neighbor lists."""
    if protein_hits.size == 0 or peptide_hits.size == 0:
        return 0
    protein_set = set(protein_hits.tolist())
    peptide_set = set(peptide_hits.tolist())
    return len(protein_set & peptide_set)


def stream_water_counts(
    stride: int = 10,
    cutoff_nm: float = 0.35,
    chunk_size: int = 250,
    topology_path: Optional[Path] = None,
    trajectory_paths: Optional[List[Path]] = None,
    ps_per_frame: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Stream trajectories chunk-by-chunk and compute water metrics."""
    top_path = Path(topology_path) if topology_path else Path(TOPOLOGY_PDB)
    traj_paths = [Path(p) for p in trajectory_paths] if trajectory_paths is not None else [Path(p) for p in TRAJECTORIES]

    selections = _prepare_selection_indices(top_path)
    interface_counts: List[int] = []
    rfxv_counts: List[int] = []

    total_frames = 0
    frame_dt_ps = None
    for traj_path in traj_paths:
        if not traj_path.exists():
            raise FileNotFoundError(f"Trajectory not found: {traj_path}")
        print(f"  Streaming {traj_path.name} (stride={stride})...")
        for chunk in md.iterload(
            str(traj_path),
            top=str(top_path),
            chunk=chunk_size,
            stride=stride
        ):
            if frame_dt_ps is None:
                frame_dt_ps = _resolve_frame_dt_ps(
                    chunk,
                    stride=stride,
                    override_time_step_ps=ps_per_frame,
                )
            prot_hits = md.compute_neighbors(
                chunk,
                cutoff_nm,
                selections['protein_heavy'],
                haystack_indices=selections['water_oxygens']
            )
            pep_hits = md.compute_neighbors(
                chunk,
                cutoff_nm,
                selections['peptide_heavy'],
                haystack_indices=selections['water_oxygens']
            )
            if selections['rfxv_atoms'].size > 0:
                motif_hits = md.compute_neighbors(
                    chunk,
                    cutoff_nm,
                    selections['rfxv_atoms'],
                    haystack_indices=selections['water_oxygens']
                )
            else:
                motif_hits = [np.array([], dtype=int) for _ in range(chunk.n_frames)]
            for frame_idx in range(chunk.n_frames):
                interface_counts.append(
                    _count_bridging_waters(prot_hits[frame_idx], pep_hits[frame_idx])
                )
                rfxv_counts.append(len(motif_hits[frame_idx]))
            total_frames += chunk.n_frames

    if frame_dt_ps is None or frame_dt_ps <= 0:
        warnings.warn("SAVE_FREQ_PS is undefined/invalid; using 1 ps/frame fallback for water time axis", RuntimeWarning)
        frame_dt_ps = float(max(1, int(stride)))

    time_ps = np.arange(total_frames, dtype=float) * frame_dt_ps
    return time_ps, np.array(interface_counts, dtype=int), np.array(rfxv_counts, dtype=int)


def plot_interface_waters(
    interface_waters: np.ndarray,
    rfxv_waters: np.ndarray,
    time_ps: np.ndarray,
    output_html: Path = None,
    output_png: Path = None
) -> go.Figure:
    """
    Plot interfacial water counts over time.
    """
    time_ns = time_ps / 1000
    
    interface_smooth = gaussian_filter1d(interface_waters.astype(float), sigma=20)
    rfxv_smooth = gaussian_filter1d(rfxv_waters.astype(float), sigma=20)
    
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=(
            'Waters at Protein-Peptide Interface',
            'Waters around RFxV Motif'
        ),
        vertical_spacing=0.15
    )
    
    # Interface waters
    fig.add_trace(
        go.Scatter(
            x=time_ns, y=interface_waters,
            name='Interface Waters',
            line=dict(color=COLOR_SCHEME['water'], width=1),
            opacity=0.3
        ),
        row=1, col=1
    )
    fig.add_trace(
        go.Scatter(
            x=time_ns, y=interface_smooth,
            name='Interface (smoothed)',
            line=dict(color=COLOR_SCHEME['water'], width=2)
        ),
        row=1, col=1
    )
    
    # RFxV waters
    fig.add_trace(
        go.Scatter(
            x=time_ns, y=rfxv_waters,
            name='RFxV Waters',
            line=dict(color=COLOR_SCHEME['rfxv'], width=1),
            opacity=0.3
        ),
        row=2, col=1
    )
    fig.add_trace(
        go.Scatter(
            x=time_ns, y=rfxv_smooth,
            name='RFxV (smoothed)',
            line=dict(color=COLOR_SCHEME['rfxv'], width=2)
        ),
        row=2, col=1
    )
    
    fig.update_xaxes(title_text="Time (ns)", row=2, col=1)
    fig.update_yaxes(title_text="Water Count", row=1, col=1)
    fig.update_yaxes(title_text="Water Count", row=2, col=1)
    
    fig.update_layout(
        title="Water Molecules at Binding Interface",
        template=PLOTLY_TEMPLATE,
        width=FIG_WIDTH,
        height=FIG_HEIGHT_TALL,
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="center",
            x=0.5
        )
    )
    
    if output_html:
        _safe_write_plot_outputs(fig, output_html=output_html, output_png=output_png)
    
    return fig


def save_water_data(
    interface_waters: np.ndarray,
    rfxv_waters: np.ndarray,
    time_ps: np.ndarray,
    output_path: Path = None
):
    """Save water analysis data to CSV."""
    if output_path is None:
        output_path = DATA_OUTPUT_DIR / "water_analysis.csv"
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    df = pd.DataFrame({
        'time_ps': time_ps,
        'time_ns': time_ps / 1000,
        'interface_waters': interface_waters,
        'rfxv_waters': rfxv_waters
    })
    
    df.to_csv(output_path, index=False)
    print(f"Saved: {output_path}")


# ============================================================================
# WATER RESIDENCE TIME CALCULATION
# ============================================================================

def compute_water_residence_time(
    traj: md.Trajectory,
    selection_indices: np.ndarray,
    water_oxygens: np.ndarray,
    cutoff_nm: float = 0.35,
    dt_ps: float = None,
    method: str = 'autocorrelation'
) -> Tuple[float, np.ndarray, np.ndarray, dict]:
    """
    Compute water residence time near a set of atoms.
    
    Uses the survival probability autocorrelation function:
    
    C(t) = <h(0) * h(t)> / <h(0)>
    
    where h(t) = 1 if water is within cutoff at time t, 0 otherwise.
    
    The residence time τ is obtained by fitting C(t) to:
    - Single exponential: C(t) = exp(-t/τ)
    - Bi-exponential: C(t) = A*exp(-t/τ_fast) + (1-A)*exp(-t/τ_slow)
    
    Based on:
    - Impey, Madden & McDonald, J. Phys. Chem. 87, 5071 (1983)
    
    Parameters
    ----------
    traj : md.Trajectory
        Trajectory to analyze
    selection_indices : np.ndarray
        Atom indices for the region of interest (e.g., binding site)
    water_oxygens : np.ndarray
        Water oxygen atom indices
    cutoff_nm : float
        Distance cutoff (nm)
    dt_ps : float
        Time step between frames (ps). If None, uses SAVE_FREQ_PS.
    method : str
        'autocorrelation' (default) or 'mean_duration'
        
    Returns
    -------
    tau_residence : float
        Residence time (ps)
    time_lags : np.ndarray
        Time lags for autocorrelation
    survival_prob : np.ndarray
        Survival probability C(t)
    fit_params : dict
        Fit metadata from ``_fit_residence_time`` (model type and coefficients).
    """
    if dt_ps is None:
        dt_ps = SAVE_FREQ_PS
    
    n_frames = traj.n_frames
    n_waters = len(water_oxygens)
    
    if n_frames < 10:
        logger.warning("Trajectory too short for residence time calculation")
        return 0.0, np.array([0]), np.array([1.0])
    
    # Compute which waters are within cutoff at each frame
    presence = np.zeros((n_frames, n_waters), dtype=bool)
    water_index_map = {int(atom_idx): col_idx for col_idx, atom_idx in enumerate(water_oxygens)}
    
    logger.info(f"Computing water presence matrix ({n_frames} frames, {n_waters} waters)...")
    
    for frame_idx in range(n_frames):
        # Get waters near the selection at this frame
        neighbors = md.compute_neighbors(
            traj[frame_idx],
            cutoff_nm,
            selection_indices,
            haystack_indices=water_oxygens
        )
        if len(neighbors) > 0:
            neighbor_list = neighbors[0]  # Single frame
            for atom_idx in neighbor_list:
                col_idx = water_index_map.get(int(atom_idx))
                if col_idx is not None:
                    presence[frame_idx, col_idx] = True
    
    # Compute survival probability
    # C(t) = <h(0) * h(t)> / <h(0)>
    max_lag = min(n_frames // 2, 500)  # Max lag of 500 frames or half trajectory
    
    survival_prob = np.zeros(max_lag)
    time_lags = np.arange(max_lag) * dt_ps
    
    logger.info("Computing survival probability (continuous/Impey method)...")
    # Continuous survival probability: C(τ) = <h(t0) · ∏_{δ=0}^{τ} h(t0+δ)> / <h(t0)>
    # Water must be present at EVERY frame from t0 to t0+τ (Impey-Madden-McDonald 1983).
    # Vectorised via prefix sum: window_sum == (τ+1) iff all frames are occupied.
    cum = np.cumsum(presence.astype(np.int32), axis=0)   # shape (n_frames, n_waters)

    for lag in range(max_lag):
        n_t0 = n_frames - lag
        # sum(presence[t0 : t0+lag+1]) for every t0 and water using prefix sums
        end_sums    = cum[lag:lag + n_t0, :]         # cum[t0 + lag]
        start_corr  = cum[:n_t0, :]                  # cum[t0]
        pres0       = presence[:n_t0, :]             # presence[t0]
        window_sums = end_sums - start_corr + pres0  # Σ_{δ=0}^{lag} presence[t0+δ]
        h_continuous = (window_sums == (lag + 1))    # True iff water present at all lag+1 frames
        numerator   = float(np.sum(pres0 & h_continuous))
        denominator = float(np.sum(pres0))
        survival_prob[lag] = (numerator / denominator) if denominator > 0 else 0.0

    # Normalize (should be 1.0 at lag=0 by construction; guard against float drift)
    if survival_prob[0] > 0:
        survival_prob = survival_prob / survival_prob[0]

    # Fit to extract residence time
    tau_residence, fit_params = _fit_residence_time(time_lags, survival_prob)

    logger.info(f"Residence time: τ = {tau_residence:.2f} ps ({tau_residence/1000:.3f} ns)")

    return tau_residence, time_lags, survival_prob, fit_params


def _fit_residence_time(
    time_lags: np.ndarray,
    survival_prob: np.ndarray
) -> Tuple[float, dict]:
    """
    Fit survival probability to extract residence time.

    Tries bi-exponential fit first, falls back to single exponential.

    Returns
    -------
    tau_eff : float
        Effective residence time (ps) — integral mean of C(t).
    fit_params : dict
        Fit metadata: 'model', 'tau_eff', and model-specific coefficients.
    """
    # Filter out NaN and negative values
    valid = ~np.isnan(survival_prob) & (survival_prob > 0)
    t_valid = time_lags[valid]
    c_valid = survival_prob[valid]
    
    if len(c_valid) < 5:
        logger.warning("Not enough valid points for fitting")
        return 0.0, {'model': 'none', 'tau_eff': 0.0}

    # Try bi-exponential first
    # C(t) = A * exp(-t/τ_fast) + (1-A) * exp(-t/τ_slow)
    try:
        def bi_exp(t, A, tau_fast, tau_slow):
            return A * np.exp(-t / tau_fast) + (1 - A) * np.exp(-t / tau_slow)
        
        # Initial guess
        p0 = [0.5, 10.0, 100.0]  # A, tau_fast, tau_slow
        bounds = ([0, 0.1, 1], [1, 1000, 10000])
        
        popt, _ = curve_fit(bi_exp, t_valid, c_valid, p0=p0, bounds=bounds, maxfev=5000)
        A, tau_fast, tau_slow = popt

        # τ_eff = ∫₀^∞ C(t) dt = A·τ_fast + (1−A)·τ_slow  (integral mean)
        tau_eff = A * tau_fast + (1 - A) * tau_slow
        logger.info(f"Bi-exponential fit: τ_fast={tau_fast:.1f} ps, τ_slow={tau_slow:.1f} ps, A={A:.2f}")
        return tau_eff, {
            'model': 'biexponential',
            'A': float(A),
            'tau_fast': float(tau_fast),
            'tau_slow': float(tau_slow),
            'tau_eff': float(tau_eff),
        }

    except Exception as e:
        logger.warning(f"Bi-exponential fit failed: {e}, trying single exponential")

    # Single exponential fallback
    try:
        def single_exp(t, tau):
            return np.exp(-t / tau)

        popt, _ = curve_fit(single_exp, t_valid, c_valid, p0=[50.0], bounds=([0.1], [10000]), maxfev=5000)
        tau = popt[0]
        logger.info(f"Single exponential fit: τ = {tau:.1f} ps")
        return tau, {'model': 'monoexponential', 'tau': float(tau), 'tau_eff': float(tau)}

    except Exception as e:
        logger.warning(f"Single exponential fit failed: {e}")
        # Fall back to 1/e crossing time
        crossing = np.where(c_valid < np.exp(-1))[0]
        if len(crossing) > 0:
            fallback = float(t_valid[crossing[0]])
        else:
            fallback = float(t_valid[-1])
        return fallback, {'model': 'crossing', 'tau_eff': fallback}


def plot_water_residence(
    time_lags: np.ndarray,
    survival_prob: np.ndarray,
    tau_residence: float,
    output_html: Optional[Path] = None,
    output_png: Optional[Path] = None,
    title: str = "Water Residence Time",
    fit_params: Optional[dict] = None,
) -> go.Figure:
    """
    Plot water survival probability and residence time fit.

    Parameters
    ----------
    fit_params : dict, optional
        Fit metadata from ``_fit_residence_time``. When provided, the plotted
        curve reflects the actual fitted model (bi-exponential or mono-exponential).
        Falls back to a single exponential using ``tau_residence`` if absent.
    """
    fig = go.Figure()

    # Survival probability
    fig.add_trace(go.Scatter(
        x=time_lags, y=survival_prob,
        mode='lines',
        name='C(t)',
        line=dict(color='blue', width=2)
    ))

    # Fitted curve — render the actual model returned by _fit_residence_time
    t_fit = np.linspace(0, time_lags[-1], 200)
    model = (fit_params or {}).get('model', 'monoexponential')
    if model == 'biexponential' and fit_params is not None:
        A, tf, ts = fit_params['A'], fit_params['tau_fast'], fit_params['tau_slow']
        c_fit = A * np.exp(-t_fit / tf) + (1 - A) * np.exp(-t_fit / ts)
        fit_label = f'Bi-exp fit (τ_eff={tau_residence:.1f} ps)'
    else:
        tau_plot = max(tau_residence, 1e-6)  # guard against zero
        c_fit = np.exp(-t_fit / tau_plot)
        fit_label = f'Fit (τ = {tau_residence:.1f} ps)'

    fig.add_trace(go.Scatter(
        x=t_fit, y=c_fit,
        mode='lines',
        name=fit_label,
        line=dict(color='red', width=2, dash='dash')
    ))
    
    # 1/e line
    fig.add_hline(y=np.exp(-1), line_dash="dot", line_color="gray",
                  annotation_text="1/e")
    
    # τ marker
    fig.add_vline(x=tau_residence, line_dash="dash", line_color="orange",
                  annotation_text=f"τ = {tau_residence:.1f} ps")
    
    fig.update_layout(
        title=title,
        xaxis_title="Time lag (ps)",
        yaxis_title="Survival probability C(t)",
        yaxis=dict(range=[0, 1.05]),
        template=PLOTLY_TEMPLATE,
        width=FIG_WIDTH,
        height=FIG_HEIGHT,
        legend=dict(
            orientation="h",
            yanchor="bottom", y=1.02,
            xanchor="center", x=0.5
        )
    )
    
    if output_html:
        _safe_write_plot_outputs(fig, output_html=output_html, output_png=output_png)
    
    return fig


# ============================================================================
# MAIN
# ============================================================================
def run_water_analysis(
    stride: int = 10,
    chunk_size: int = 250,
    compute_residence: bool = True,
    residence_stride: int = 100,  # Use higher stride for residence time (expensive)
    topology_path: Optional[Path] = None,
    trajectory_paths: Optional[List[Path]] = None,
    ps_per_frame: Optional[float] = None,
    output_dir: Optional[Path] = None,
):
    """
    Run water interface analysis using streaming neighbor search.
    
    Parameters
    ----------
    stride : int
        Frame stride for water count analysis
    chunk_size : int
        Chunk size for streaming
    compute_residence : bool
        Whether to compute water residence time (expensive)
    residence_stride : int
        Stride for residence time calculation (higher = faster)
    """
    print("=" * 60)
    print("WATER INTERFACE ANALYSIS")
    print("=" * 60)
    print(f"Stride: {stride}")
    print(f"Chunk size: {chunk_size} frames")

    top_path = Path(topology_path) if topology_path else Path(TOPOLOGY_PDB)
    traj_paths = [Path(p) for p in trajectory_paths] if trajectory_paths is not None else [Path(p) for p in TRAJECTORIES]
    html_dir, figures_dir, data_dir = _resolve_output_directories(output_dir)
    
    print("\nStreaming trajectories and computing neighbor lists...")
    time_ps, interface_waters, rfxv_waters = stream_water_counts(
        stride=stride,
        chunk_size=chunk_size,
        topology_path=top_path,
        trajectory_paths=traj_paths,
        ps_per_frame=ps_per_frame,
    )
    total_ns = time_ps[-1] / 1000 if len(time_ps) else 0.0
    print(f"  Frames processed: {len(time_ps)} ({total_ns:.1f} ns)")
    print(f"  Mean interface waters: {interface_waters.mean():.1f}")
    print(f"  Mean RFxV waters: {rfxv_waters.mean():.1f}")
    
    print("\nSaving data...")
    save_water_data(interface_waters, rfxv_waters, time_ps, output_path=data_dir / 'water_analysis.csv')
    
    print("\nCreating plots...")
    plot_interface_waters(
        interface_waters,
        rfxv_waters,
        time_ps,
        output_html=html_dir / "08_water_interface.html",
        output_png=figures_dir / "08_water_interface.png"
    )
    
    # Water residence time calculation
    residence_results = {}
    if compute_residence:
        print("\n" + "-" * 40)
        print("WATER RESIDENCE TIME CALCULATION")
        print("-" * 40)
        
        # Load a trajectory chunk for residence time
        # Use a subset for efficiency
        print(f"Loading trajectory subset (stride={residence_stride})...")
        traj_path = traj_paths[0] if len(traj_paths) > 0 else None
        
        if traj_path and Path(traj_path).exists():
            traj = md.load(str(traj_path), top=str(top_path), stride=residence_stride)
            print(f"  Loaded {traj.n_frames} frames")
            
            # Get selections
            selections = _prepare_selection_indices(top_path)
            
            # Compute residence time at interface — use all protein heavy atoms (not an arbitrary subset)
            print("\nComputing residence time at interface...")
            interface_atoms = np.concatenate([
                selections['protein_heavy'],
                selections['peptide_heavy']
            ])
            # Resolve effective dt for the strided trajectory
            _dt_res = (float(ps_per_frame) * residence_stride) if ps_per_frame else (
                (float(traj.timestep)) if (traj.timestep and traj.timestep > 0) else (
                    (float(SAVE_FREQ_PS) * residence_stride) if SAVE_FREQ_PS else None
                )
            )
            tau_interface, t_lags, c_interface, _fp_interface = compute_water_residence_time(
                traj,
                interface_atoms,
                selections['water_oxygens'],
                dt_ps=_dt_res,
            )
            residence_results['interface'] = {
                'tau_ps': tau_interface,
                'tau_ns': tau_interface / 1000
            }

            # Plot residence time with actual fitted model
            plot_water_residence(
                t_lags, c_interface, tau_interface,
                output_html=html_dir / "09_water_residence.html",
                output_png=figures_dir / "09_water_residence.png",
                title="Water Residence Time at Binding Interface",
                fit_params=_fp_interface,
            )
            
            # Save residence data
            res_df = pd.DataFrame({
                'time_lag_ps': t_lags,
                'survival_prob': c_interface
            })
            res_df.to_csv(data_dir / "water_residence.csv", index=False)
            print(f"Saved: {data_dir / 'water_residence.csv'}")
            
            # Also compute for RFxV motif if atoms exist
            if selections['rfxv_atoms'].size > 0:
                print("\nComputing residence time at RFxV motif...")
                tau_rfxv, t_lags_rfxv, c_rfxv, _fp_rfxv = compute_water_residence_time(
                    traj,
                    selections['rfxv_atoms'],
                    selections['water_oxygens'],
                    dt_ps=_dt_res,
                )
                residence_results['rfxv'] = {
                    'tau_ps': tau_rfxv,
                    'tau_ns': tau_rfxv / 1000
                }

                plot_water_residence(
                    t_lags_rfxv, c_rfxv, tau_rfxv,
                    output_html=html_dir / "10_water_residence_rfxv.html",
                    output_png=figures_dir / "10_water_residence_rfxv.png",
                    title="Water Residence Time at RFxV Motif",
                    fit_params=_fp_rfxv,
                )
        else:
            print("Warning: Trajectory not found for residence time calculation")
    
    print("\n✓ Water analysis complete!")
    
    return interface_waters, rfxv_waters, residence_results


if __name__ == "__main__":
    run_water_analysis(stride=10, compute_residence=True)

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.io import to_html


OUTPUT_ROOT = Path(r"c:\Users\busta\Downloads\MICA\_replot_md_controls_20260327")

GREEN_DARK = "#6FAA55"
GREEN_LIGHT = "#9ED681"
BLUE_DARK = "#6DB7E8"
BLUE_LIGHT = "#A8D8F5"
GRID_COLOR = "rgba(18, 52, 86, 0.10)"


@dataclass(frozen=True)
class RmsfChainGroup:
    chain_ids: Tuple[int, ...]
    chain_key: str
    display_name: str
    color: str


@dataclass(frozen=True)
class ReplicaSpec:
    system_label: str
    replica_label: str
    source_root: Path
    preferred_prefix: Optional[str] = None
    runtime_end_ns: Optional[float] = None
    runtime_source: str = "rmsd_csv_native"
    topology_path: Optional[Path] = None
    trajectory_paths: Sequence[Path] = ()
    rmsf_chain_groups: Sequence[RmsfChainGroup] = ()


REPLICA_SPECS: Sequence[ReplicaSpec] = (
    ReplicaSpec(
        system_label="OSR1 CCT",
        replica_label="R1",
        source_root=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_1_bound_OSR1(cct)-WNK4(RfXV)-PDB_2v3s\outputs\replica_1"),
        runtime_end_ns=202.4999997795523,
        runtime_source="raw_openmm_energy_csv:replica_1_2v3s_energy.csv",
        topology_path=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_1_bound_OSR1(cct)-WNK4(RfXV)-PDB_2v3s\outputs\replica_1\rmsd\_prepared_topology\replica_1_equilibrated_2vs3__sanres_nocryst.pdb"),
        trajectory_paths=(
            Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_1_bound_OSR1(cct)-WNK4(RfXV)-PDB_2v3s\replica_1\replica_1_inicio.dcd"),
            Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_1_bound_OSR1(cct)-WNK4(RfXV)-PDB_2v3s\replica_1\continuacion\replica_1_2vs3_.dcd"),
        ),
        rmsf_chain_groups=(
            RmsfChainGroup((0,), "cct_a_ca", "CCT A (pair A-D)", GREEN_DARK),
            RmsfChainGroup((3,), "rfxv_d_ca", "RFxV D (pair A-D)", BLUE_DARK),
            RmsfChainGroup((1,), "cct_b_ca", "CCT B (pair B-C)", GREEN_LIGHT),
            RmsfChainGroup((2,), "rfxv_c_ca", "RFxV C (pair B-C)", BLUE_LIGHT),
        ),
    ),
    ReplicaSpec(
        system_label="OSR1 CCT",
        replica_label="R2",
        source_root=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_1_bound_OSR1(cct)-WNK4(RfXV)-PDB_2v3s\outputs\replica_2"),
        runtime_end_ns=176.43999986904892,
        runtime_source="raw_openmm_energy_csv:replica_2_2v3s_energy.csv",
        topology_path=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_1_bound_OSR1(cct)-WNK4(RfXV)-PDB_2v3s\outputs\replica_2\rmsd\_prepared_topology\replica_2_equilibrated_inicio_2vs3__sanres_nocryst.pdb"),
        trajectory_paths=(
            Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_1_bound_OSR1(cct)-WNK4(RfXV)-PDB_2v3s\replica_2\replica_2_full_combined.dcd"),
        ),
        rmsf_chain_groups=(
            RmsfChainGroup((0,), "cct_a_ca", "CCT A (pair A-D)", GREEN_DARK),
            RmsfChainGroup((3,), "rfxv_d_ca", "RFxV D (pair A-D)", BLUE_DARK),
            RmsfChainGroup((1,), "cct_b_ca", "CCT B (pair B-C)", GREEN_LIGHT),
            RmsfChainGroup((2,), "rfxv_c_ca", "RFxV C (pair B-C)", BLUE_LIGHT),
        ),
    ),
    ReplicaSpec(
        system_label="SPAK CCT",
        replica_label="R1",
        source_root=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_9_bound_SPAK-WNK4\outputs\control_9_r1_rmsd"),
        preferred_prefix="control_9_r1",
        topology_path=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_9_bound_SPAK-WNK4\outputs\control_9_r1_rmsd\control_9_r1\rmsd\_prepared_topology\replica_1_equilibrated__nocryst.pdb"),
        trajectory_paths=(
            Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_9_bound_SPAK-WNK4\replica_1-20260314T222347Z-1-002\replica_1-001.dcd"),
        ),
    ),
    ReplicaSpec(
        system_label="SPAK CCT",
        replica_label="R2",
        source_root=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_9_bound_SPAK-WNK4\outputs\control_9_r2_rmsd"),
        preferred_prefix="control_9_r2",
        topology_path=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_9_bound_SPAK-WNK4\outputs\control_9_r2_rmsd\control_9_r2\rmsd\_prepared_topology\replica_2_equilibrated__nocryst.pdb"),
        trajectory_paths=(
            Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_9_bound_SPAK-WNK4\replica_2-20260314T225006Z-1-002\replica_2-001.dcd"),
        ),
    ),
    ReplicaSpec(
        system_label="WNK2 CCTL1",
        replica_label="R1",
        source_root=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_bound_2_WNK2(cctl1)-WNK1(17aa-RFxV)-PDB_6FBK\outputs\replica_1\rmsd"),
        runtime_end_ns=184.1799998425,
        runtime_source="raw_openmm_continuation_csv:replica_1_6FBK_cont.csv",
        topology_path=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_bound_2_WNK2(cctl1)-WNK1(17aa-RFxV)-PDB_6FBK\outputs\replica_1\rmsd\_prepared_topology\replica_1_6fbk_equilibrated__nocryst.pdb"),
        trajectory_paths=(
            Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_bound_2_WNK2(cctl1)-WNK1(17aa-RFxV)-PDB_6FBK\results_6FBK-20251130T005949Z-1-002\results_6FBK\replica_1_6FBK.dcd"),
            Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_bound_2_WNK2(cctl1)-WNK1(17aa-RFxV)-PDB_6FBK\replica_1\replica_1_6FBK_cont-003.dcd"),
        ),
    ),
    ReplicaSpec(
        system_label="WNK2 CCTL1",
        replica_label="R2",
        source_root=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_bound_2_WNK2(cctl1)-WNK1(17aa-RFxV)-PDB_6FBK\outputs\replica_2\rmsd"),
        runtime_end_ns=182.6799998476,
        runtime_source="raw_openmm_continuation_csv:replica_2_6FBK_cont.csv",
        topology_path=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_bound_2_WNK2(cctl1)-WNK1(17aa-RFxV)-PDB_6FBK\outputs\replica_2\rmsd\_prepared_topology\replica_2_6FBK_equilibrated__nocryst.pdb"),
        trajectory_paths=(
            Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_bound_2_WNK2(cctl1)-WNK1(17aa-RFxV)-PDB_6FBK\results_6FBK-20251130T005949Z-1-002\results_6FBK\replica_2_6FBK.dcd"),
            Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\Control_bound_2_WNK2(cctl1)-WNK1(17aa-RFxV)-PDB_6FBK\replica_2\replica_2_6FBK_cont-001.dcd"),
        ),
    ),
    ReplicaSpec(
        system_label="NRBP1 CCT",
        replica_label="R1",
        source_root=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\NRBP1_replica1\NRBP1_rep1\smic_rmsd_rmsf_outputs\nrbp1_replica1"),
        preferred_prefix="nrbp1_replica1",
        topology_path=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\NRBP1_replica1\NRBP1_rep1\smic_rmsd_rmsf_outputs\nrbp1_replica1\rmsd\_prepared_topology\replica_1_equilibrated__nocryst.pdb"),
        trajectory_paths=(
            Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\NRBP1_replica1\NRBP1_rep1\replica_1-001.dcd"),
        ),
    ),
    ReplicaSpec(
        system_label="NRBP1 CCT",
        replica_label="R2",
        source_root=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\NRBP1_replica2\results_nrbp1pf2\smic_rmsd_rmsf_outputs\nrbp1_replica2"),
        preferred_prefix="nrbp1_replica2",
        topology_path=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\NRBP1_replica2\results_nrbp1pf2\smic_rmsd_rmsf_outputs\nrbp1_replica2\rmsd\_prepared_topology\replica_1_equilibrated__nocryst.pdb"),
        trajectory_paths=(
            Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\NRBP1_replica2\results_nrbp1pf2\replica_1-001 (1).dcd"),
        ),
    ),
    ReplicaSpec(
        system_label="FERRY CCTL",
        replica_label="R1",
        source_root=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\novel_ccts_and_WNK4(erendira)\novel_ccts\results_bound\C12OSR4_R1"),
        preferred_prefix="c12osr4_r1",
        topology_path=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\novel_ccts_and_WNK4(erendira)\novel_ccts\results_bound\C12OSR4_R1\smic_rmsd_rmsf_outputs_20260224\c12osr4_r1_20260224\rmsd\_prepared_topology\replica_1_equilibrated__C12ORF4__sanres_nocryst.pdb"),
        trajectory_paths=(
            Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\novel_ccts_and_WNK4(erendira)\novel_ccts\results_bound\C12OSR4_R1\replica_1__C12ORF4.dcd"),
        ),
    ),
    ReplicaSpec(
        system_label="FERRY CCTL",
        replica_label="R2",
        source_root=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\novel_ccts_and_WNK4(erendira)\novel_ccts\results_bound\C12OSR4_R2"),
        preferred_prefix="c12osr4_r2",
        topology_path=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\novel_ccts_and_WNK4(erendira)\novel_ccts\results_bound\C12OSR4_R2\smic_rmsd_rmsf_outputs_20260224\c12osr4_r2_20260224\rmsd\_prepared_topology\replica_2_equilibrated__C12ORF4__sanres_nocryst.pdb"),
        trajectory_paths=(
            Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\novel_ccts_and_WNK4(erendira)\novel_ccts\results_bound\C12OSR4_R2\replica_2_C12ORF4.dcd"),
        ),
    ),
    ReplicaSpec(
        system_label="WNK4 CCTL2",
        replica_label="R1",
        source_root=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\wnk4ROUND2_run\runs\replica_1\smic_rmsd_rmsf_outputs_20260223\wnk4round2_r1_20260223"),
        preferred_prefix="wnk4round2_r1_20260223",
        topology_path=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\wnk4ROUND2_run\runs\replica_1\smic_rmsd_rmsf_outputs_20260223\wnk4round2_r1_20260223\rmsd\_prepared_topology\replica_1_equilibrated__sanres_nocryst.pdb"),
        trajectory_paths=(
            Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\wnk4ROUND2_run\runs\replica_1\replica_1-003.dcd"),
        ),
    ),
    ReplicaSpec(
        system_label="WNK4 CCTL2",
        replica_label="R2",
        source_root=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\wnk4ROUND2_run\runs\replica_2\smic_rmsd_rmsf_outputs_20260223\wnk4round2_r2_20260223"),
        preferred_prefix="wnk4round2_r2_20260223",
        topology_path=Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\wnk4ROUND2_run\runs\replica_2\smic_rmsd_rmsf_outputs_20260223\wnk4round2_r2_20260223\rmsd\_prepared_topology\replica_2_equilibrated__sanres_nocryst.pdb"),
        trajectory_paths=(
            Path(r"d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\wnk4ROUND2_run\runs\replica_2\replica_2-002.dcd"),
        ),
    ),
)


def _slugify(value: str) -> str:
    return "_".join(value.lower().replace("/", " ").split())


def _moving_average(values: Sequence[float], window: int = 31) -> List[float]:
    series = pd.Series(values, dtype=float)
    return series.rolling(window=window, center=True, min_periods=1).mean().tolist()


def _rgba(hex_color: str, alpha: float) -> str:
    hex_color = hex_color.lstrip("#")
    red = int(hex_color[0:2], 16)
    green = int(hex_color[2:4], 16)
    blue = int(hex_color[4:6], 16)
    return f"rgba({red},{green},{blue},{alpha})"


def _find_csv(data_root: Path, suffix: str, preferred_prefix: Optional[str]) -> Path:
    direct_data = data_root / "rmsd" / "data"
    if direct_data.exists():
        candidates = sorted(direct_data.glob(f"*_{suffix}.csv"))
    else:
        candidates = sorted(data_root.rglob(f"*_{suffix}.csv"))

    if preferred_prefix:
        for candidate in candidates:
            if candidate.name.lower().startswith(preferred_prefix.lower()):
                return candidate

    if not candidates:
        raise FileNotFoundError(f"No {suffix} CSV found under {data_root}")

    canonical = [candidate for candidate in candidates if "patch" not in candidate.name.lower()]
    return canonical[0] if canonical else candidates[0]


def _resolve_time_ns(rmsd_df: pd.DataFrame) -> pd.Series:
    if "time_ns" in rmsd_df.columns:
        return pd.to_numeric(rmsd_df["time_ns"], errors="coerce")
    if "time_ps" in rmsd_df.columns:
        return pd.to_numeric(rmsd_df["time_ps"], errors="coerce") / 1000.0
    raise KeyError("RMSD CSV lacks time_ns/time_ps columns")


def _linear_time_axis(length: int, end_time_ns: float) -> List[float]:
    if length <= 1:
        return [0.0] * max(length, 1)
    return [index * end_time_ns / (length - 1) for index in range(length)]


def _prepare_rmsd_time_axis(spec: ReplicaSpec, rmsd_df: pd.DataFrame) -> Tuple[pd.DataFrame, str, float]:
    work = rmsd_df.copy()
    time_ns = _resolve_time_ns(work)
    native_end_ns = float(time_ns.max())

    if spec.runtime_end_ns is None:
        work["_plot_time_ns"] = time_ns
        return work, "native_rmsd_csv_time_axis", native_end_ns

    authoritative_end_ns = float(spec.runtime_end_ns)
    if native_end_ns <= 0:
        work["_plot_time_ns"] = pd.Series(_linear_time_axis(len(work), authoritative_end_ns), index=work.index)
        return work, "authoritative_runtime_linearized_missing_native_axis", authoritative_end_ns

    scale = authoritative_end_ns / native_end_ns
    work["_plot_time_ns"] = time_ns * scale
    if abs(scale - 1.0) <= 0.01:
        return work, "authoritative_runtime_matches_native_axis", authoritative_end_ns
    return work, "authoritative_runtime_rescaled_from_derived_axis", authoritative_end_ns


def _is_protein_residue(residue: Any) -> bool:
    aa3 = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "ASX", "GLX", "SEC", "PYL",
    }
    solvent_like = {
        "HOH", "WAT", "TIP3", "SOL", "H2O",
        "NA", "NA+", "CL", "CL-", "K", "K+", "CA", "CA2+", "MG", "MG2+",
    }
    name = (getattr(residue, "name", "") or "").strip().upper()
    if not name or name in solvent_like:
        return False
    flag = getattr(residue, "is_protein", None)
    if flag is not None:
        return bool(flag)
    return name in aa3


def _infer_rmsf_chain_roles(topology: md.Topology) -> Tuple[List[int], List[int]]:
    chain_stats: List[Tuple[int, int, int]] = []
    for chain in topology.chains:
        residues = list(chain.residues)
        protein_res = sum(1 for residue in residues if _is_protein_residue(residue))
        if protein_res <= 0:
            continue
        chain_stats.append((chain.index, protein_res, len(residues)))

    chain_stats.sort(key=lambda item: (-item[1], -item[2], item[0]))
    if not chain_stats:
        return [], []

    protein_chains = [chain_stats[0][0]]

    peptide_chains = [
        chain_id
        for chain_id, protein_res, _ in chain_stats
        if 2 <= protein_res <= 30 and chain_id not in protein_chains
    ]
    if not peptide_chains:
        peptide_chains = [
            chain_id for chain_id, _, _ in chain_stats if chain_id not in protein_chains
        ][:1]

    return protein_chains, peptide_chains


def _default_rmsf_chain_groups(topology: md.Topology) -> List[RmsfChainGroup]:
    protein_chains, peptide_chains = _infer_rmsf_chain_roles(topology)
    groups: List[RmsfChainGroup] = []
    if protein_chains:
        groups.append(RmsfChainGroup(tuple(protein_chains), "protein_ca", "Protein Cα", GREEN_DARK))
    if peptide_chains:
        groups.append(RmsfChainGroup(tuple(peptide_chains), "peptide_ca", "Peptide Cα", BLUE_DARK))
    return groups


def _select_chain_ca_atoms(topology: md.Topology, chain_ids: Sequence[int]) -> np.ndarray:
    selected: List[int] = []
    chain_id_set = set(chain_ids)
    for atom in topology.atoms:
        if atom.residue.chain.index not in chain_id_set:
            continue
        if atom.name != "CA":
            continue
        if not _is_protein_residue(atom.residue):
            continue
        selected.append(atom.index)
    return np.array(selected, dtype=int)


def _count_dcd_frames(path: Path) -> int:
    with md.formats.DCDTrajectoryFile(str(path), "r") as handle:
        return int(handle.__len__())


def _compute_self_aligned_ca_rmsf(
    trajectory_paths: Sequence[Path],
    topology_path: Path,
    atom_indices: np.ndarray,
    *,
    target_frames: int = 900,
) -> Tuple[np.ndarray, int, int]:
    if atom_indices.size == 0:
        return np.array([], dtype=float), 1, 0

    total_frames = sum(_count_dcd_frames(path) for path in trajectory_paths)
    stride = max(1, (total_frames + target_frames - 1) // target_frames)

    reference: Optional[md.Trajectory] = None
    blocks: List[Any] = []
    sampled_frames = 0
    for trajectory_path in trajectory_paths:
        for chunk in md.iterload(str(trajectory_path), top=str(topology_path), atom_indices=atom_indices, stride=stride, chunk=250):
            if reference is None:
                reference = chunk[0]
            chunk.superpose(reference, frame=0)
            blocks.append(chunk.xyz.copy())
            sampled_frames += chunk.n_frames

    if not blocks:
        return np.array([], dtype=float), stride, sampled_frames

    xyz = np.concatenate(blocks, axis=0)
    mean_xyz = xyz.mean(axis=0)
    diff = xyz - mean_xyz
    rmsf_nm = (np.sqrt(np.mean(np.sum(diff ** 2, axis=2), axis=0))).astype(float)
    return rmsf_nm, stride, sampled_frames


def _build_corrected_rmsf_dataframe(spec: ReplicaSpec) -> Tuple[pd.DataFrame, str, Dict[str, Any]]:
    if spec.topology_path is None or not spec.trajectory_paths:
        raise RuntimeError(f"Missing raw topology/trajectory mapping for {spec.system_label} {spec.replica_label}")

    frame0 = md.load_frame(str(spec.trajectory_paths[0]), 0, top=str(spec.topology_path))
    topology = frame0.topology
    chain_groups = list(spec.rmsf_chain_groups) or _default_rmsf_chain_groups(topology)

    rows: List[Dict[str, Any]] = []
    stride = 1
    sampled_frames = 0
    chain_group_meta: List[Dict[str, Any]] = []
    for panel_order, group in enumerate(chain_groups, start=1):
        ca_atoms = _select_chain_ca_atoms(topology, group.chain_ids)
        rmsf_nm, group_stride, group_frames = _compute_self_aligned_ca_rmsf(
            spec.trajectory_paths,
            spec.topology_path,
            ca_atoms,
        )
        stride = max(stride, group_stride)
        sampled_frames = max(sampled_frames, group_frames)
        chain_group_meta.append(
            {
                "chain_key": group.chain_key,
                "display_name": group.display_name,
                "chain_ids": list(group.chain_ids),
            }
        )

        for atom_idx, value in zip(ca_atoms, rmsf_nm):
            residue = topology.atom(int(atom_idx)).residue
            rows.append(
                {
                    "chain": group.chain_key,
                    "chain_display": group.display_name,
                    "panel_order": panel_order,
                    "color_hex": group.color,
                    "residue": int(residue.resSeq),
                    "residue_name": residue.name,
                    "residue_label": f"{residue.name}{int(residue.resSeq)}",
                    "rmsf_nm": float(value),
                    "rmsf_A": float(value) * 10.0,
                }
            )

    metadata = {
        "chain_groups": chain_group_meta,
        "rmsf_sampling_stride": stride,
        "rmsf_sampled_frames": sampled_frames,
    }
    return pd.DataFrame(rows), "corrected_raw_self_aligned_ca", metadata


def _rmsf_ca_views(rmsf_df: pd.DataFrame) -> Tuple[pd.DataFrame, str]:
    if {"panel_order", "chain_display"}.issubset(rmsf_df.columns):
        filtered = rmsf_df.copy()
        filtered["panel_order"] = pd.to_numeric(filtered["panel_order"], errors="coerce").fillna(999).astype(int)
        return filtered, "explicit_panel_rows"

    chain_values = set(rmsf_df["chain"].astype(str).str.lower().unique())
    if {"protein_ca", "peptide_ca"}.intersection(chain_values):
        filtered = rmsf_df[
            rmsf_df["chain"].astype(str).str.lower().isin({"protein_ca", "peptide_ca"})
        ].copy()
        return filtered, "explicit_ca_rows"

    filtered = rmsf_df[
        rmsf_df["chain"].astype(str).str.lower().isin({"protein", "peptide"})
    ].copy()
    return filtered, "legacy_chain_rows"


def _system_title(system_label: str, replica_label: str) -> str:
    return f"{system_label} {replica_label}"


def _plot_rmsd(system_title: str, rmsd_df: pd.DataFrame) -> go.Figure:
    time_ns = pd.to_numeric(rmsd_df["_plot_time_ns"], errors="coerce")

    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.10,
        subplot_titles=(
            "Protein Stability (Cα + Backbone)",
            "Peptide Global Structure (Backbone vs All Atoms)",
        ),
    )

    series_specs = (
        ("rmsd_protein_backbone_A", "Protein backbone", GREEN_DARK, 1, 1.8),
        ("rmsd_protein_ca_A", "Protein Cα", GREEN_LIGHT, 1, 2.6),
        ("rmsd_peptide_backbone_A", "Peptide backbone", BLUE_DARK, 2, 1.8),
        ("rmsd_peptide_all_A", "Peptide all atoms", BLUE_LIGHT, 2, 2.6),
    )

    for column, trace_name, color, row, width in series_specs:
        if column not in rmsd_df.columns:
            continue
        values = pd.to_numeric(rmsd_df[column], errors="coerce")
        fig.add_trace(
            go.Scatter(
                x=time_ns,
                y=values,
                mode="lines",
                name=trace_name,
                line=dict(color=color, width=width),
                hovertemplate="t = %{x:.2f} ns<br>RMSD = %{y:.2f} Å<extra></extra>",
            ),
            row=row,
            col=1,
        )
        fig.add_trace(
            go.Scatter(
                x=time_ns,
                y=_moving_average(values.ffill().bfill(), 41),
                mode="lines",
                showlegend=False,
                line=dict(color=color, width=4),
                opacity=0.22,
                hoverinfo="skip",
            ),
            row=row,
            col=1,
        )

    for row, ref_col, color in ((1, "rmsd_protein_ca_A", GREEN_LIGHT), (2, "rmsd_peptide_all_A", BLUE_LIGHT)):
        if ref_col not in rmsd_df.columns:
            continue
        values = pd.to_numeric(rmsd_df[ref_col], errors="coerce")
        mu = float(values.mean())
        sigma = float(values.std())
        fig.add_hrect(
            y0=mu - sigma,
            y1=mu + sigma,
            fillcolor=_rgba(color, 0.12),
            line_width=0,
            row=row,
            col=1,
        )
        fig.add_hline(
            y=mu,
            line_dash="dash",
            line_color=color,
            annotation_text=f"μ={mu:.2f} Å",
            annotation_position="top left",
            row=row,
            col=1,
        )

    time_end = float(time_ns.max())
    fig.update_layout(
        template="plotly_white",
        width=1800,
        height=1000,
        title=dict(
            text=f"<b>{system_title}</b> · RMSD over real simulation time<br><sub>Window covered: 0-{time_end:.1f} ns</sub>",
            x=0.5,
            xanchor="center",
        ),
        margin=dict(t=110, l=70, r=50, b=70),
        legend=dict(
            orientation="h",
            x=0.5,
            xanchor="center",
            y=1.08,
            yanchor="bottom",
            bgcolor="rgba(255,255,255,0.85)",
        ),
        paper_bgcolor="white",
        plot_bgcolor="white",
    )
    fig.update_xaxes(title_text="Time (ns)", row=2, col=1, gridcolor=GRID_COLOR)
    fig.update_yaxes(title_text="RMSD (Å)", row=1, col=1, gridcolor=GRID_COLOR)
    fig.update_yaxes(title_text="RMSD (Å)", row=2, col=1, gridcolor=GRID_COLOR)
    return fig


def _plot_rmsf(system_title: str, rmsf_df: pd.DataFrame, rmsf_mode: str) -> go.Figure:
    filtered, mode = _rmsf_ca_views(rmsf_df)
    panels: List[Tuple[str, str, str]] = []
    if mode == "explicit_panel_rows":
        panel_frame = (
            filtered[["chain", "chain_display", "panel_order", "color_hex"]]
            .drop_duplicates()
            .sort_values(["panel_order", "chain_display"])
        )
        panels = [
            (str(row.chain), str(row.chain_display), str(row.color_hex))
            for row in panel_frame.itertuples(index=False)
        ]
    else:
        for chain_name, display_name, color in (
            ("protein_ca", "Protein Cα", GREEN_DARK),
            ("peptide_ca", "Peptide Cα", BLUE_DARK),
            ("protein", "Protein", GREEN_DARK),
            ("peptide", "Peptide", BLUE_DARK),
        ):
            subset = filtered[filtered["chain"].astype(str).str.lower() == chain_name]
            if not subset.empty:
                panels.append((chain_name, display_name, color))

    if not panels:
        raise RuntimeError(f"No RMSF rows available for plotting in mode {rmsf_mode}")

    fig = make_subplots(
        rows=1,
        cols=len(panels),
        horizontal_spacing=0.08,
        subplot_titles=tuple(f"{display_name} RMSF" for _, display_name, _ in panels),
    )

    for col, (chain_name, display_name, color) in enumerate(panels, start=1):
        subset = filtered[filtered["chain"].astype(str).str.lower() == chain_name].copy()
        subset["residue"] = pd.to_numeric(subset["residue"], errors="coerce")
        subset["rmsf_A"] = pd.to_numeric(subset["rmsf_A"], errors="coerce")
        subset = subset.dropna(subset=["residue", "rmsf_A"]).sort_values("residue")
        labels = (
            subset["residue_label"]
            if "residue_label" in subset.columns
            else subset["residue_name"].astype(str) + subset["residue"].astype(int).astype(str)
        )

        fig.add_trace(
            go.Bar(
                x=subset["residue"],
                y=subset["rmsf_A"],
                marker=dict(color=color, line=dict(color=color, width=0.6)),
                opacity=0.88,
                customdata=list(labels),
                name=display_name,
                hovertemplate="%{customdata}<br>RMSF = %{y:.2f} Å<extra></extra>",
                showlegend=False,
            ),
            row=1,
            col=col,
        )
        fig.add_trace(
            go.Scatter(
                x=subset["residue"],
                y=_moving_average(subset["rmsf_A"], 11),
                mode="lines",
                line=dict(color=color, width=3),
                opacity=0.35,
                hoverinfo="skip",
                showlegend=False,
            ),
            row=1,
            col=col,
        )
        mu = float(subset["rmsf_A"].mean())
        sigma = float(subset["rmsf_A"].std())
        fig.add_hrect(
            y0=mu - sigma,
            y1=mu + sigma,
            fillcolor=_rgba(color, 0.12),
            line_width=0,
            row=1,
            col=col,
        )
        fig.add_hline(
            y=mu,
            line_dash="dash",
            line_color=color,
            annotation_text=f"μ={mu:.2f} Å",
            annotation_position="top left",
            row=1,
            col=col,
        )
        fig.update_xaxes(title_text="Residue", gridcolor=GRID_COLOR, row=1, col=col)
        fig.update_yaxes(title_text="RMSF (Å)", gridcolor=GRID_COLOR, row=1, col=col)

    if rmsf_mode == "corrected_raw_self_aligned_ca":
        mode_text = "corrected self-aligned Cα RMSF from raw DCD topology"
    elif mode == "explicit_panel_rows":
        mode_text = "explicit chain-aware Cα RMSF panels"
    elif mode == "explicit_ca_rows":
        mode_text = "explicit Cα rows"
    else:
        mode_text = "legacy RMSF rows used as best available precomputed residue RMSF"
    fig.update_layout(
        template="plotly_white",
        width=900 * len(panels),
        height=780,
        title=dict(
            text=f"<b>{system_title}</b> · RMSF (Cα-focused replot)<br><sub>Source mode: {mode_text}</sub>",
            x=0.5,
            xanchor="center",
        ),
        margin=dict(t=110, l=70, r=40, b=80),
        paper_bgcolor="white",
        plot_bgcolor="white",
    )
    return fig


def _write_combined_html(system_title: str, rmsd_fig: go.Figure, rmsf_fig: go.Figure, output_html: Path) -> None:
    rmsd_html = to_html(rmsd_fig, include_plotlyjs="cdn", full_html=False, div_id=f"{_slugify(system_title)}_rmsd")
    rmsf_html = to_html(rmsf_fig, include_plotlyjs=False, full_html=False, div_id=f"{_slugify(system_title)}_rmsf")
    html = f"""
<!DOCTYPE html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\" />
  <title>{system_title} · RMSD/RMSF</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 0; background: #f5f8fb; color: #10324a; }}
    main {{ max-width: 1880px; margin: 0 auto; padding: 28px 24px 40px; }}
    h1 {{ margin: 0 0 10px; font-size: 34px; }}
    p {{ margin: 0 0 20px; font-size: 16px; line-height: 1.5; }}
    section {{ background: white; padding: 18px; border-radius: 18px; box-shadow: 0 10px 35px rgba(16, 50, 74, 0.08); margin-bottom: 22px; }}
  </style>
</head>
<body>
  <main>
    <h1>{system_title}</h1>
    <p>Replot generado desde CSV precomputados de SMIC con RMSD en tiempo real y RMSF Cα-focused.</p>
    <section>{rmsd_html}</section>
    <section>{rmsf_html}</section>
  </main>
</body>
</html>
"""
    output_html.write_text(html, encoding="utf-8")


def _write_rmsd_png(system_title: str, rmsd_df: pd.DataFrame, output_png: Path) -> None:
    time_ns = pd.to_numeric(rmsd_df["_plot_time_ns"], errors="coerce")
    fig, axes = plt.subplots(2, 1, figsize=(16, 9), sharex=True, constrained_layout=True)
    fig.patch.set_facecolor("white")

    panel_specs = (
        (axes[0], "rmsd_protein_backbone_A", "Protein backbone", GREEN_DARK, "rmsd_protein_ca_A", "Protein Cα", GREEN_LIGHT),
        (axes[1], "rmsd_peptide_backbone_A", "Peptide backbone", BLUE_DARK, "rmsd_peptide_all_A", "Peptide all atoms", BLUE_LIGHT),
    )

    for axis, col1, label1, color1, col2, label2, color2 in panel_specs:
        if col1 in rmsd_df.columns:
            values1 = pd.to_numeric(rmsd_df[col1], errors="coerce")
            axis.plot(time_ns, values1, color=color1, linewidth=1.8, label=label1)
            axis.plot(time_ns, _moving_average(values1.ffill().bfill(), 41), color=color1, linewidth=3.2, alpha=0.22)
        if col2 in rmsd_df.columns:
            values2 = pd.to_numeric(rmsd_df[col2], errors="coerce")
            axis.plot(time_ns, values2, color=color2, linewidth=2.4, label=label2)
            axis.plot(time_ns, _moving_average(values2.ffill().bfill(), 41), color=color2, linewidth=3.6, alpha=0.22)
            mu = float(values2.mean())
            sigma = float(values2.std())
            axis.axhspan(mu - sigma, mu + sigma, color=color2, alpha=0.10)
            axis.axhline(mu, color=color2, linestyle="--", linewidth=1.2)
        axis.set_ylabel("RMSD (Å)")
        axis.grid(True, color="#D9E5EE", linewidth=0.8)
        axis.legend(frameon=False, loc="upper right")

    axes[0].set_title(system_title + " · RMSD over real simulation time", fontsize=17, pad=14)
    axes[1].set_xlabel("Time (ns)")
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def _write_rmsf_png(system_title: str, rmsf_df: pd.DataFrame, output_png: Path) -> None:
    filtered, mode = _rmsf_ca_views(rmsf_df)
    panel_specs: List[Tuple[pd.DataFrame, str, str]] = []
    if mode == "explicit_panel_rows":
        panel_frame = (
            filtered[["chain", "chain_display", "panel_order", "color_hex"]]
            .drop_duplicates()
            .sort_values(["panel_order", "chain_display"])
        )
        for row in panel_frame.itertuples(index=False):
            subset = filtered[filtered["chain"].astype(str) == str(row.chain)].copy()
            if not subset.empty:
                panel_specs.append((subset, str(row.chain_display), str(row.color_hex)))
    else:
        for chain_name, display_name, color in (
            ("protein_ca", "Protein Cα", GREEN_DARK),
            ("peptide_ca", "Peptide Cα", BLUE_DARK),
            ("protein", "Protein", GREEN_DARK),
            ("peptide", "Peptide", BLUE_DARK),
        ):
            subset = filtered[filtered["chain"].astype(str).str.lower() == chain_name].copy()
            if not subset.empty:
                panel_specs.append((subset, display_name, color))

    fig, axes = plt.subplots(1, len(panel_specs), figsize=(8 * len(panel_specs), 6), constrained_layout=True)
    fig.patch.set_facecolor("white")
    if len(panel_specs) == 1:
        axes = [axes]

    for axis, (subset, display_name, color) in zip(axes, panel_specs):
        subset["residue"] = pd.to_numeric(subset["residue"], errors="coerce")
        subset["rmsf_A"] = pd.to_numeric(subset["rmsf_A"], errors="coerce")
        subset = subset.dropna(subset=["residue", "rmsf_A"]).sort_values("residue")
        axis.bar(subset["residue"], subset["rmsf_A"], color=color, alpha=0.85, width=0.9)
        axis.plot(subset["residue"], _moving_average(subset["rmsf_A"], 11), color=color, linewidth=2.8, alpha=0.38)
        mu = float(subset["rmsf_A"].mean())
        sigma = float(subset["rmsf_A"].std())
        axis.axhspan(mu - sigma, mu + sigma, color=color, alpha=0.10)
        axis.axhline(mu, color=color, linestyle="--", linewidth=1.2)
        axis.set_title(display_name)
        axis.set_xlabel("Residue")
        axis.set_ylabel("RMSF (Å)")
        axis.grid(True, color="#D9E5EE", linewidth=0.8)

    title_suffix = "explicit Cα rows" if mode == "explicit_ca_rows" else "corrected or fallback RMSF rows"
    fig.suptitle(system_title + f" · RMSF (Cα-focused replot, {title_suffix})", fontsize=16)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def _dashboard_html(records: Sequence[Dict[str, Any]], output_html: Path) -> None:
    cards = []
    for record in records:
        rel_html = Path(record["combined_html_name"]).name
        rel_rmsd_png = Path(record["rmsd_png_name"]).name
        rel_rmsf_png = Path(record["rmsf_png_name"]).name
        cards.append(
            f"""
            <article>
              <h2>{record['system_title']}</h2>
              <p><strong>Input RMSD:</strong> {record['rmsd_csv']}</p>
              <p><strong>Input RMSF:</strong> {record['rmsf_csv']}</p>
              <p><strong>RMSF mode:</strong> {record['rmsf_mode']}</p>
              <p><strong>Runtime source:</strong> {record['runtime_source']}</p>
              <p><strong>Duration:</strong> {record['time_end_ns']:.2f} ns</p>
              <p><a href=\"{rel_html}\">Combined interactive HTML</a></p>
              <div class=\"grid\">
                <figure>
                  <img src=\"{rel_rmsd_png}\" alt=\"{record['system_title']} RMSD\" />
                  <figcaption>RMSD</figcaption>
                </figure>
                <figure>
                  <img src=\"{rel_rmsf_png}\" alt=\"{record['system_title']} RMSF\" />
                  <figcaption>RMSF</figcaption>
                </figure>
              </div>
            </article>
            """
        )

    html = f"""
<!DOCTYPE html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\" />
  <title>Control MD RMSD/RMSF Dashboard</title>
  <style>
    body {{ font-family: Arial, sans-serif; background: linear-gradient(180deg, #eef6fb 0%, #f8fbfd 100%); color: #14364c; margin: 0; }}
    main {{ max-width: 1900px; margin: 0 auto; padding: 28px 24px 40px; }}
    h1 {{ margin: 0 0 10px; font-size: 38px; }}
    p {{ line-height: 1.55; }}
    article {{ background: white; border-radius: 20px; padding: 22px; box-shadow: 0 12px 34px rgba(20, 54, 76, 0.08); margin-top: 22px; }}
    .grid {{ display: grid; grid-template-columns: repeat(2, minmax(280px, 1fr)); gap: 16px; margin-top: 16px; }}
    figure {{ margin: 0; background: #f7fafc; border-radius: 14px; padding: 12px; }}
    img {{ width: 100%; height: auto; border-radius: 10px; display: block; }}
    figcaption {{ padding-top: 10px; font-size: 14px; }}
    a {{ color: #1e6ea1; }}
  </style>
</head>
<body>
  <main>
    <h1>Unified MD RMSD/RMSF Replot Dashboard</h1>
    <p>Research-grade replot of precomputed SMIC RMSD plus regenerated self-aligned Cα RMSF from raw trajectories for all requested control systems and replicas. RMSD uses real simulation time on the x-axis.</p>
    {''.join(cards)}
  </main>
</body>
</html>
"""
    output_html.write_text(html, encoding="utf-8")


def _process_replica(spec: ReplicaSpec, output_dir: Path) -> Dict[str, Any]:
    rmsd_csv = _find_csv(spec.source_root, "rmsd", spec.preferred_prefix)
    rmsf_csv = _find_csv(spec.source_root, "rmsf", spec.preferred_prefix)

    rmsd_df = pd.read_csv(rmsd_csv)
    rmsf_df, rmsf_mode, rmsf_meta = _build_corrected_rmsf_dataframe(spec)
    rmsd_window_df, time_mode, time_end_ns = _prepare_rmsd_time_axis(spec, rmsd_df)

    system_title = _system_title(spec.system_label, spec.replica_label)
    base_name = _slugify(system_title)

    rmsd_fig = _plot_rmsd(system_title, rmsd_window_df)
    rmsf_fig = _plot_rmsf(system_title, rmsf_df, rmsf_mode)

    rmsd_png = output_dir / f"{base_name}_rmsd.png"
    rmsf_png = output_dir / f"{base_name}_rmsf_ca.png"
    combined_html = output_dir / f"{base_name}_rmsd_rmsf.html"
    rmsd_html = output_dir / f"{base_name}_rmsd.html"
    rmsf_html = output_dir / f"{base_name}_rmsf_ca.html"
    rmsf_data_csv = output_dir / f"{base_name}_rmsf_ca_corrected.csv"

    rmsd_fig.write_html(str(rmsd_html))
    rmsf_fig.write_html(str(rmsf_html))
    _write_combined_html(system_title, rmsd_fig, rmsf_fig, combined_html)
    _write_rmsd_png(system_title, rmsd_window_df, rmsd_png)
    _write_rmsf_png(system_title, rmsf_df, rmsf_png)
    rmsf_df.to_csv(rmsf_data_csv, index=False)

    return {
        "system_label": spec.system_label,
        "replica_label": spec.replica_label,
        "system_title": system_title,
        "source_root": str(spec.source_root),
        "rmsd_csv": str(rmsd_csv),
        "rmsf_csv": str(rmsf_csv),
        "rmsf_mode": rmsf_mode,
        "rmsf_runtime_source": "raw_dcd_self_aligned_ca",
        "rmsf_sampling_stride": rmsf_meta["rmsf_sampling_stride"],
        "rmsf_sampled_frames": rmsf_meta["rmsf_sampled_frames"],
        "time_mode": time_mode,
        "runtime_source": spec.runtime_source,
        "time_end_ns": time_end_ns,
        "rmsd_html_name": rmsd_html.name,
        "rmsf_html_name": rmsf_html.name,
        "combined_html_name": combined_html.name,
        "rmsd_png_name": rmsd_png.name,
        "rmsf_png_name": rmsf_png.name,
        "rmsf_data_csv_name": rmsf_data_csv.name,
    }


def _write_summary(records: Sequence[Dict[str, Any]], output_dir: Path) -> None:
    summary_df = pd.DataFrame(records)
    summary_df.to_csv(output_dir / "replot_manifest.csv", index=False)

    (output_dir / "replot_manifest.json").write_text(
        json.dumps(records, indent=2),
        encoding="utf-8",
    )

    lines = [
        "# Unified MD RMSD/RMSF Replot Manifest",
        "",
        "This document records each resolved source CSV pair and the generated artifacts.",
        "",
    ]
    for record in records:
        lines.extend(
            [
                f"## {record['system_title']}",
                "",
                f"- Source root: {record['source_root']}",
                f"- RMSD CSV: {record['rmsd_csv']}",
                f"- RMSF CSV: {record['rmsf_csv']}",
                f"- RMSF mode: {record['rmsf_mode']}",
                f"- Corrected RMSF CSV: {record['rmsf_data_csv_name']}",
                f"- RMSF source: {record['rmsf_runtime_source']}",
                f"- RMSF sampling stride: {record['rmsf_sampling_stride']}",
                f"- RMSF sampled frames: {record['rmsf_sampled_frames']}",
                f"- Time mode: {record['time_mode']}",
                f"- Runtime source: {record['runtime_source']}",
                f"- Time covered: {record['time_end_ns']:.2f} ns",
                f"- Combined HTML: {record['combined_html_name']}",
                f"- RMSD PNG: {record['rmsd_png_name']}",
                f"- RMSF PNG: {record['rmsf_png_name']}",
                "",
            ]
        )
    (output_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    records = [_process_replica(spec, OUTPUT_ROOT) for spec in REPLICA_SPECS]
    _write_summary(records, OUTPUT_ROOT)
    _dashboard_html(records, OUTPUT_ROOT / "all_systems_rmsd_rmsf_dashboard.html")
    print(f"[OK] Wrote unified replots to: {OUTPUT_ROOT}")


if __name__ == "__main__":
    main()
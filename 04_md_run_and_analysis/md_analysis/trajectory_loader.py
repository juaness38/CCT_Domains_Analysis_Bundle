"""
CCT-OSR1 Binding Analysis - Trajectory Loader
==============================================
Handles loading, validation, and concatenation of MD trajectories.
"""

import mdtraj as md
import numpy as np
from pathlib import Path
from typing import Iterable, List, Optional, Tuple, Sequence
import warnings
import os

try:
    from . import config as _config  # Legacy defaults (optional)
except ImportError:  # pragma: no cover - optional dependency
    try:  # pragma: no cover
        import config as _config
    except ImportError:  # pragma: no cover
        _config = None

# ---------------------------------------------------------------------------
# Legacy defaults (used when runtime configuration is not provided)
# ---------------------------------------------------------------------------
DEFAULT_TOPOLOGY = Path(getattr(_config, "TOPOLOGY_PDB", "")) if _config else None
DEFAULT_TRAJECTORIES = (
    [Path(p) for p in getattr(_config, "TRAJECTORIES", [])]
    if _config
    else []
)
DEFAULT_SAVE_FREQ_PS = getattr(_config, "SAVE_FREQ_PS", None)
DEFAULT_DATA_DIR = Path(getattr(_config, "DATA_OUTPUT_DIR", "data")) if _config else Path("outputs/data")

# Keep full in-memory trajectories at a conservative size to avoid OOMs when
# joining chunked loads. Can be overridden with env var SMIC_MAX_TRAJ_GB.
DEFAULT_MAX_TRAJ_GB = 2.0


def _get_max_traj_gb() -> float:
    raw = os.getenv("SMIC_MAX_TRAJ_GB", "").strip()
    if not raw:
        return DEFAULT_MAX_TRAJ_GB
    try:
        value = float(raw)
        return value if np.isfinite(value) and value > 0 else DEFAULT_MAX_TRAJ_GB
    except (TypeError, ValueError):
        return DEFAULT_MAX_TRAJ_GB


def _estimate_frame_count(traj_path: Path) -> Optional[int]:
    """Best-effort frame count from trajectory header without loading coordinates."""
    try:
        with md.open(str(traj_path)) as handle:
            count = len(handle)
            return int(count) if count and count > 0 else None
    except Exception:
        return None


def _estimate_topology_atom_count(topology: Optional[Path]) -> Optional[int]:
    if topology is None:
        return None
    try:
        top = md.load_topology(str(topology))
        n_atoms = int(getattr(top, "n_atoms", 0))
        return n_atoms if n_atoms > 0 else None
    except Exception:
        return None


def _adaptive_stride_for_memory(
    traj_path: Path,
    topology: Optional[Path],
    requested_stride: int,
) -> int:
    """Increase stride automatically when expected in-memory trajectory is too large."""
    requested_stride = max(1, int(requested_stride))
    n_frames = _estimate_frame_count(traj_path)
    n_atoms = _estimate_topology_atom_count(topology)
    if n_frames is None or n_atoms is None:
        return requested_stride

    # xyz uses float32 in MDTraj => 4 bytes per scalar
    bytes_per_frame = float(n_atoms) * 3.0 * 4.0
    target_bytes = float(_get_max_traj_gb()) * (1024.0 ** 3)
    if target_bytes <= 0:
        return requested_stride

    expected_stride1 = bytes_per_frame * float(n_frames)
    if expected_stride1 <= target_bytes:
        return requested_stride

    min_stride_from_memory = int(np.ceil(expected_stride1 / target_bytes))
    return max(requested_stride, min_stride_from_memory)


def _infer_saved_frame_dt_ps_with_mdanalysis(
    topology: Optional[Path],
    traj_path: Path,
) -> Optional[float]:
    """Best-effort: infer time between *saved* frames (ps) from trajectory metadata.

    This is intentionally lightweight: MDAnalysis reads the trajectory header to
    expose `trajectory.dt` without loading coordinates for all frames.
    """

    if topology is None:
        return None

    try:
        import MDAnalysis as mda  # type: ignore

        u = mda.Universe(str(topology), str(traj_path))
        dt = float(getattr(u.trajectory, "dt", 0.0))
        if np.isfinite(dt) and dt > 0:
            return dt
    except Exception:  # pragma: no cover
        return None

    return None


def _resolve_frame_dt_ps(
    traj: md.Trajectory,
    *,
    stride: int,
    override_time_step_ps: Optional[float] = None,
) -> Optional[float]:
    """Resolve the time between *loaded* frames (ps).

    Notes
    -----
    - `override_time_step_ps` (from CLI) and `DEFAULT_SAVE_FREQ_PS` represent the
      time between *saved* frames in the original trajectory. We scale by
      `stride` to get the time between frames *after* loading with stride.
    - When MDTraj provides a per-frame `traj.time` array (ps), we prefer it
      because it's the most direct representation of the data.
    """

    stride = max(1, int(stride))

    # CLI override (or MDAnalysis-inferred saved-frame dt): per-saved-frame (ps)
    # Always prefer this when provided because some DCDs yield incorrect time
    # arrays/timesteps via MDTraj.
    if override_time_step_ps is not None:
        return float(override_time_step_ps) * float(stride)

    # Prefer explicit per-frame time stamps if available.
    try:
        times = getattr(traj, "time", None)
        if times is not None and len(times) >= 2:
            dt = float(times[1] - times[0])
            if np.isfinite(dt) and dt > 0:
                return dt
    except Exception:  # pragma: no cover
        pass

    # MDTraj-reported timestep is assumed to be per-loaded-frame (ps)
    candidate = getattr(traj, "timestep", None)
    if candidate is not None:
        try:
            dt = float(candidate)
            if np.isfinite(dt) and dt > 0:
                return dt
        except (TypeError, ValueError):  # pragma: no cover - rare types
            pass

    # Legacy default: per-saved-frame (ps)
    if DEFAULT_SAVE_FREQ_PS is not None:
        return float(DEFAULT_SAVE_FREQ_PS) * float(stride)

    return None


def load_topology(pdb_path: Optional[Path] = None) -> md.Topology:
    """Load and validate topology from PDB file."""
    pdb_path = Path(pdb_path) if pdb_path else DEFAULT_TOPOLOGY
    if pdb_path is None or (isinstance(pdb_path, Path) and not str(pdb_path)):
        raise ValueError("Topology path must be provided")

    print(f"Loading topology from: {pdb_path}")
    traj = md.load(str(pdb_path))
    top = traj.topology
    
    # Validate structure
    n_chains = top.n_chains
    n_residues = top.n_residues
    n_atoms = top.n_atoms
    
    print(f"  Chains: {n_chains}")
    print(f"  Residues: {n_residues}")
    print(f"  Atoms: {n_atoms}")
    
    # List chains
    for chain in top.chains:
        residues = list(chain.residues)
        if len(residues) > 0:
            first_res = residues[0]
            last_res = residues[-1]
            print(f"  Chain {chain.index}: {first_res.name}{first_res.resSeq} - {last_res.name}{last_res.resSeq} ({len(residues)} residues)")
    
    return top


def load_single_trajectory(
    traj_path: Path,
    topology: Optional[Path] = None,
    stride: int = 1,
    chunk_size: int = 250,
    time_step_ps: Optional[float] = None,
) -> md.Trajectory:
    """Load a single trajectory file using chunked iterload to limit memory usage."""
    print(f"Loading trajectory: {traj_path.name}")

    if chunk_size <= 0:
        raise ValueError("chunk_size must be a positive integer")

    effective_stride = _adaptive_stride_for_memory(traj_path, Path(topology) if topology else DEFAULT_TOPOLOGY, stride)
    if effective_stride > max(1, int(stride)):
        n_frames = _estimate_frame_count(traj_path)
        n_atoms = _estimate_topology_atom_count(Path(topology) if topology else DEFAULT_TOPOLOGY)
        est_gib = None
        if n_frames is not None and n_atoms is not None:
            est_bytes = float(n_frames) * float(n_atoms) * 3.0 * 4.0 / float(effective_stride)
            est_gib = est_bytes / (1024.0 ** 3)
        target_gib = _get_max_traj_gb()
        if est_gib is not None:
            print(
                f"  [WARN] Trajectory too large for safe in-memory join with stride={stride}; "
                f"using adaptive stride={effective_stride} (~{est_gib:.2f} GiB target <= {target_gib:.2f} GiB)."
            )
        else:
            print(
                f"  [WARN] Trajectory too large for safe in-memory join with stride={stride}; "
                f"using adaptive stride={effective_stride}."
            )

    chunks: List[md.Trajectory] = []
    total_frames = 0

    for chunk_idx, chunk in enumerate(
        md.iterload(str(traj_path), top=str(topology or DEFAULT_TOPOLOGY), stride=effective_stride, chunk=chunk_size)
    ):
        chunks.append(chunk)
        total_frames += chunk.n_frames
        print(f"    chunk {chunk_idx + 1:02d}: {chunk.n_frames} frames")

    if not chunks:
        raise RuntimeError(f"No frames loaded from {traj_path}")

    traj = md.join(chunks) if len(chunks) > 1 else chunks[0]

    effective_time_step_ps = time_step_ps
    if effective_time_step_ps is None:
        inferred = _infer_saved_frame_dt_ps_with_mdanalysis(Path(topology or DEFAULT_TOPOLOGY), traj_path)
        if inferred is not None:
            effective_time_step_ps = inferred

    frame_dt_ps = _resolve_frame_dt_ps(traj, stride=effective_stride, override_time_step_ps=effective_time_step_ps)
    if frame_dt_ps is not None:
        try:
            traj.time = np.arange(traj.n_frames, dtype=float) * float(frame_dt_ps)
        except Exception:  # pragma: no cover
            pass
        time_ns = total_frames * frame_dt_ps / 1000
        print(f"  Frames: {total_frames} ({time_ns:.1f} ns)")
    else:
        print(f"  Frames: {total_frames}")

    return traj


def load_and_concatenate_trajectories(
    traj_paths: Optional[Iterable[Path]] = None,
    topology: Optional[Path] = None,
    stride: int = 1,
    validate: bool = True,
    time_step_ps: Optional[float] = None,
    box_lengths_nm: Optional[Sequence[float]] = None,
    autoimage_selection: Optional[str] = "chainid 0 and name CA",
    apply_autoimage: bool = True,
) -> md.Trajectory:
    """
    Load and concatenate multiple trajectory files.
    
    Parameters
    ----------
    traj_paths : list of Path
        Ordered list of trajectory files to concatenate
    topology : Path
        PDB file with topology
    stride : int
        Load every nth frame (for memory efficiency)
    validate : bool
        Check trajectory continuity
        
    Returns
    -------
    md.Trajectory
        Concatenated trajectory
    """
    print("=" * 60)
    print("LOADING TRAJECTORIES")
    print("=" * 60)
    
    trajectories = []
    total_frames = 0
    
    resolved_paths = list(traj_paths) if traj_paths is not None else list(DEFAULT_TRAJECTORIES)
    if not resolved_paths:
        raise ValueError("No trajectory files provided")

    topo_path = Path(topology) if topology else DEFAULT_TOPOLOGY
    box_tuple = tuple(float(v) for v in box_lengths_nm) if box_lengths_nm is not None else None

    expected_n_atoms: Optional[int] = None
    skipped: List[Tuple[Path, str]] = []

    # Track time offset (ps) so concatenated trajectories have a continuous time axis.
    time_offset_ps: float = 0.0

    for i, traj_path in enumerate(resolved_paths):
        if not traj_path.exists():
            raise FileNotFoundError(f"Trajectory not found: {traj_path}")

        try:
            traj = load_single_trajectory(traj_path, topo_path, stride, time_step_ps=time_step_ps)

            if expected_n_atoms is None:
                expected_n_atoms = int(traj.n_atoms)
            elif int(traj.n_atoms) != expected_n_atoms:
                skipped.append((traj_path, f"atom_count_mismatch expected={expected_n_atoms} got={traj.n_atoms}"))
                continue

            # Ensure each segment's time axis is monotonic and continuous.
            # MDTraj/MDAnalysis commonly start each file at t=0; for analysis/UI,
            # we want a single continuous axis across concatenated segments.
            try:
                times = getattr(traj, "time", None)
                if times is not None and len(times) == traj.n_frames and traj.n_frames > 0:
                    t = np.asarray(times, dtype=float)
                    if np.all(np.isfinite(t)):
                        dt = float(t[1] - t[0]) if len(t) > 1 else 0.0
                        t = (t - float(t[0])) + float(time_offset_ps)
                        traj.time = t
                        time_offset_ps = float(t[-1]) + (dt if (np.isfinite(dt) and dt > 0) else 0.0)
            except Exception:
                pass

            trajectories.append(traj)
            total_frames += traj.n_frames
        except Exception as exc:
            # Inventories often include mixed artifacts (stripped vs solvated,
            # partial copies, incompatible formats). Skip per-file and continue.
            skipped.append((traj_path, f"load_failed: {type(exc).__name__}: {exc}"))
            continue
    
    if not trajectories:
        details = "\n".join([f"- {p}: {reason}" for p, reason in skipped[:12]])
        more = "" if len(skipped) <= 12 else f"\n... ({len(skipped) - 12} more skipped)"
        raise RuntimeError(
            "No compatible trajectories could be loaded with the provided topology.\n"
            "Skipped files:\n" + details + more
        )

    if skipped:
        print(f"\n[WARN] Skipped {len(skipped)} incompatible trajectory file(s):")
        for p, reason in skipped[:12]:
            print(f"  - {p.name}: {reason}")
        if len(skipped) > 12:
            print(f"  ... ({len(skipped) - 12} more)")

    # Concatenate
    print(f"\nConcatenating {len(trajectories)} trajectories...")
    combined = md.join(trajectories)
    
    # Assign proper time
    def _best_loaded_dt_ps(trajs: Sequence[md.Trajectory]) -> Optional[float]:
        for candidate_traj in trajs:
            try:
                times = getattr(candidate_traj, "time", None)
                if times is not None and len(times) >= 2:
                    dt = float(times[1] - times[0])
                    if np.isfinite(dt) and dt > 0:
                        return dt
            except Exception:  # pragma: no cover
                pass
            try:
                ts = getattr(candidate_traj, "timestep", None)
                if ts is not None:
                    dt = float(ts)
                    if np.isfinite(dt) and dt > 0:
                        return dt
            except Exception:  # pragma: no cover
                pass
        return None

    frame_dt_ps = _best_loaded_dt_ps(trajectories)
    if frame_dt_ps is None:
        frame_dt_ps = _resolve_frame_dt_ps(combined, stride=stride, override_time_step_ps=time_step_ps)
    if frame_dt_ps is not None:
        time_ps = np.arange(combined.n_frames, dtype=float) * float(frame_dt_ps)
        combined.time = time_ps
        total_time_ns = combined.n_frames * float(frame_dt_ps) / 1000
        print(f"Combined trajectory: {combined.n_frames} frames ({total_time_ns:.1f} ns)")
    else:
        print(f"Combined trajectory: {combined.n_frames} frames")
    
    if validate:
        validate_trajectory(combined)

    if box_tuple is not None:
        _assign_box_vectors(combined, box_tuple)
        angstroms = tuple(val * 10.0 for val in box_tuple)
        print(f"Applying periodic box: {angstroms[0]:.2f} A x {angstroms[1]:.2f} A x {angstroms[2]:.2f} A")
        if apply_autoimage and autoimage_selection:
            _autoimage_with_anchor(combined, autoimage_selection, box_tuple)
    
    return combined


def validate_trajectory(traj: md.Trajectory) -> bool:
    """Validate trajectory for common issues."""
    print("\nValidating trajectory...")
    
    issues = []
    
    # Check for NaN positions
    chunk = 250
    has_nan = False
    nan_count = 0
    for start in range(0, traj.n_frames, chunk):
        stop = min(traj.n_frames, start + chunk)
        block = traj.xyz[start:stop]
        if np.isnan(block).any():
            has_nan = True
            nan_count += np.isnan(block).sum()
            break
    if has_nan:
        issues.append(f"NaN positions detected: {nan_count} values")
    
    # Check box vectors
    if traj.unitcell_vectors is not None:
        box_lengths = np.linalg.norm(traj.unitcell_vectors, axis=2)
        if np.any(box_lengths < 1.0) or np.any(box_lengths > 50.0):
            issues.append("Unusual box dimensions detected")
    
    # Check for exploded system (atoms too far from COM)
    com = md.compute_center_of_mass(traj)
    
    if len(issues) == 0:
        print("  [OK] Trajectory validation passed")
        return True
    else:
        for issue in issues:
            print(f"  [WARN] {issue}")
        return False


def _assign_box_vectors(traj: md.Trajectory, box_lengths_nm: Sequence[float]) -> None:
    """Populate unitcell vectors for a trajectory given box lengths in nm."""
    box = np.asarray(box_lengths_nm, dtype=float)
    if box.shape != (3,):
        raise ValueError("box_lengths_nm must contain three values")
    vectors = np.zeros((traj.n_frames, 3, 3), dtype=float)
    for dim in range(3):
        vectors[:, dim, dim] = box[dim]
    traj.unitcell_vectors = vectors


def _autoimage_with_anchor(
    traj: md.Trajectory,
    selection: str,
    box_lengths_nm: Sequence[float],
) -> None:
    """Remove PBC jumps by keeping an anchor selection continuous."""
    anchor_idx = traj.topology.select(selection)
    if anchor_idx.size == 0:
        print(f"  [WARN] Autoimage skipped: selection '{selection}' empty")
        return
    box = np.asarray(box_lengths_nm, dtype=float)
    coords = traj.xyz
    anchor_positions = coords[:, anchor_idx, :].mean(axis=1)
    prev = anchor_positions[0]
    half_box = box / 2.0

    for frame in range(1, traj.n_frames):
        current = anchor_positions[frame]
        delta = current - prev
        shift_units = np.round(delta / box)
        shift_units[np.abs(delta) < half_box] = 0.0
        shift = shift_units * box
        if np.any(shift_units != 0):
            coords[frame] -= shift
            anchor_positions[frame] = current - shift
            prev = anchor_positions[frame]
        else:
            prev = current

    traj.xyz = coords
    print(f"  ✓ Applied autoimage using '{selection}' as anchor")


def get_atom_indices(traj: md.Trajectory, selection: str) -> np.ndarray:
    """Get atom indices from MDTraj selection string."""
    return traj.topology.select(selection)


def get_chain_indices(traj: md.Trajectory) -> dict:
    """Get atom indices for each chain."""
    top = traj.topology
    indices = {}
    
    for chain in top.chains:
        chain_atoms = [atom.index for atom in chain.atoms]
        indices[chain.index] = np.array(chain_atoms)
    
    return indices


def get_protein_peptide_indices(traj: md.Trajectory) -> Tuple[np.ndarray, np.ndarray]:
    """Get atom indices for protein and peptide chains.

    Robust heuristic:
    - Identify "protein-like" chains using residue classification (exclude solvent/ions).
    - Pick the two chains with the most protein residues.
    - Designate the longer protein-like chain as protein, the shorter as peptide.

    This keeps downstream selections stable (protein=chainid 0, peptide=chainid 1)
    after slicing in :func:`extract_protein_peptide`.
    """
    top = traj.topology

    aa3 = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "ASX", "GLX", "SEC", "PYL",
    }
    solvent_like = {
        "HOH", "WAT", "TIP3", "SOL", "H2O",
        "NA", "NA+", "CL", "CL-", "K", "K+", "CA", "CA2+", "MG", "MG2+",
    }

    def _is_protein_residue(res) -> bool:
        name = (getattr(res, "name", "") or "").strip().upper()
        if not name or name in solvent_like:
            return False
        # mdtraj residues usually have `.is_protein`
        flag = getattr(res, "is_protein", None)
        if flag is not None:
            return bool(flag)
        return name in aa3

    chain_stats = []
    for chain in top.chains:
        residues = list(chain.residues)
        protein_res = sum(1 for r in residues if _is_protein_residue(r))
        if protein_res <= 0:
            continue
        chain_stats.append((chain, protein_res, len(residues)))

    # Fallback: keep legacy behavior if we can't confidently detect.
    if len(chain_stats) < 2:
        protein_atoms = []
        peptide_atoms = []
        for chain in top.chains:
            atoms = [atom.index for atom in chain.atoms]
            if chain.index == 0:
                protein_atoms = atoms
            elif chain.index == 1:
                peptide_atoms = atoms
        return np.array(protein_atoms), np.array(peptide_atoms)

    chain_stats.sort(key=lambda x: (x[1], x[2]), reverse=True)
    c_a, c_b = chain_stats[0], chain_stats[1]

    chain_a, a_prot_res, a_len = c_a
    chain_b, b_prot_res, b_len = c_b

    # Protein is the longer protein-like chain (by protein residue count, then total residues).
    if (a_prot_res, a_len) >= (b_prot_res, b_len):
        protein_chain = chain_a
        peptide_chain = chain_b
    else:
        protein_chain = chain_b
        peptide_chain = chain_a

    protein_atoms = [atom.index for atom in protein_chain.atoms]
    peptide_atoms = [atom.index for atom in peptide_chain.atoms]
    return np.array(protein_atoms), np.array(peptide_atoms)


def extract_protein_peptide(traj: md.Trajectory) -> md.Trajectory:
    """Extract only protein and peptide atoms (remove water/ions)."""
    protein_idx, peptide_idx = get_protein_peptide_indices(traj)
    all_idx = np.concatenate([protein_idx, peptide_idx])
    return traj.atom_slice(all_idx)


def save_trajectory_info(
    traj: md.Trajectory,
    output_path: Optional[Path] = None,
    time_step_ps: Optional[float] = None,
) -> None:
    """Save trajectory information to file."""
    if output_path is None:
        output_path = DEFAULT_DATA_DIR / "trajectory_info.txt"
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        f.write("CCT-OSR1 Trajectory Information\n")
        f.write("=" * 50 + "\n\n")
        timestep = _resolve_time_step(traj, override=time_step_ps)
        f.write(f"Total frames: {traj.n_frames}\n")
        if timestep is not None:
            total_time_ns = traj.n_frames * timestep / 1000
            f.write(f"Total time: {total_time_ns:.2f} ns\n")
            f.write(f"Time step: {timestep:.3f} ps/frame\n")
        else:
            f.write("Total time: N/A\n")
            f.write("Time step: N/A\n")
        f.write(f"Atoms: {traj.n_atoms}\n")
        f.write(f"Residues: {traj.n_residues}\n\n")
        
        f.write("Chain Information:\n")
        f.write("-" * 30 + "\n")
        for chain in traj.topology.chains:
            residues = list(chain.residues)
            if len(residues) > 0:
                f.write(f"Chain {chain.index}: {len(residues)} residues, ")
                f.write(f"{len([a for a in chain.atoms])} atoms\n")
    
    print(f"Trajectory info saved to: {output_path}")


# ============================================================================
# MAIN - Test loading
# ============================================================================
if __name__ == "__main__":
    # Test loading
    top = load_topology()
    
    print("\n")
    traj = load_and_concatenate_trajectories(stride=1)
    
    # Save info
    save_trajectory_info(traj)
    
    # Test selections
    protein_idx, peptide_idx = get_protein_peptide_indices(traj)
    print(f"\nProtein atoms: {len(protein_idx)}")
    print(f"Peptide atoms: {len(peptide_idx)}")

"""Utilities for runtime configuration of analysis scripts.

Provides a lightweight layer to resolve topology / trajectory paths and
output directories either from command line arguments, environment
variables or legacy configuration defaults. This keeps the scientific
logic of the analysis modules untouched while allowing them to be reused
in new systems without editing source code.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import argparse
import math
import os
import xml.etree.ElementTree as ET

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ANGSTROM_TO_NM = 0.1
DEFAULT_ENV_PREFIX = "CCT"

@dataclass
class AnalysisConfig:
    """Container for run-time configuration of an analysis task."""

    topology: Path
    trajectories: List[Path]
    output_root: Path
    label: str = "analysis"
    stride: int = 1
    time_step_ps: Optional[float] = None
    env_prefix: str = DEFAULT_ENV_PREFIX
    box_lengths_nm: Optional[Tuple[float, float, float]] = None
    # SLICE-F1 (D2): full 3x3 row-vector cell matrix in nm. Correct for
    # triclinic / truncated-octahedron cells where diagonal components are NOT
    # the edge lengths. When set, takes precedence over box_lengths_nm.
    box_vectors_nm: Optional[Tuple[Tuple[float, float, float], ...]] = None
    box_source: Optional[Path] = None
    autoimage_selection: Optional[str] = "chainid 0 and name CA"
    autoimage_enabled: bool = True

    figures_dir: Path = field(init=False)
    html_dir: Path = field(init=False)
    data_dir: Path = field(init=False)

    def __post_init__(self) -> None:
        self.topology = Path(self.topology)
        self.trajectories = [Path(t) for t in self.trajectories]
        self.output_root = Path(self.output_root)
        if self.box_source:
            self.box_source = Path(self.box_source)
        if self.box_lengths_nm is not None:
            self.box_lengths_nm = tuple(
                float(v) for v in self.box_lengths_nm
            )
        if self.box_vectors_nm is not None:
            vecs = tuple(
                tuple(float(x) for x in row) for row in self.box_vectors_nm
            )
            if len(vecs) != 3 or any(len(row) != 3 for row in vecs):
                raise ValueError("box_vectors_nm must be a 3x3 matrix")
            self.box_vectors_nm = vecs
            # Keep lengths consistent as the vector NORMS (D2): diagonal
            # components are NOT edge lengths for triclinic cells.
            norms = tuple(
                math.sqrt(sum(x * x for x in row)) for row in vecs
            )
            if self.box_lengths_nm is None:
                self.box_lengths_nm = norms

        self.figures_dir = self.output_root / "figures"
        self.html_dir = self.output_root / "html"
        self.data_dir = self.output_root / "data"

    def ensure_directories(self) -> None:
        """Create output folder structure if necessary."""

        self.output_root.mkdir(parents=True, exist_ok=True)
        self.figures_dir.mkdir(parents=True, exist_ok=True)
        self.html_dir.mkdir(parents=True, exist_ok=True)
        self.data_dir.mkdir(parents=True, exist_ok=True)

    @property
    def box_lengths_angstrom(self) -> Optional[Tuple[float, float, float]]:
        if self.box_lengths_nm is None:
            return None
        return tuple(val * 10.0 for val in self.box_lengths_nm)

    def should_autoimage(self) -> bool:
        return bool(
            self.autoimage_enabled
            and self.box_lengths_nm is not None
            and self.autoimage_selection
        )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _split_paths(raw: Optional[str]) -> List[Path]:
    if not raw:
        return []
    parts = []
    for token in raw.replace(",", os.pathsep).split(os.pathsep):
        token = token.strip()
        if token:
            parts.append(token)
    return [Path(p) for p in parts]


def _env(name: str, prefix: str) -> Optional[str]:
    return os.environ.get(f"{prefix}_{name}")


# ---------------------------------------------------------------------------
# CLI parser
# ---------------------------------------------------------------------------

def build_cli_parser(
    description: str,
    env_prefix: str = DEFAULT_ENV_PREFIX,
    add_label: bool = True,
) -> argparse.ArgumentParser:
    """Create a standard argparse parser for analysis scripts."""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--topology", type=Path, help="Topology file (PDB/PRMTOP/...)")
    parser.add_argument(
        "--trajectory",
        "--trajectories",
        nargs="+",
        dest="trajectories",
        type=Path,
        help="One or more trajectory files (DCD/XTC/...)",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        help="Root directory for analysis outputs",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=None,
        help="Frame stride when loading trajectories",
    )
    parser.add_argument(
        "--time-step",
        type=float,
        default=None,
        help="Time per saved frame in picoseconds (overrides metadata)",
    )

    if add_label:
        parser.add_argument(
            "--label",
            type=str,
            default=None,
            help="Label used in summary files/plots",
        )

    parser.add_argument(
        "--use-env",
        action="store_true",
        help=f"Pull defaults from environment variables with prefix {env_prefix}_*",
    )
    parser.add_argument(
        "--env-prefix",
        type=str,
        default=env_prefix,
        help="Custom environment variable prefix (default: CCT)",
    )
    parser.add_argument(
        "--box-xml",
        type=Path,
        help="OpenMM state XML with PeriodicBoxVectors (nm)",
    )
    parser.add_argument(
        "--box-length",
        type=float,
        default=None,
        help="Cubic box length in Angstrom (applies equally to X/Y/Z)",
    )
    parser.add_argument(
        "--box-lengths",
        type=float,
        nargs=3,
        metavar=("LX", "LY", "LZ"),
        help="Box lengths in Angstrom for X/Y/Z",
    )
    parser.add_argument(
        "--autoimage-selection",
        type=str,
        default=None,
        help="MDTraj selection string used as anchor for autoimage (default: protein C\u03b1)",
    )
    parser.add_argument(
        "--disable-autoimage",
        action="store_true",
        help="Disable anchor-based autoimage correction even if box is provided",
    )

    return parser


# ---------------------------------------------------------------------------
# Config resolution
# ---------------------------------------------------------------------------

def resolve_config(
    args: argparse.Namespace,
    *,
    default_topology: Optional[Path] = None,
    default_trajectories: Optional[Iterable[Path]] = None,
    default_output: Optional[Path] = None,
    default_stride: int = 1,
    default_time_step_ps: Optional[float] = None,
    default_label: Optional[str] = None,
    default_box_lengths_nm: Optional[Sequence[float]] = None,
    default_autoimage_selection: Optional[str] = "chainid 0 and name CA",
    default_autoimage_enabled: bool = True,
) -> AnalysisConfig:
    """Resolve an :class:`AnalysisConfig` from CLI arguments / environment."""

    prefix = args.env_prefix if getattr(args, "env_prefix", None) else DEFAULT_ENV_PREFIX

    env_top = _env("TOPOLOGY", prefix) if getattr(args, "use_env", False) else None
    env_traj = _env("TRAJECTORIES", prefix) if getattr(args, "use_env", False) else None
    env_out = _env("OUTPUT", prefix) if getattr(args, "use_env", False) else None
    env_stride = _env("STRIDE", prefix) if getattr(args, "use_env", False) else None
    env_time_step = _env("TIME_STEP_PS", prefix) if getattr(args, "use_env", False) else None
    env_label = _env("LABEL", prefix) if getattr(args, "use_env", False) else None
    env_box_xml = _env("BOX_XML", prefix) if getattr(args, "use_env", False) else None
    env_box_length = _env("BOX_LENGTH_A", prefix) if getattr(args, "use_env", False) else None
    env_box_lengths = _env("BOX_LENGTHS_A", prefix) if getattr(args, "use_env", False) else None
    env_autoimage = _env("AUTOIMAGE_SELECTION", prefix) if getattr(args, "use_env", False) else None
    env_autoimage_disable = _env("AUTOIMAGE_DISABLED", prefix) if getattr(args, "use_env", False) else None

    # Topology -----------------------------------------------------------
    topology = args.topology or env_top or default_topology
    if topology is None:
        raise ValueError("Topology file must be provided via --topology or environment")

    # Trajectories -------------------------------------------------------
    if args.trajectories:
        traj_values = [Path(t) for t in args.trajectories]
    elif env_traj:
        traj_values = _split_paths(env_traj)
    elif default_trajectories:
        traj_values = [Path(t) for t in default_trajectories]
    else:
        raise ValueError("At least one trajectory must be provided")

    # Output root --------------------------------------------------------
    output_root = (
        args.output_root
        or (Path(env_out) if env_out else None)
        or default_output
        or Path("outputs")
    )

    # Stride -------------------------------------------------------------
    stride = args.stride or (int(env_stride) if env_stride else None) or default_stride

    # Time step ----------------------------------------------------------
    if getattr(args, "time_step", None) is not None:
        time_step_ps = args.time_step
    elif env_time_step:
        time_step_ps = float(env_time_step)
    else:
        time_step_ps = default_time_step_ps

    # Label --------------------------------------------------------------
    label = args.label or env_label or default_label

    # Box vectors --------------------------------------------------------
    box_lengths_nm = None
    box_vectors_nm = None
    box_source = None

    def _set_box_from_xml(
        xml_path: Path,
    ) -> Tuple[
        Tuple[float, float, float],
        Tuple[Tuple[float, float, float], ...],
    ]:
        """Read PeriodicBoxVectors from an OpenMM state XML.

        Returns ``(lengths_nm, vectors_nm)``. ``lengths_nm`` are the vector
        NORMS (true edge lengths for orthorhombic cells; the correct minimum-
        image reference for triclinic/truncated-octahedron cells). ``vectors_nm``
        is the full 3x3 row-vector matrix in nm.

        SLICE-F1 (D2): the previous implementation read only the diagonal
        components (A.x, B.y, C.z), which is wrong for any non-orthogonal cell.
        """
        if not xml_path:
            return None
        if not xml_path.exists():
            raise FileNotFoundError(f"Box XML not found: {xml_path}")
        tree = ET.parse(str(xml_path))
        root = tree.getroot()
        vectors = root.find("PeriodicBoxVectors")
        if vectors is None:
            raise ValueError(f"PeriodicBoxVectors not found in {xml_path}")
        matrix: List[List[float]] = []
        for axis in ("A", "B", "C"):
            node = vectors.find(axis)
            if node is None:
                raise ValueError(
                    f"Axis {axis} missing from PeriodicBoxVectors in {xml_path}"
                )
            row = [
                float(node.attrib.get(comp, "0"))
                for comp in ("x", "y", "z")
            ]
            norm = math.sqrt(sum(v * v for v in row))
            if not math.isfinite(norm) or norm <= 0:
                raise ValueError(
                    f"Axis {axis} has non-positive/invalid norm in {xml_path}"
                )
            matrix.append(row)
        lengths = tuple(
            math.sqrt(sum(v * v for v in row)) for row in matrix
        )
        vectors_nm = tuple(tuple(row) for row in matrix)
        return lengths, vectors_nm

    def _guess_state_xml_from_traj(paths: Sequence[Path]) -> Optional[Path]:
        for traj_path in paths:
            folder = traj_path.parent
            if not folder.exists():
                continue
            direct = folder / f"{traj_path.stem}_state.xml"
            if direct.exists():
                return direct
            candidates = sorted(folder.glob("*state*.xml"))
            for candidate in candidates:
                if candidate.is_file():
                    return candidate
        return None

    if getattr(args, "box_xml", None):
        box_lengths_nm, box_vectors_nm = _set_box_from_xml(Path(args.box_xml))
        box_source = Path(args.box_xml)
    elif env_box_xml:
        box_lengths_nm, box_vectors_nm = _set_box_from_xml(Path(env_box_xml))
        box_source = Path(env_box_xml)

    if box_lengths_nm is None:
        guessed_xml = _guess_state_xml_from_traj(traj_values)
        if guessed_xml:
            try:
                box_lengths_nm, box_vectors_nm = _set_box_from_xml(guessed_xml)
                box_source = guessed_xml
                print(f"Detected box from {guessed_xml}")
            except Exception:
                pass

    def _parse_lengths(raw: Sequence[float]) -> Tuple[float, float, float]:
        vals = tuple(float(v) for v in raw)
        if len(vals) != 3:
            raise ValueError("Box lengths must have 3 values")
        return tuple(v * ANGSTROM_TO_NM for v in vals)

    if box_lengths_nm is None and getattr(args, "box_lengths", None):
        box_lengths_nm = _parse_lengths(args.box_lengths)
    elif box_lengths_nm is None and env_box_lengths:
        parts = [p.strip() for p in env_box_lengths.replace(",", " ").split() if p.strip()]
        if len(parts) != 3:
            raise ValueError("Environment BOX_LENGTHS_A must have three values")
        box_lengths_nm = _parse_lengths(parts)

    if box_lengths_nm is None and getattr(args, "box_length", None):
        length_nm = float(args.box_length) * ANGSTROM_TO_NM
        box_lengths_nm = (length_nm, length_nm, length_nm)
    elif box_lengths_nm is None and env_box_length:
        length_nm = float(env_box_length) * ANGSTROM_TO_NM
        box_lengths_nm = (length_nm, length_nm, length_nm)

    if box_lengths_nm is None and default_box_lengths_nm is not None:
        box_lengths_nm = tuple(float(v) for v in default_box_lengths_nm)

    # Autoimage ----------------------------------------------------------
    autoimage_selection = (
        args.autoimage_selection or env_autoimage or default_autoimage_selection
    )

    autoimage_enabled = not (
        getattr(args, "disable_autoimage", False)
        or (
            env_autoimage_disable
            and env_autoimage_disable.strip() not in {"", "false", "0", "False"}
        )
    )
    if not default_autoimage_enabled:
        autoimage_enabled = False

    cfg = AnalysisConfig(
        topology=Path(topology),
        trajectories=traj_values,
        output_root=Path(output_root),
        label=label,
        stride=stride,
        time_step_ps=time_step_ps,
        env_prefix=prefix,
        box_lengths_nm=box_lengths_nm,
        box_vectors_nm=box_vectors_nm,
        box_source=box_source,
        autoimage_selection=autoimage_selection,
        autoimage_enabled=autoimage_enabled,
    )

    cfg.ensure_directories()
    return cfg


def configure_from_cli(
    description,
    default_topology=None,
    default_trajectories=None,
    default_output=None,
    default_stride=1,
    default_time_step_ps=None,
    default_label=None,
) -> AnalysisConfig:
    parser = build_cli_parser(description)
    args = parser.parse_args()
    return resolve_config(
        args,
        default_topology=default_topology,
        default_trajectories=default_trajectories,
        default_output=default_output,
        default_stride=default_stride,
        default_time_step_ps=default_time_step_ps,
        default_label=default_label,
    )

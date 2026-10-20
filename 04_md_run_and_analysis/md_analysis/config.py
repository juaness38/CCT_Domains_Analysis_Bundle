"""
CCT-OSR1 Binding Analysis - Configuration
==========================================
WNK1-OSR1 Binding Project V1.0

System: CCT domain (OSR1) + WNK1 peptide (RFxV motif)
Trajectories: replica_1.dcd + replica_1_cont_cuda.dcd
"""

from pathlib import Path
import numpy as np

# ============================================================================
# PATHS
# ============================================================================
PROJECT_ROOT = Path(r"D:\DESCARGAS_PESADAS\CCT-OSR1_ANALYSIS")
DATA_DIR = PROJECT_ROOT / "data"
OUTPUT_DIR = PROJECT_ROOT / "outputs"
FIGURES_DIR = OUTPUT_DIR / "figures"
HTML_DIR = OUTPUT_DIR / "html"
DATA_OUTPUT_DIR = OUTPUT_DIR / "data"

# Topology
TOPOLOGY_PDB = DATA_DIR / "topology.pdb"

# Trajectories (in order)
TRAJECTORIES = [
    Path(r"D:\DESCARGAS_PESADAS\replica_1.dcd"),
    Path(r"D:\DESCARGAS_PESADAS\replica_1_cont_cuda.dcd"),
]

# ============================================================================
# SYSTEM DEFINITION
# ============================================================================
# Chain A: CCT domain of OSR1 (96 residues, 1470 atoms)
# Chain B: WNK1 peptide with RFxV motif (17 residues, 274 atoms)
# Chain C: Water (HOH)
# Chain D: Ions (Na+, Cl-)

CHAIN_DEFINITIONS = {
    'protein': {'chain': 'A', 'resid_range': (1, 96), 'name': 'CCT-OSR1'},
    'peptide': {'chain': 'B', 'resid_range': (1, 17), 'name': 'WNK1-peptide'},
}

# RFxV motif in peptide (based on sequence LTQVVHSAGRRFIVSPV)
RFXV_MOTIF = {
    'residues': [10, 11, 12, 13, 14],  # R-R-F-I-V
    'core': [11, 12, 13, 14],  # R-F-I-V (canonical)
    'sequence': 'RRFIV',
}

# Key residues for analysis (from hypothesis)
# Phase 1: Electrostatic steering to Glu/Ser 37-40
ELECTROSTATIC_ZONE = list(range(37, 41))  # residues 37-40 in CCT

# Hydrophobic pocket residues (for Phe insertion)
HYDROPHOBIC_POCKET = [15, 17, 19, 45, 47, 49, 51]  # approximate, refine with structure

# ============================================================================
# SIMULATION PARAMETERS
# ============================================================================
TIMESTEP_FS = 4.0  # femtoseconds (HMR enabled)
SAVE_FREQ_STEPS = 12500  # frames saved every 12500 steps
SAVE_FREQ_PS = SAVE_FREQ_STEPS * TIMESTEP_FS / 1000  # = 50 ps per frame

TEMPERATURE_K = 300.0
BOX_SIZE_NM = 10.3  # approximate cubic box

# ============================================================================
# ANALYSIS PARAMETERS
# ============================================================================
# RMSD/RMSF
RMSD_SELECTION = {
    'protein_backbone': 'protein and backbone',
    'protein_ca': 'protein and name CA',
    'peptide_backbone': 'chainid 1 and backbone',
    'peptide_all': 'chainid 1',
    'rfxv_motif': 'chainid 1 and resid 10 to 14',
}

# Contact analysis
CONTACT_CUTOFF_NM = 0.45  # 4.5 Angstroms
CONTACT_HEAVY_ONLY = True

# Water shell analysis
WATER_SHELL_CUTOFFS = [0.35, 0.5, 0.8]  # nm - first, second, extended shell

# Clustering
CLUSTER_RMSD_CUTOFF = 0.2  # nm
MIN_CLUSTER_SIZE = 10  # frames

# Binding thresholds
BINDING_DISTANCE_NM = 2.0
UNBOUND_DISTANCE_NM = 3.0
MIN_BOUND_DURATION_NS = 1.0

# ============================================================================
# VISUALIZATION
# ============================================================================
PLOTLY_TEMPLATE = 'plotly_white'
COLOR_SCHEME = {
    'protein': '#1f77b4',  # blue
    'peptide': '#ff7f0e',  # orange
    'rfxv': '#d62728',     # red
    'water': '#17becf',    # cyan
    'contact': '#2ca02c',  # green
    'highlight': '#9467bd', # purple
}

# Figure sizes
FIG_WIDTH = 1200
FIG_HEIGHT = 600
FIG_HEIGHT_TALL = 800

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================
def time_to_ns(frame_idx: int) -> float:
    """Convert frame index to time in nanoseconds."""
    return frame_idx * SAVE_FREQ_PS / 1000

def frames_to_ns(n_frames: int) -> float:
    """Convert number of frames to nanoseconds."""
    return n_frames * SAVE_FREQ_PS / 1000

def ns_to_frames(ns: float) -> int:
    """Convert nanoseconds to frame count."""
    return int(ns * 1000 / SAVE_FREQ_PS)

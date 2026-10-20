"""
SMIC Visualization Module
==========================

Publication-quality visualization for IFP analysis results. 

Author: MICA Team
License: MIT
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, List, Tuple, Union, TYPE_CHECKING
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)

# Type checking imports
if TYPE_CHECKING:
    from analysis_ifp import IFPTrajectoryResult

# Optional imports with availability flags
try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.gridspec import GridSpec
    import matplotlib.ticker as ticker
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    logger.warning("matplotlib not available")

try:
    import seaborn as sns
    SEABORN_AVAILABLE = True
except ImportError:
    SEABORN_AVAILABLE = False

try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    logger.warning("plotly not available - HTML export disabled")


# =============================================================================
# STYLE CONFIGURATION
# =============================================================================

@dataclass
class PlotStyle:
    """Publication-ready plot style configuration."""
    # Figure defaults
    figsize_single: Tuple[float, float] = (8, 6)
    figsize_wide: Tuple[float, float] = (12, 6)
    figsize_tall: Tuple[float, float] = (8, 10)
    figsize_dashboard: Tuple[float, float] = (16, 12)
    
    # Resolution
    dpi_screen: int = 100
    dpi_print: int = 300
    
    # Fonts (Nature/Science style)
    font_family: str = 'sans-serif'
    font_size_title: int = 12
    font_size_label: int = 10
    font_size_tick: int = 9
    font_size_legend: int = 9
    font_size_annotation: int = 8
    
    # Colors
    background: str = 'white'
    grid_color: str = '#E0E0E0'
    grid_alpha: float = 0.5
    
    # Line widths
    line_width: float = 1.5
    spine_width: float = 0.8
    
    def apply(self):
        """Apply style to matplotlib."""
        if not MATPLOTLIB_AVAILABLE:
            return
        
        plt.rcParams.update({
            'font.family': self.font_family,
            'font.size': self.font_size_tick,
            'axes.titlesize': self.font_size_title,
            'axes.labelsize': self.font_size_label,
            'xtick.labelsize': self.font_size_tick,
            'ytick.labelsize': self.font_size_tick,
            'legend.fontsize': self.font_size_legend,
            'figure.facecolor': self.background,
            'axes.facecolor': self.background,
            'axes.linewidth': self.spine_width,
            'axes.grid': True,
            'grid.alpha': self.grid_alpha,
            'grid.color': self.grid_color,
            'savefig.dpi': self.dpi_print,
            'savefig.bbox': 'tight',
            'savefig.pad_inches': 0.1,
        })


# Default style instance
STYLE = PlotStyle()

# IFP color scheme (consistent across all plots)
IFP_COLORS = {
    'AR': '#1abc9c',  # Teal - Aromatic
    'HY': '#27ae60',  # Green - Hydrophobic
    'HD': '#e74c3c',  # Red - H-bond Donor (protein donates)
    'HA': '#3498db',  # Blue - H-bond Acceptor (protein accepts)
    'WB': '#00bcd4',  # Cyan - Water Bridge
    'IP': '#f39c12',  # Orange - Salt Bridge (+)
    'IN': '#9b59b6',  # Purple - Salt Bridge (-)
    'IO': '#e91e63',  # Pink - Metal Ion
    'HL': '#ff5722',  # Deep Orange - Halogen
}

IFP_LABELS = {
    'AR': 'Aromatic (π)',
    'HY': 'Hydrophobic',
    'HD': 'H-bond (P→L)',  # Protein donor
    'HA': 'H-bond (L→P)',  # Ligand/peptide donor
    'WB': 'Water bridge',
    'IP': 'Salt bridge (+)',
    'IN': 'Salt bridge (−)',
    'IO': 'Metal coord.',
    'HL': 'Halogen bond',
}

# Order for consistent plotting
IFP_ORDER = ['HD', 'HA', 'WB', 'HY', 'AR', 'IP', 'IN', 'HL', 'IO']


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def _check_matplotlib():
    """Check matplotlib availability."""
    if not MATPLOTLIB_AVAILABLE:
        raise ImportError(
            "matplotlib required for visualization. "
            "Install: pip install matplotlib"
        )


def _get_ifp_type(column: str) -> str:
    """Extract IFP type from column name."""
    return column.split('_')[0] if '_' in column else column[:2]


def _get_residue(column: str) -> str:
    """Extract residue from column name."""
    parts = column.split('_')
    return '_'.join(parts[1:]) if len(parts) > 1 else column


def _sort_ifp_columns(columns: List[str], occupancies: Dict[str, float] = None) -> List[str]:
    """Sort IFP columns by type order, then by residue number or occupancy."""
    def sort_key(col):
        ifp_type = _get_ifp_type(col)
        type_order = IFP_ORDER.index(ifp_type) if ifp_type in IFP_ORDER else 99
        
        # Try to extract residue number
        residue = _get_residue(col)
        try:
            resnum = int(''.join(filter(str.isdigit, residue)))
        except ValueError:
            resnum = 999
        
        # Use negative occupancy for descending sort
        occ = -(occupancies.get(col, 0) if occupancies else 0)
        
        return (type_order, resnum, occ)
    
    return sorted(columns, key=sort_key)


def _save_figure(
    fig: "plt.Figure",
    save_path: Union[str, Path],
    formats: List[str] = None,
    dpi: int = None,
) -> List[Path]:
    """
    Save figure in multiple formats.
    
    Parameters
    ----------
    fig : matplotlib.Figure
    save_path : Path
        Base path (extension will be replaced)
    formats : list
        Output formats (default: ['png', 'pdf'])
    dpi : int
        Resolution for raster formats
    
    Returns
    -------
    list of Path
        Saved file paths
    """
    if formats is None:
        formats = ['png', 'pdf']
    
    if dpi is None:
        dpi = STYLE.dpi_print
    
    save_path = Path(save_path)
    base = save_path.parent / save_path.stem
    
    saved = []
    for fmt in formats:
        out_path = base.with_suffix(f'.{fmt}')
        fig.savefig(out_path, format=fmt, dpi=dpi, bbox_inches='tight')
        saved.append(out_path)
        logger.debug(f"Saved {out_path}")
    
    return saved


# =============================================================================
# INDIVIDUAL PLOT FUNCTIONS
# =============================================================================

def plot_ifp_heatmap(
    ifp_result: "IFPTrajectoryResult",
    ifp_types: Optional[List[str]] = None,
    min_occupancy: float = 0.05,
    max_residues: int = 40,
    time_unit: str = 'ns',
    cmap: str = 'YlGnBu',
    figsize: Tuple[float, float] = None,
    title: Optional[str] = None,
    ax: Optional["plt.Axes"] = None,
    save_path: Optional[Union[str, Path]] = None,
    show: bool = True,
) -> Optional["plt.Figure"]:
    """
    Plot IFP heatmap (IFP columns × time frames).
    
    Parameters
    ----------
    ifp_result : IFPTrajectoryResult
        Analysis results from analysis_ifp
    ifp_types : list, optional
        Filter to specific IFP types (e.g., ['HD', 'HA'])
    min_occupancy : float
        Minimum occupancy to include (0-1)
    max_residues : int
        Maximum number of IFP rows to show
    time_unit : str
        'ns' or 'ps' for x-axis
    cmap : str
        Matplotlib colormap
    figsize : tuple, optional
        Figure size (auto-calculated if None)
    title : str, optional
        Plot title
    ax : matplotlib.Axes, optional
        Existing axes to plot on
    save_path : Path, optional
        Save path (saves both PNG and PDF)
    show : bool
        Display figure
    
    Returns
    -------
    fig : matplotlib.Figure or None
    """
    _check_matplotlib()
    STYLE.apply()
    
    # Get data from result
    df = ifp_result.ifp_dataframe. copy()
    ifp_cols = [c for c in ifp_result.ifp_columns if c in df. columns]
    
    # Filter by IFP type
    if ifp_types:
        ifp_cols = [c for c in ifp_cols 
                    if _get_ifp_type(c) in ifp_types]
    
    # Filter by occupancy
    occupancies = ifp_result.ifp_occupancies
    ifp_cols = [c for c in ifp_cols 
                if occupancies.get(c, 0) >= min_occupancy]
    
    if not ifp_cols:
        logger.warning(f"No IFPs pass filters (min_occ={min_occupancy})")
        return None
    
    # Sort and limit
    ifp_cols = _sort_ifp_columns(ifp_cols, occupancies)[:max_residues]
    
    # Build matrix
    matrix = df[ifp_cols].values. T  # IFPs × frames
    
    # Time axis
    time = ifp_result.time_ps
    if time_unit == 'ns':
        time = time / 1000
        time_label = 'Time (ns)'
    else:
        time_label = 'Time (ps)'
    
    # Create figure if needed
    created_fig = ax is None
    if created_fig:
        n_rows = len(ifp_cols)
        height = max(4, min(12, n_rows * 0.25 + 2))
        figsize = figsize or (10, height)
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    
    # Plot heatmap
    im = ax.imshow(
        matrix,
        aspect='auto',
        cmap=cmap,
        interpolation='nearest',
        extent=[time[0], time[-1], len(ifp_cols) - 0.5, -0.5],
        vmin=0, vmax=1,
    )
    
    # Y-axis: IFP labels with colored text
    ax.set_yticks(range(len(ifp_cols)))
    ax.set_yticklabels(ifp_cols, fontsize=STYLE.font_size_tick - 1)
    
    for i, col in enumerate(ifp_cols):
        color = IFP_COLORS.get(_get_ifp_type(col), '#666666')
        ax.get_yticklabels()[i].set_color(color)
    
    ax.set_xlabel(time_label)
    ax.set_ylabel('Interaction Fingerprint')
    
    if title:
        ax.set_title(title)
    
    # Colorbar
    if created_fig:
        cbar = plt.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
        cbar.set_label('Present', fontsize=STYLE.font_size_label)
        cbar.set_ticks([0, 1])
    
    if created_fig:
        plt.tight_layout()
    
    # Save
    if save_path and created_fig:
        _save_figure(fig, save_path)
        logger.info(f"Saved heatmap to {save_path}")
    
    if show and created_fig:
        plt.show()
    
    return fig if created_fig else None


def plot_ifp_timeseries(
    ifp_result: "IFPTrajectoryResult",
    ifp_types: Optional[List[str]] = None,
    smoothing: int = 10,
    stacked: bool = False,
    time_unit: str = 'ns',
    figsize: Tuple[float, float] = None,
    title: Optional[str] = None,
    ax: Optional["plt.Axes"] = None,
    save_path: Optional[Union[str, Path]] = None,
    show: bool = True,
) -> Optional["plt.Figure"]:
    """
    Plot IFP count time series by type.
    
    Parameters
    ----------
    ifp_result : IFPTrajectoryResult
    ifp_types : list, optional
        IFP types to plot (default: all detected)
    smoothing : int
        Rolling average window (frames)
    stacked : bool
        Use stacked area plot instead of lines
    time_unit : str
        'ns' or 'ps'
    figsize : tuple
    title : str
    ax : matplotlib. Axes
    save_path : Path
    show : bool
    
    Returns
    -------
    fig : matplotlib.Figure or None
    """
    _check_matplotlib()
    STYLE.apply()
    
    # Determine IFP types to plot
    if ifp_types is None:
        detected = set(_get_ifp_type(c) for c in ifp_result.ifp_columns)
        ifp_types = [t for t in IFP_ORDER if t in detected]
    
    # Time axis
    time = ifp_result.time_ps
    if time_unit == 'ns':
        time = time / 1000
        time_label = 'Time (ns)'
    else:
        time_label = 'Time (ps)'
    
    # Create figure
    created_fig = ax is None
    if created_fig:
        figsize = figsize or STYLE.figsize_wide
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    
    # Collect time series
    series_data = {}
    for ifp_type in ifp_types:
        counts = ifp_result.get_ifp_timeseries(ifp_type)
        if counts.max() == 0:
            continue
        
        # Smooth
        if smoothing > 1 and len(counts) >= smoothing:
            counts = pd.Series(counts). rolling(
                window=smoothing, center=True, min_periods=1
            ).mean(). values
        
        series_data[ifp_type] = counts
    
    if not series_data:
        logger.warning("No IFP data to plot")
        return None
    
    # Plot
    if stacked:
        labels = list(series_data.keys())
        colors = [IFP_COLORS.get(t, '#808080') for t in labels]
        ax.stackplot(
            time, 
            *series_data.values(),
            labels=[IFP_LABELS.get(t, t) for t in labels],
            colors=colors,
            alpha=0.8,
        )
    else:
        for ifp_type, counts in series_data.items():
            ax.plot(
                time, counts,
                label=IFP_LABELS.get(ifp_type, ifp_type),
                color=IFP_COLORS.get(ifp_type, '#808080'),
                linewidth=STYLE.line_width,
            )
    
    ax.set_xlabel(time_label)
    ax.set_ylabel('Number of interactions')
    ax.legend(loc='upper right', framealpha=0.9)
    ax.set_xlim(time[0], time[-1])
    ax.set_ylim(bottom=0)
    
    if title:
        ax.set_title(title)
    
    if created_fig:
        plt.tight_layout()
    
    if save_path and created_fig:
        _save_figure(fig, save_path)
    
    if show and created_fig:
        plt.show()
    
    return fig if created_fig else None


def plot_ifp_occupancy(
    ifp_result: "IFPTrajectoryResult",
    min_occupancy: float = 0.05,
    max_bars: int = 25,
    orientation: str = 'horizontal',
    show_values: bool = True,
    figsize: Tuple[float, float] = None,
    title: Optional[str] = None,
    ax: Optional["plt.Axes"] = None,
    save_path: Optional[Union[str, Path]] = None,
    show: bool = True,
) -> Optional["plt.Figure"]:
    """
    Plot IFP occupancy bar chart.
    
    Parameters
    ----------
    ifp_result : IFPTrajectoryResult
    min_occupancy : float
        Minimum occupancy to show (0-1)
    max_bars : int
        Maximum bars to display
    orientation : str
        'horizontal' or 'vertical'
    show_values : bool
        Show percentage labels on bars
    figsize : tuple
    title : str
    ax : matplotlib.Axes
    save_path : Path
    show : bool
    
    Returns
    -------
    fig : matplotlib.Figure or None
    """
    _check_matplotlib()
    STYLE.apply()
    
    # Filter and sort
    occs = ifp_result.ifp_occupancies
    filtered = [(k, v) for k, v in occs.items() if v >= min_occupancy]
    filtered.sort(key=lambda x: -x[1])
    filtered = filtered[:max_bars]
    
    if not filtered:
        logger.warning(f"No IFPs with occupancy >= {min_occupancy*100:. 0f}%")
        return None
    
    labels = [x[0] for x in filtered]
    values = [x[1] * 100 for x in filtered]  # Percentage
    colors = [IFP_COLORS.get(_get_ifp_type(l), '#808080') for l in labels]
    
    # Create figure
    created_fig = ax is None
    if created_fig:
        if orientation == 'horizontal':
            height = max(4, len(labels) * 0.3 + 1)
            figsize = figsize or (8, height)
        else:
            figsize = figsize or (max(6, len(labels) * 0.4), 6)
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    
    # Plot
    if orientation == 'horizontal':
        y_pos = range(len(labels))
        bars = ax.barh(y_pos, values, color=colors, edgecolor='white', linewidth=0.5)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels)
        ax.invert_yaxis()
        ax.set_xlabel('Occupancy (%)')
        ax.set_xlim(0, 105)
        
        if show_values:
            for bar, val in zip(bars, values):
                ax.text(
                    val + 1, bar.get_y() + bar.get_height()/2,
                    f'{val:.0f}%', va='center', fontsize=STYLE.font_size_annotation
                )
    else:
        x_pos = range(len(labels))
        bars = ax.bar(x_pos, values, color=colors, edgecolor='white', linewidth=0.5)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(labels, rotation=45, ha='right')
        ax.set_ylabel('Occupancy (%)')
        ax.set_ylim(0, 105)
        
        if show_values:
            for bar, val in zip(bars, values):
                ax.text(
                    bar.get_x() + bar.get_width()/2, val + 1,
                    f'{val:.0f}%', ha='center', fontsize=STYLE.font_size_annotation
                )
    
    if title:
        ax.set_title(title)
    
    if created_fig:
        plt.tight_layout()
    
    if save_path and created_fig:
        _save_figure(fig, save_path)
    
    if show and created_fig:
        plt.show()
    
    return fig if created_fig else None


def plot_residue_profile(
    ifp_result: "IFPTrajectoryResult",
    min_occupancy: float = 0.1,
    max_residues: int = 20,
    figsize: Tuple[float, float] = None,
    title: Optional[str] = None,
    ax: Optional["plt.Axes"] = None,
    save_path: Optional[Union[str, Path]] = None,
    show: bool = True,
) -> Optional["plt.Figure"]:
    """
    Plot per-residue IFP profile (stacked bars).
    
    Shows all IFP types for each residue. 
    """
    _check_matplotlib()
    STYLE.apply()
    
    # Group by residue
    residue_data = {}
    for col, occ in ifp_result.ifp_occupancies.items():
        ifp_type = _get_ifp_type(col)
        residue = _get_residue(col)
        
        if residue not in residue_data:
            residue_data[residue] = {}
        residue_data[residue][ifp_type] = occ * 100
    
    # Filter and sort
    residue_totals = {r: sum(d.values()) for r, d in residue_data.items()}
    filtered = [(r, residue_totals[r]) for r in residue_data 
                if residue_totals[r] >= min_occupancy * 100]
    filtered.sort(key=lambda x: -x[1])
    residues = [x[0] for x in filtered[:max_residues]]
    
    if not residues:
        logger.warning("No residues pass filter")
        return None
    
    # Create figure
    created_fig = ax is None
    if created_fig:
        figsize = figsize or (max(8, len(residues) * 0.5), 6)
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    
    # Stacked bars
    x = np.arange(len(residues))
    width = 0.7
    bottom = np.zeros(len(residues))
    
    # Determine which types are present
    all_types = set()
    for r in residues:
        all_types.update(residue_data[r].keys())
    ordered_types = [t for t in IFP_ORDER if t in all_types]
    
    for ifp_type in ordered_types:
        values = [residue_data[r].get(ifp_type, 0) for r in residues]
        ax.bar(
            x, values, width, bottom=bottom,
            color=IFP_COLORS.get(ifp_type, '#808080'),
            label=IFP_LABELS.get(ifp_type, ifp_type),
            edgecolor='white', linewidth=0.5,
        )
        bottom += values
    
    ax.set_xticks(x)
    ax.set_xticklabels(residues, rotation=45, ha='right')
    ax.set_ylabel('Occupancy (%)')
    ax.set_xlabel('Residue')
    ax.legend(loc='upper right', fontsize=STYLE.font_size_legend)
    
    if title:
        ax.set_title(title)
    
    if created_fig:
        plt.tight_layout()
    
    if save_path and created_fig:
        _save_figure(fig, save_path)
    
    if show and created_fig:
        plt.show()
    
    return fig if created_fig else None


# =============================================================================
# DASHBOARD / MULTIPANEL FIGURES
# =============================================================================

def create_ifp_figure(
    ifp_result: "IFPTrajectoryResult",
    layout: str = 'default',
    figsize: Tuple[float, float] = None,
    title: Optional[str] = None,
    save_path: Optional[Union[str, Path]] = None,
    show: bool = True,
) -> Optional["plt.Figure"]:
    """
    Create publication-ready multipanel IFP figure. 
    
    Parameters
    ----------
    ifp_result : IFPTrajectoryResult
    layout : str
        'default' : 2×2 grid (heatmap, timeseries, occupancy, residue)
        'wide' : 1×3 (timeseries, occupancy, residue)
        'minimal' : 1×2 (heatmap, occupancy)
    figsize : tuple
    title : str
        Supertitle
    save_path : Path
    show : bool
    
    Returns
    -------
    fig : matplotlib.Figure
    """
    _check_matplotlib()
    STYLE.apply()
    
    if layout == 'default':
        figsize = figsize or (14, 10)
        fig = plt.figure(figsize=figsize)
        gs = GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.25)
        
        # A: Heatmap (top-left, spans both rows height-wise ideally, but here top-left)
        ax_heatmap = fig.add_subplot(gs[0, 0])
        plot_ifp_heatmap(ifp_result, ax=ax_heatmap, show=False, max_residues=25)
        ax_heatmap.set_title('A) IFP Heatmap', loc='left', fontweight='bold')
        
        # B: Time series (top-right)
        ax_ts = fig.add_subplot(gs[0, 1])
        plot_ifp_timeseries(ifp_result, ax=ax_ts, show=False)
        ax_ts.set_title('B) Time Series', loc='left', fontweight='bold')
        
        # C: Occupancy (bottom-left)
        ax_occ = fig.add_subplot(gs[1, 0])
        plot_ifp_occupancy(ifp_result, ax=ax_occ, show=False, max_bars=15)
        ax_occ.set_title('C) Occupancy', loc='left', fontweight='bold')
        
        # D: Residue profile (bottom-right)
        ax_res = fig.add_subplot(gs[1, 1])
        plot_residue_profile(ifp_result, ax=ax_res, show=False, max_residues=15)
        ax_res.set_title('D) Residue Profile', loc='left', fontweight='bold')
        
    elif layout == 'wide':
        figsize = figsize or (16, 5)
        fig, axes = plt.subplots(1, 3, figsize=figsize)
        
        plot_ifp_timeseries(ifp_result, ax=axes[0], show=False)
        axes[0].set_title('A) Time Series', loc='left', fontweight='bold')
        
        plot_ifp_occupancy(ifp_result, ax=axes[1], show=False, 
                          orientation='vertical', max_bars=12)
        axes[1].set_title('B) Occupancy', loc='left', fontweight='bold')
        
        plot_residue_profile(ifp_result, ax=axes[2], show=False, max_residues=12)
        axes[2].set_title('C) Residue Profile', loc='left', fontweight='bold')
        
    elif layout == 'minimal':
        figsize = figsize or (12, 5)
        fig, axes = plt.subplots(1, 2, figsize=figsize)
        
        plot_ifp_heatmap(ifp_result, ax=axes[0], show=False, max_residues=20)
        axes[0].set_title('A) IFP Heatmap', loc='left', fontweight='bold')
        
        plot_ifp_occupancy(ifp_result, ax=axes[1], show=False, max_bars=15)
        axes[1].set_title('B) Occupancy', loc='left', fontweight='bold')
    
    else:
        raise ValueError(f"Unknown layout: {layout}")
    
    if title:
        fig.suptitle(title, fontsize=STYLE.font_size_title + 2, fontweight='bold', y=1.02)
    
    plt.tight_layout()
    
    if save_path:
        _save_figure(fig, save_path)
        logger.info(f"Saved figure to {save_path}")
    
    if show:
        plt.show()
    
    return fig


def generate_all_plots(
    ifp_result: "IFPTrajectoryResult",
    output_dir: Union[str, Path],
    prefix: str = "ifp",
    formats: List[str] = None,
    min_occupancy: float = 0.05,
) -> Dict[str, List[Path]]:
    """
    Generate all individual plots and dashboard. 
    
    Parameters
    ----------
    ifp_result : IFPTrajectoryResult
    output_dir : Path
    prefix : str
    formats : list
        Output formats (default: ['png', 'pdf'])
    min_occupancy : float
    
    Returns
    -------
    dict
        {plot_type: [saved_paths]}
    """
    _check_matplotlib()
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    if formats is None:
        formats = ['png', 'pdf']
    
    results = {}
    
    # Individual plots
    plots = [
        ('heatmap', plot_ifp_heatmap, {'min_occupancy': min_occupancy}),
        ('timeseries', plot_ifp_timeseries, {}),
        ('timeseries_stacked', plot_ifp_timeseries, {'stacked': True}),
        ('occupancy', plot_ifp_occupancy, {'min_occupancy': min_occupancy}),
        ('residue_profile', plot_residue_profile, {'min_occupancy': min_occupancy}),
    ]
    
    for name, func, kwargs in plots:
        try:
            save_path = output_dir / f"{prefix}_{name}.png"
            fig = func(ifp_result, save_path=save_path, show=False, **kwargs)
            if fig:
                results[name] = _save_figure(fig, save_path, formats)
                plt.close(fig)
        except Exception as e:
            logger.warning(f"Failed to create {name}: {e}")
    
    # Dashboard layouts
    for layout in ['default', 'wide', 'minimal']:
        try:
            save_path = output_dir / f"{prefix}_figure_{layout}.png"
            fig = create_ifp_figure(ifp_result, layout=layout, save_path=save_path, show=False)
            if fig:
                results[f'figure_{layout}'] = _save_figure(fig, save_path, formats)
                plt.close(fig)
        except Exception as e:
            logger.warning(f"Failed to create {layout} figure: {e}")
    
    logger.info(f"Generated {len(results)} plot sets in {output_dir}")
    return results


# =============================================================================
# PLOTLY HTML INTERACTIVE VISUALIZATIONS
# =============================================================================

def _check_plotly():
    """Check Plotly availability."""
    if not PLOTLY_AVAILABLE:
        raise ImportError(
            "Plotly required for HTML export. "
            "Install: pip install plotly"
        )


def plot_ifp_heatmap_html(
    ifp_result: "IFPTrajectoryResult",
    ifp_types: Optional[List[str]] = None,
    min_occupancy: float = 0.05,
    max_residues: int = 50,
    time_unit: str = 'ns',
    title: Optional[str] = None,
    save_path: Optional[Union[str, Path]] = None,
) -> "go.Figure":
    """
    Create interactive IFP heatmap with Plotly (HTML export).
    
    Parameters
    ----------
    ifp_result : IFPTrajectoryResult
        Analysis results from analysis_ifp
    ifp_types : list, optional
        Filter to specific IFP types (e.g., ['HD', 'HA'])
    min_occupancy : float
        Minimum occupancy to include (0-1)
    max_residues : int
        Maximum number of IFP rows
    time_unit : str
        'ns' or 'ps' for x-axis
    title : str, optional
        Plot title
    save_path : Path, optional
        Path to save HTML file
        
    Returns
    -------
    go.Figure
        Plotly figure object
    """
    _check_plotly()
    
    # Get data
    df = ifp_result.ifp_dataframe.copy()
    ifp_cols = [c for c in ifp_result.ifp_columns if c in df.columns]
    
    # Filter by IFP type
    if ifp_types:
        ifp_cols = [c for c in ifp_cols if _get_ifp_type(c) in ifp_types]
    
    # Filter by occupancy
    occupancies = ifp_result.ifp_occupancies
    ifp_cols = [c for c in ifp_cols if occupancies.get(c, 0) >= min_occupancy]
    
    if not ifp_cols:
        logger.warning(f"No IFPs pass filters (min_occ={min_occupancy})")
        return None
    
    # Sort and limit
    ifp_cols = _sort_ifp_columns(ifp_cols, occupancies)[:max_residues]
    
    # Build matrix
    matrix = df[ifp_cols].values.T  # IFPs × frames
    
    # Time axis
    time = ifp_result.time_ps
    if time_unit == 'ns':
        time = time / 1000
        time_label = 'Time (ns)'
    else:
        time_label = 'Time (ps)'
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=matrix,
        x=time,
        y=ifp_cols,
        colorscale='YlGnBu',
        zmin=0,
        zmax=1,
        hovertemplate='<b>%{y}</b><br>Time: %{x:.2f}<br>Value: %{z}<extra></extra>',
        colorbar=dict(title='Present', tickvals=[0, 1])
    ))
    
    fig.update_layout(
        title=title or 'IFP Heatmap (Interactive)',
        xaxis_title=time_label,
        yaxis_title='Interaction Fingerprint',
        height=max(400, len(ifp_cols) * 20 + 100),
        template='plotly_white',
        hoverlabel=dict(bgcolor="white"),
    )
    
    # Save if path provided
    if save_path:
        save_path = Path(save_path)
        if not save_path.suffix:
            save_path = save_path.with_suffix('.html')
        fig.write_html(str(save_path), include_plotlyjs=True)
        logger.info(f"Saved interactive heatmap to {save_path}")
    
    return fig


def plot_ifp_timeseries_html(
    ifp_result: "IFPTrajectoryResult",
    ifp_types: Optional[List[str]] = None,
    smoothing: int = 10,
    time_unit: str = 'ns',
    title: Optional[str] = None,
    save_path: Optional[Union[str, Path]] = None,
) -> "go.Figure":
    """
    Create interactive IFP time series with Plotly.
    
    Parameters
    ----------
    ifp_result : IFPTrajectoryResult
    ifp_types : list, optional
        IFP types to plot (default: all detected)
    smoothing : int
        Rolling average window (frames)
    time_unit : str
        'ns' or 'ps'
    title : str
    save_path : Path
        
    Returns
    -------
    go.Figure
    """
    _check_plotly()
    
    # Determine IFP types to plot
    if ifp_types is None:
        detected = set(_get_ifp_type(c) for c in ifp_result.ifp_columns)
        ifp_types = [t for t in IFP_ORDER if t in detected]
    
    # Time axis
    time = ifp_result.time_ps
    if time_unit == 'ns':
        time = time / 1000
        time_label = 'Time (ns)'
    else:
        time_label = 'Time (ps)'
    
    fig = go.Figure()
    
    for ifp_type in ifp_types:
        counts = ifp_result.get_ifp_timeseries(ifp_type)
        if counts.max() == 0:
            continue
        
        # Smooth
        if smoothing > 1 and len(counts) >= smoothing:
            counts = pd.Series(counts).rolling(
                window=smoothing, center=True, min_periods=1
            ).mean().values
        
        fig.add_trace(go.Scatter(
            x=time,
            y=counts,
            mode='lines',
            name=IFP_LABELS.get(ifp_type, ifp_type),
            line=dict(color=IFP_COLORS.get(ifp_type, '#808080'), width=2),
            hovertemplate=f'{IFP_LABELS.get(ifp_type, ifp_type)}<br>' +
                          'Time: %{x:.2f}<br>Count: %{y:.1f}<extra></extra>'
        ))
    
    fig.update_layout(
        title=title or 'IFP Time Series (Interactive)',
        xaxis_title=time_label,
        yaxis_title='Number of interactions',
        template='plotly_white',
        hovermode='x unified',
        legend=dict(yanchor="top", y=0.99, xanchor="right", x=0.99),
    )
    
    if save_path:
        save_path = Path(save_path)
        if not save_path.suffix:
            save_path = save_path.with_suffix('.html')
        fig.write_html(str(save_path), include_plotlyjs=True)
        logger.info(f"Saved interactive timeseries to {save_path}")
    
    return fig


def plot_ifp_occupancy_html(
    ifp_result: "IFPTrajectoryResult",
    min_occupancy: float = 0.05,
    max_bars: int = 30,
    title: Optional[str] = None,
    save_path: Optional[Union[str, Path]] = None,
) -> "go.Figure":
    """
    Create interactive IFP occupancy bar chart with Plotly.
    
    Parameters
    ----------
    ifp_result : IFPTrajectoryResult
    min_occupancy : float
    max_bars : int
    title : str
    save_path : Path
        
    Returns
    -------
    go.Figure
    """
    _check_plotly()
    
    # Filter and sort
    occs = ifp_result.ifp_occupancies
    filtered = [(k, v) for k, v in occs.items() if v >= min_occupancy]
    filtered.sort(key=lambda x: -x[1])
    filtered = filtered[:max_bars]
    
    if not filtered:
        logger.warning(f"No IFPs with occupancy >= {min_occupancy*100:.0f}%")
        return None
    
    labels = [x[0] for x in filtered]
    values = [x[1] * 100 for x in filtered]
    colors = [IFP_COLORS.get(_get_ifp_type(l), '#808080') for l in labels]
    
    fig = go.Figure(data=go.Bar(
        y=labels,
        x=values,
        orientation='h',
        marker_color=colors,
        text=[f'{v:.0f}%' for v in values],
        textposition='outside',
        hovertemplate='<b>%{y}</b><br>Occupancy: %{x:.1f}%<extra></extra>'
    ))
    
    fig.update_layout(
        title=title or 'IFP Occupancy (Interactive)',
        xaxis_title='Occupancy (%)',
        yaxis=dict(autorange="reversed"),
        template='plotly_white',
        height=max(400, len(labels) * 25 + 100),
        xaxis=dict(range=[0, 105]),
    )
    
    if save_path:
        save_path = Path(save_path)
        if not save_path.suffix:
            save_path = save_path.with_suffix('.html')
        fig.write_html(str(save_path), include_plotlyjs=True)
        logger.info(f"Saved interactive occupancy to {save_path}")
    
    return fig


def plot_residue_profile_html(
    ifp_result: "IFPTrajectoryResult",
    min_occupancy: float = 0.1,
    max_residues: int = 25,
    title: Optional[str] = None,
    save_path: Optional[Union[str, Path]] = None,
) -> "go.Figure":
    """
    Create interactive per-residue IFP profile with Plotly.
    
    Parameters
    ----------
    ifp_result : IFPTrajectoryResult
    min_occupancy : float
    max_residues : int
    title : str
    save_path : Path
        
    Returns
    -------
    go.Figure
    """
    _check_plotly()
    
    # Group by residue
    residue_data = {}
    for col, occ in ifp_result.ifp_occupancies.items():
        ifp_type = _get_ifp_type(col)
        residue = _get_residue(col)
        
        if residue not in residue_data:
            residue_data[residue] = {}
        residue_data[residue][ifp_type] = occ * 100
    
    # Filter and sort
    residue_totals = {r: sum(d.values()) for r, d in residue_data.items()}
    filtered = [(r, residue_totals[r]) for r in residue_data 
                if residue_totals[r] >= min_occupancy * 100]
    filtered.sort(key=lambda x: -x[1])
    residues = [x[0] for x in filtered[:max_residues]]
    
    if not residues:
        logger.warning("No residues pass filter")
        return None
    
    # Determine which types are present
    all_types = set()
    for r in residues:
        all_types.update(residue_data[r].keys())
    ordered_types = [t for t in IFP_ORDER if t in all_types]
    
    fig = go.Figure()
    
    for ifp_type in ordered_types:
        values = [residue_data[r].get(ifp_type, 0) for r in residues]
        fig.add_trace(go.Bar(
            name=IFP_LABELS.get(ifp_type, ifp_type),
            x=residues,
            y=values,
            marker_color=IFP_COLORS.get(ifp_type, '#808080'),
            hovertemplate=f'{IFP_LABELS.get(ifp_type, ifp_type)}<br>' +
                          '<b>%{x}</b><br>Occupancy: %{y:.1f}%<extra></extra>'
        ))
    
    fig.update_layout(
        title=title or 'Residue IFP Profile (Interactive)',
        xaxis_title='Residue',
        yaxis_title='Occupancy (%)',
        barmode='stack',
        template='plotly_white',
        legend=dict(yanchor="top", y=0.99, xanchor="right", x=0.99),
    )
    
    if save_path:
        save_path = Path(save_path)
        if not save_path.suffix:
            save_path = save_path.with_suffix('.html')
        fig.write_html(str(save_path), include_plotlyjs=True)
        logger.info(f"Saved interactive residue profile to {save_path}")
    
    return fig


def create_ifp_dashboard_html(
    ifp_result: "IFPTrajectoryResult",
    title: Optional[str] = None,
    save_path: Optional[Union[str, Path]] = None,
) -> "go.Figure":
    """
    Create comprehensive interactive IFP dashboard with Plotly.
    
    Combines: heatmap, timeseries, occupancy bar chart, and residue profile
    in a single interactive HTML file.
    
    Parameters
    ----------
    ifp_result : IFPTrajectoryResult
    title : str
        Dashboard title
    save_path : Path
        HTML output path
        
    Returns
    -------
    go.Figure
        Plotly figure with subplots
    """
    _check_plotly()
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            'A) IFP Heatmap', 
            'B) Time Series',
            'C) Occupancy', 
            'D) Residue Profile'
        ),
        specs=[
            [{"type": "heatmap"}, {"type": "scatter"}],
            [{"type": "bar"}, {"type": "bar"}]
        ],
        vertical_spacing=0.12,
        horizontal_spacing=0.1,
    )
    
    # Get data
    df = ifp_result.ifp_dataframe.copy()
    time = ifp_result.time_ps / 1000  # ns
    occupancies = ifp_result.ifp_occupancies
    
    # Filter and sort IFP columns
    ifp_cols = [c for c in ifp_result.ifp_columns 
                if c in df.columns and occupancies.get(c, 0) >= 0.05]
    ifp_cols = _sort_ifp_columns(ifp_cols, occupancies)[:30]
    
    # A) Heatmap
    if ifp_cols:
        matrix = df[ifp_cols].values.T
        fig.add_trace(
            go.Heatmap(
                z=matrix,
                x=time,
                y=ifp_cols,
                colorscale='YlGnBu',
                zmin=0, zmax=1,
                showscale=True,
                colorbar=dict(x=0.45, len=0.4, y=0.8),
                hovertemplate='<b>%{y}</b><br>Time: %{x:.2f} ns<br>Value: %{z}<extra></extra>',
            ),
            row=1, col=1
        )
    
    # B) Time series
    detected = set(_get_ifp_type(c) for c in ifp_result.ifp_columns)
    ifp_types = [t for t in IFP_ORDER if t in detected]
    
    for ifp_type in ifp_types:
        counts = ifp_result.get_ifp_timeseries(ifp_type)
        if counts.max() == 0:
            continue
        
        # Smooth
        counts = pd.Series(counts).rolling(window=10, center=True, min_periods=1).mean().values
        
        fig.add_trace(
            go.Scatter(
                x=time,
                y=counts,
                mode='lines',
                name=IFP_LABELS.get(ifp_type, ifp_type),
                line=dict(color=IFP_COLORS.get(ifp_type, '#808080'), width=1.5),
                legendgroup=ifp_type,
                showlegend=True,
            ),
            row=1, col=2
        )
    
    # C) Occupancy
    filtered_occ = [(k, v) for k, v in occupancies.items() if v >= 0.05]
    filtered_occ.sort(key=lambda x: -x[1])
    filtered_occ = filtered_occ[:15]
    
    if filtered_occ:
        labels_c = [x[0] for x in filtered_occ]
        values_c = [x[1] * 100 for x in filtered_occ]
        colors_c = [IFP_COLORS.get(_get_ifp_type(l), '#808080') for l in labels_c]
        
        fig.add_trace(
            go.Bar(
                y=labels_c,
                x=values_c,
                orientation='h',
                marker_color=colors_c,
                showlegend=False,
                hovertemplate='<b>%{y}</b><br>Occupancy: %{x:.1f}%<extra></extra>',
            ),
            row=2, col=1
        )
    
    # D) Residue profile
    residue_data = {}
    for col, occ in occupancies.items():
        ifp_type = _get_ifp_type(col)
        residue = _get_residue(col)
        if residue not in residue_data:
            residue_data[residue] = {}
        residue_data[residue][ifp_type] = occ * 100
    
    residue_totals = {r: sum(d.values()) for r, d in residue_data.items()}
    filtered_res = [(r, residue_totals[r]) for r in residue_data if residue_totals[r] >= 10]
    filtered_res.sort(key=lambda x: -x[1])
    residues = [x[0] for x in filtered_res[:12]]
    
    if residues:
        all_types = set()
        for r in residues:
            all_types.update(residue_data[r].keys())
        ordered_types = [t for t in IFP_ORDER if t in all_types]
        
        for ifp_type in ordered_types:
            values_d = [residue_data[r].get(ifp_type, 0) for r in residues]
            fig.add_trace(
                go.Bar(
                    name=IFP_LABELS.get(ifp_type, ifp_type),
                    x=residues,
                    y=values_d,
                    marker_color=IFP_COLORS.get(ifp_type, '#808080'),
                    showlegend=False,
                ),
                row=2, col=2
            )
    
    # Update layout
    fig.update_layout(
        title=dict(text=title or 'IFP Analysis Dashboard', x=0.5, font=dict(size=16)),
        height=900,
        template='plotly_white',
        barmode='stack',
        legend=dict(yanchor="top", y=0.95, xanchor="right", x=1.15),
    )
    
    # Axis labels
    fig.update_xaxes(title_text='Time (ns)', row=1, col=1)
    fig.update_yaxes(title_text='IFP', row=1, col=1)
    fig.update_xaxes(title_text='Time (ns)', row=1, col=2)
    fig.update_yaxes(title_text='Count', row=1, col=2)
    fig.update_xaxes(title_text='Occupancy (%)', row=2, col=1)
    fig.update_yaxes(autorange="reversed", row=2, col=1)
    fig.update_xaxes(title_text='Residue', row=2, col=2)
    fig.update_yaxes(title_text='Occupancy (%)', row=2, col=2)
    
    # Save
    if save_path:
        save_path = Path(save_path)
        if not save_path.suffix:
            save_path = save_path.with_suffix('.html')
        fig.write_html(str(save_path), include_plotlyjs=True)
        logger.info(f"Saved interactive dashboard to {save_path}")
    
    return fig


def generate_all_plots_html(
    ifp_result: "IFPTrajectoryResult",
    output_dir: Union[str, Path],
    prefix: str = "ifp",
    min_occupancy: float = 0.05,
) -> Dict[str, Path]:
    """
    Generate all Plotly HTML interactive visualizations.
    
    Parameters
    ----------
    ifp_result : IFPTrajectoryResult
    output_dir : Path
    prefix : str
    min_occupancy : float
    
    Returns
    -------
    dict
        {plot_type: saved_path}
    """
    _check_plotly()
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    results = {}
    
    # Individual plots
    plots = [
        ('heatmap', plot_ifp_heatmap_html, {'min_occupancy': min_occupancy}),
        ('timeseries', plot_ifp_timeseries_html, {}),
        ('occupancy', plot_ifp_occupancy_html, {'min_occupancy': min_occupancy}),
        ('residue_profile', plot_residue_profile_html, {'min_occupancy': min_occupancy}),
        ('dashboard', create_ifp_dashboard_html, {}),
    ]
    
    for name, func, kwargs in plots:
        try:
            save_path = output_dir / f"{prefix}_{name}.html"
            fig = func(ifp_result, save_path=save_path, **kwargs)
            if fig:
                results[name] = save_path
        except Exception as e:
            logger.warning(f"Failed to create {name} HTML: {e}")
    
    logger.info(f"Generated {len(results)} HTML plots in {output_dir}")
    return results


# ==============================================================================
# MM-GBSA VISUALIZATION FUNCTIONS
# ==============================================================================

# MM-GBSA Color scheme
MMGBSA_COLORS = {
    'total': '#2E86AB',      # Blue - total binding energy
    'vdw': '#A23B72',        # Magenta - van der Waals
    'elec': '#F18F01',       # Orange - electrostatic
    'gb': '#C73E1D',         # Red - GB solvation
    'favorable': '#2E7D32',  # Green - favorable contributions
    'unfavorable': '#C62828', # Red - unfavorable contributions
    'neutral': '#757575',    # Gray - near-zero contributions
}


def plot_mmgbsa_energy_evolution(
    frames: List[int],
    delta_g_series: List[float],
    system_name: str = "System",
    save_path: Optional[Union[str, Path]] = None,
    figsize: Tuple[float, float] = (12, 6),
    show_rolling_avg: bool = True,
    window: int = 50
) -> Optional[plt.Figure]:
    """
    Plot MM-GBSA binding energy evolution over trajectory frames.
    
    Parameters
    ----------
    frames : List[int]
        Frame numbers
    delta_g_series : List[float]
        ΔG values for each frame (kcal/mol)
    system_name : str
        System name for title
    save_path : Path, optional
        Path to save figure
    figsize : tuple
        Figure size
    show_rolling_avg : bool
        Whether to show rolling average
    window : int
        Rolling average window size
        
    Returns
    -------
    plt.Figure or None
    """
    _check_matplotlib()
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Convert to numpy for calculations
    frames = np.array(frames)
    delta_g = np.array(delta_g_series)
    
    # Plot raw data
    ax.plot(frames, delta_g, 'o', alpha=0.3, markersize=3, 
            color=MMGBSA_COLORS['total'], label='Per-frame ΔG')
    
    # Rolling average
    if show_rolling_avg and len(delta_g) >= window:
        rolling_avg = pd.Series(delta_g).rolling(window=window, center=True).mean()
        ax.plot(frames, rolling_avg, '-', linewidth=2.5, 
                color='#1A535C', label=f'Rolling avg (n={window})')
    
    # Mean and std
    mean_dg = np.mean(delta_g)
    std_dg = np.std(delta_g)
    
    ax.axhline(mean_dg, color='red', linestyle='--', linewidth=1.5, 
               label=f'Mean: {mean_dg:.2f} ± {std_dg:.2f} kcal/mol')
    ax.fill_between(frames, mean_dg - std_dg, mean_dg + std_dg, 
                    alpha=0.2, color='red')
    
    # Styling
    ax.set_xlabel('Frame', fontsize=12, fontweight='bold')
    ax.set_ylabel('ΔG (kcal/mol)', fontsize=12, fontweight='bold')
    ax.set_title(f'MM-GBSA Binding Energy Evolution: {system_name}', 
                 fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Add convergence indicator
    if len(delta_g) >= window * 2:
        first_half_mean = np.mean(delta_g[:len(delta_g)//2])
        second_half_mean = np.mean(delta_g[len(delta_g)//2:])
        drift = second_half_mean - first_half_mean
        convergence_text = f'Drift: {drift:+.2f} kcal/mol'
        ax.text(0.02, 0.98, convergence_text, transform=ax.transAxes, 
                fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved MM-GBSA energy evolution to {save_path}")
    
    return fig


def plot_mmgbsa_per_residue(
    per_residue: Dict[str, Dict[str, float]],
    system_name: str = "System",
    save_path: Optional[Union[str, Path]] = None,
    top_n: int = 40,
    figsize: Tuple[float, float] = (14, 8),
    show_components: bool = True
) -> Optional[plt.Figure]:
    """
    Plot per-residue MM-GBSA energy decomposition.
    
    Parameters
    ----------
    per_residue : Dict[str, Dict[str, float]]
        Dictionary mapping residue names to energy components.
        Each residue has keys: 'total', 'vdw', 'elec', 'gb'
    system_name : str
        System name for title
    save_path : Path, optional
        Path to save figure
    top_n : int
        Number of top contributing residues to show
    figsize : tuple
        Figure size
    show_components : bool
        Whether to show individual energy components
        
    Returns
    -------
    plt.Figure or None
    """
    _check_matplotlib()
    
    if not per_residue:
        logger.warning("No per-residue data provided")
        return None
    
    # Sort by absolute total contribution
    sorted_residues = sorted(per_residue.items(), 
                            key=lambda x: abs(x[1].get('total', 0)), 
                            reverse=True)[:top_n]
    
    residue_names = [r[0] for r in sorted_residues]
    totals = [r[1].get('total', 0) for r in sorted_residues]
    
    if show_components:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, 
                                        gridspec_kw={'height_ratios': [2, 1]})
        
        # Top panel: Total contributions (bar chart)
        colors = [MMGBSA_COLORS['favorable'] if t < 0 else MMGBSA_COLORS['unfavorable'] 
                  for t in totals]
        bars = ax1.bar(range(len(residue_names)), totals, color=colors, edgecolor='black', 
                       linewidth=0.5)
        
        ax1.set_xticks(range(len(residue_names)))
        ax1.set_xticklabels(residue_names, rotation=45, ha='right', fontsize=9)
        ax1.set_ylabel('ΔG contribution (kcal/mol)', fontsize=11, fontweight='bold')
        ax1.set_title(f'Per-Residue MM-GBSA Decomposition: {system_name}\n(Top {top_n} contributors)', 
                     fontsize=13, fontweight='bold')
        ax1.axhline(0, color='black', linewidth=1)
        ax1.grid(axis='y', alpha=0.3)
        
        # Add value labels on significant bars
        for i, (bar, val) in enumerate(zip(bars, totals)):
            if abs(val) > 0.5:
                y_pos = val + 0.1 if val > 0 else val - 0.3
                ax1.text(i, y_pos, f'{val:.1f}', ha='center', va='bottom' if val > 0 else 'top',
                        fontsize=7, fontweight='bold')
        
        # Bottom panel: Stacked components
        vdw = [r[1].get('vdw', 0) for r in sorted_residues]
        elec = [r[1].get('elec', 0) for r in sorted_residues]
        gb = [r[1].get('gb', 0) for r in sorted_residues]
        
        x = np.arange(len(residue_names))
        width = 0.25
        
        ax2.bar(x - width, vdw, width, label='VDW', color=MMGBSA_COLORS['vdw'], alpha=0.8)
        ax2.bar(x, elec, width, label='Electrostatic', color=MMGBSA_COLORS['elec'], alpha=0.8)
        ax2.bar(x + width, gb, width, label='GB Solvation', color=MMGBSA_COLORS['gb'], alpha=0.8)
        
        ax2.set_xticks(x)
        ax2.set_xticklabels(residue_names, rotation=45, ha='right', fontsize=9)
        ax2.set_ylabel('Energy (kcal/mol)', fontsize=11, fontweight='bold')
        ax2.set_xlabel('Residue', fontsize=11, fontweight='bold')
        ax2.legend(loc='upper right', fontsize=9)
        ax2.axhline(0, color='black', linewidth=1)
        ax2.grid(axis='y', alpha=0.3)
        
    else:
        fig, ax = plt.subplots(figsize=figsize)
        
        colors = [MMGBSA_COLORS['favorable'] if t < 0 else MMGBSA_COLORS['unfavorable'] 
                  for t in totals]
        ax.bar(range(len(residue_names)), totals, color=colors, edgecolor='black')
        
        ax.set_xticks(range(len(residue_names)))
        ax.set_xticklabels(residue_names, rotation=45, ha='right', fontsize=9)
        ax.set_ylabel('ΔG contribution (kcal/mol)', fontsize=12, fontweight='bold')
        ax.set_xlabel('Residue', fontsize=12, fontweight='bold')
        ax.set_title(f'Per-Residue MM-GBSA Decomposition: {system_name}', 
                     fontsize=14, fontweight='bold')
        ax.axhline(0, color='black', linewidth=1)
        ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved per-residue decomposition to {save_path}")
    
    return fig


def plot_mmgbsa_components_summary(
    vdw_mean: float,
    elec_mean: float,
    gb_mean: float,
    delta_g_mean: float,
    system_name: str = "System",
    save_path: Optional[Union[str, Path]] = None,
    figsize: Tuple[float, float] = (10, 6)
) -> Optional[plt.Figure]:
    """
    Plot MM-GBSA energy components summary.
    
    Parameters
    ----------
    vdw_mean : float
        Mean VDW contribution (kcal/mol)
    elec_mean : float
        Mean electrostatic contribution (kcal/mol)
    gb_mean : float
        Mean GB solvation contribution (kcal/mol)
    delta_g_mean : float
        Mean total ΔG (kcal/mol)
    system_name : str
        System name for title
    save_path : Path, optional
        Path to save figure
    figsize : tuple
        Figure size
        
    Returns
    -------
    plt.Figure or None
    """
    _check_matplotlib()
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Left: Bar chart of components
    components = ['VDW', 'Electrostatic', 'GB Solvation', 'Total ΔG']
    values = [vdw_mean, elec_mean, gb_mean, delta_g_mean]
    colors = [MMGBSA_COLORS['vdw'], MMGBSA_COLORS['elec'], 
              MMGBSA_COLORS['gb'], MMGBSA_COLORS['total']]
    
    bars = ax1.barh(components, values, color=colors, edgecolor='black', height=0.6)
    
    # Add value labels
    for bar, val in zip(bars, values):
        x_pos = val + 0.5 if val > 0 else val - 0.5
        ha = 'left' if val > 0 else 'right'
        ax1.text(x_pos, bar.get_y() + bar.get_height()/2, f'{val:.2f}', 
                va='center', ha=ha, fontsize=11, fontweight='bold')
    
    ax1.axvline(0, color='black', linewidth=1)
    ax1.set_xlabel('Energy (kcal/mol)', fontsize=12, fontweight='bold')
    ax1.set_title(f'Energy Components: {system_name}', fontsize=13, fontweight='bold')
    ax1.grid(axis='x', alpha=0.3)
    
    # Right: Pie chart of absolute contributions
    abs_values = [abs(v) for v in [vdw_mean, elec_mean, gb_mean]]
    labels = ['VDW', 'Electrostatic', 'GB Solvation']
    pie_colors = [MMGBSA_COLORS['vdw'], MMGBSA_COLORS['elec'], MMGBSA_COLORS['gb']]
    
    wedges, texts, autotexts = ax2.pie(abs_values, labels=labels, colors=pie_colors,
                                        autopct='%1.1f%%', startangle=90,
                                        explode=(0.02, 0.02, 0.02))
    for autotext in autotexts:
        autotext.set_fontsize(10)
        autotext.set_fontweight('bold')
    
    ax2.set_title('Relative Contribution (|value|)', fontsize=13, fontweight='bold')
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved energy components summary to {save_path}")
    
    return fig


def plot_mmgbsa_hotspots(
    per_residue: Dict[str, Dict[str, float]],
    threshold: float = -1.0,
    system_name: str = "System",
    save_path: Optional[Union[str, Path]] = None,
    figsize: Tuple[float, float] = (12, 8)
) -> Optional[plt.Figure]:
    """
    Plot binding hotspots (residues with significant favorable contributions).
    
    Parameters
    ----------
    per_residue : Dict[str, Dict[str, float]]
        Per-residue energy decomposition
    threshold : float
        Energy threshold for hotspot (kcal/mol, default -1.0)
    system_name : str
        System name for title
    save_path : Path, optional
        Path to save figure
    figsize : tuple
        Figure size
        
    Returns
    -------
    plt.Figure or None
    """
    _check_matplotlib()
    
    if not per_residue:
        logger.warning("No per-residue data provided")
        return None
    
    # Filter hotspots (favorable contributors)
    hotspots = {k: v for k, v in per_residue.items() 
                if v.get('total', 0) < threshold}
    
    if not hotspots:
        logger.warning(f"No hotspots found with threshold < {threshold}")
        return None
    
    # Sort by contribution (most favorable first)
    sorted_hotspots = sorted(hotspots.items(), 
                            key=lambda x: x[1].get('total', 0))
    
    fig, ax = plt.subplots(figsize=figsize)
    
    residue_names = [r[0] for r in sorted_hotspots]
    totals = [r[1].get('total', 0) for r in sorted_hotspots]
    
    # Gradient coloring based on contribution
    norm = plt.Normalize(min(totals), threshold)
    colors = plt.cm.Greens_r(norm(totals))
    
    bars = ax.barh(range(len(residue_names)), totals, color=colors, 
                   edgecolor='black', linewidth=0.5)
    
    ax.set_yticks(range(len(residue_names)))
    ax.set_yticklabels(residue_names, fontsize=10)
    ax.set_xlabel('ΔG contribution (kcal/mol)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Residue', fontsize=12, fontweight='bold')
    ax.set_title(f'Binding Hotspots: {system_name}\n(ΔG < {threshold} kcal/mol)', 
                 fontsize=14, fontweight='bold')
    ax.axvline(0, color='black', linewidth=1)
    ax.grid(axis='x', alpha=0.3)
    
    # Add value labels
    for bar, val in zip(bars, totals):
        ax.text(val - 0.2, bar.get_y() + bar.get_height()/2, f'{val:.2f}', 
                va='center', ha='right', fontsize=9, fontweight='bold', color='white')
    
    # Summary text
    total_hotspot = sum(totals)
    summary = f'Total hotspot contribution: {total_hotspot:.2f} kcal/mol\nHotspots identified: {len(hotspots)}'
    ax.text(0.98, 0.02, summary, transform=ax.transAxes, fontsize=10,
            verticalalignment='bottom', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        logger.info(f"Saved hotspots plot to {save_path}")
    
    return fig


# ==============================================================================
# MM-GBSA PLOTLY HTML FUNCTIONS
# ==============================================================================

def plot_mmgbsa_energy_evolution_html(
    frames: List[int],
    delta_g_series: List[float],
    system_name: str = "System",
    save_path: Optional[Union[str, Path]] = None,
    show_rolling_avg: bool = True,
    window: int = 50
) -> Optional[go.Figure]:
    """
    Create interactive MM-GBSA energy evolution plot with Plotly.
    
    Parameters
    ----------
    frames : List[int]
        Frame numbers
    delta_g_series : List[float]
        ΔG values for each frame
    system_name : str
        System name for title
    save_path : Path, optional
        Path to save HTML file
    show_rolling_avg : bool
        Whether to show rolling average
    window : int
        Rolling average window size
        
    Returns
    -------
    go.Figure or None
    """
    _check_plotly()
    
    frames = np.array(frames)
    delta_g = np.array(delta_g_series)
    mean_dg = np.mean(delta_g)
    std_dg = np.std(delta_g)
    
    fig = go.Figure()
    
    # Per-frame values
    fig.add_trace(go.Scatter(
        x=frames, y=delta_g,
        mode='markers',
        name='Per-frame ΔG',
        marker=dict(color=MMGBSA_COLORS['total'], size=4, opacity=0.4),
        hovertemplate='Frame: %{x}<br>ΔG: %{y:.2f} kcal/mol<extra></extra>'
    ))
    
    # Rolling average
    if show_rolling_avg and len(delta_g) >= window:
        rolling_avg = pd.Series(delta_g).rolling(window=window, center=True).mean()
        fig.add_trace(go.Scatter(
            x=frames, y=rolling_avg,
            mode='lines',
            name=f'Rolling avg (n={window})',
            line=dict(color='#1A535C', width=3),
            hovertemplate='Frame: %{x}<br>Rolling Avg: %{y:.2f} kcal/mol<extra></extra>'
        ))
    
    # Mean line
    fig.add_trace(go.Scatter(
        x=[frames[0], frames[-1]],
        y=[mean_dg, mean_dg],
        mode='lines',
        name=f'Mean: {mean_dg:.2f} ± {std_dg:.2f}',
        line=dict(color='red', dash='dash', width=2)
    ))
    
    # Std deviation band
    fig.add_trace(go.Scatter(
        x=list(frames) + list(frames[::-1]),
        y=[mean_dg + std_dg] * len(frames) + [mean_dg - std_dg] * len(frames),
        fill='toself',
        fillcolor='rgba(255, 0, 0, 0.1)',
        line=dict(color='rgba(255, 0, 0, 0)'),
        name='±1 Std Dev',
        hoverinfo='skip'
    ))
    
    fig.update_layout(
        title=dict(text=f'MM-GBSA Binding Energy Evolution: {system_name}', font=dict(size=16)),
        xaxis_title='Frame',
        yaxis_title='ΔG (kcal/mol)',
        hovermode='x unified',
        template='plotly_white',
        legend=dict(x=0.02, y=0.98, bgcolor='rgba(255,255,255,0.8)')
    )
    
    if save_path:
        fig.write_html(str(save_path))
        logger.info(f"Saved interactive energy evolution to {save_path}")
    
    return fig


def plot_mmgbsa_per_residue_html(
    per_residue: Dict[str, Dict[str, float]],
    system_name: str = "System",
    save_path: Optional[Union[str, Path]] = None,
    top_n: int = 50
) -> Optional[go.Figure]:
    """
    Create interactive per-residue MM-GBSA decomposition with Plotly.
    
    Parameters
    ----------
    per_residue : Dict[str, Dict[str, float]]
        Per-residue energy decomposition
    system_name : str
        System name for title
    save_path : Path, optional
        Path to save HTML file
    top_n : int
        Number of top contributing residues to show
        
    Returns
    -------
    go.Figure or None
    """
    _check_plotly()
    
    if not per_residue:
        logger.warning("No per-residue data provided")
        return None
    
    # Sort by absolute total contribution
    sorted_residues = sorted(per_residue.items(), 
                            key=lambda x: abs(x[1].get('total', 0)), 
                            reverse=True)[:top_n]
    
    residue_names = [r[0] for r in sorted_residues]
    totals = [r[1].get('total', 0) for r in sorted_residues]
    vdw = [r[1].get('vdw', 0) for r in sorted_residues]
    elec = [r[1].get('elec', 0) for r in sorted_residues]
    gb = [r[1].get('gb', 0) for r in sorted_residues]
    
    # Color based on favorable/unfavorable
    colors = [MMGBSA_COLORS['favorable'] if t < 0 else MMGBSA_COLORS['unfavorable'] 
              for t in totals]
    
    fig = make_subplots(rows=2, cols=1, 
                        subplot_titles=['Total ΔG Contribution', 'Energy Components'],
                        vertical_spacing=0.15)
    
    # Top: Total contributions
    fig.add_trace(go.Bar(
        x=residue_names, y=totals,
        marker_color=colors,
        name='Total ΔG',
        hovertemplate='%{x}<br>Total: %{y:.2f} kcal/mol<extra></extra>'
    ), row=1, col=1)
    
    # Bottom: Components grouped
    fig.add_trace(go.Bar(
        x=residue_names, y=vdw, name='VDW',
        marker_color=MMGBSA_COLORS['vdw'],
        hovertemplate='%{x}<br>VDW: %{y:.2f} kcal/mol<extra></extra>'
    ), row=2, col=1)
    
    fig.add_trace(go.Bar(
        x=residue_names, y=elec, name='Electrostatic',
        marker_color=MMGBSA_COLORS['elec'],
        hovertemplate='%{x}<br>Elec: %{y:.2f} kcal/mol<extra></extra>'
    ), row=2, col=1)
    
    fig.add_trace(go.Bar(
        x=residue_names, y=gb, name='GB Solvation',
        marker_color=MMGBSA_COLORS['gb'],
        hovertemplate='%{x}<br>GB: %{y:.2f} kcal/mol<extra></extra>'
    ), row=2, col=1)
    
    fig.update_layout(
        title=dict(text=f'Per-Residue MM-GBSA Decomposition: {system_name}', font=dict(size=16)),
        barmode='group',
        template='plotly_white',
        height=800,
        showlegend=True,
        legend=dict(x=1.02, y=0.5)
    )
    
    fig.update_xaxes(tickangle=45, tickfont=dict(size=9))
    fig.update_yaxes(title_text='ΔG (kcal/mol)', row=1, col=1)
    fig.update_yaxes(title_text='Energy (kcal/mol)', row=2, col=1)
    
    if save_path:
        fig.write_html(str(save_path))
        logger.info(f"Saved interactive per-residue decomposition to {save_path}")
    
    return fig


def plot_mmgbsa_components_html(
    vdw_mean: float,
    elec_mean: float,
    gb_mean: float,
    delta_g_mean: float,
    system_name: str = "System",
    save_path: Optional[Union[str, Path]] = None
) -> Optional[go.Figure]:
    """
    Create interactive MM-GBSA energy components visualization with Plotly.
    
    Parameters
    ----------
    vdw_mean : float
        Mean VDW contribution
    elec_mean : float
        Mean electrostatic contribution
    gb_mean : float
        Mean GB solvation contribution
    delta_g_mean : float
        Mean total ΔG
    system_name : str
        System name for title
    save_path : Path, optional
        Path to save HTML file
        
    Returns
    -------
    go.Figure or None
    """
    _check_plotly()
    
    fig = make_subplots(rows=1, cols=2, 
                        specs=[[{'type': 'bar'}, {'type': 'pie'}]],
                        subplot_titles=['Energy Components', 'Relative Contribution'])
    
    # Bar chart
    components = ['VDW', 'Electrostatic', 'GB Solvation', 'Total ΔG']
    values = [vdw_mean, elec_mean, gb_mean, delta_g_mean]
    colors = [MMGBSA_COLORS['vdw'], MMGBSA_COLORS['elec'], 
              MMGBSA_COLORS['gb'], MMGBSA_COLORS['total']]
    
    fig.add_trace(go.Bar(
        y=components, x=values,
        orientation='h',
        marker_color=colors,
        text=[f'{v:.2f}' for v in values],
        textposition='outside',
        hovertemplate='%{y}: %{x:.2f} kcal/mol<extra></extra>'
    ), row=1, col=1)
    
    # Pie chart
    abs_values = [abs(v) for v in [vdw_mean, elec_mean, gb_mean]]
    pie_labels = ['VDW', 'Electrostatic', 'GB Solvation']
    pie_colors = [MMGBSA_COLORS['vdw'], MMGBSA_COLORS['elec'], MMGBSA_COLORS['gb']]
    
    fig.add_trace(go.Pie(
        labels=pie_labels,
        values=abs_values,
        marker_colors=pie_colors,
        textinfo='label+percent',
        hovertemplate='%{label}: %{value:.2f} kcal/mol (%{percent})<extra></extra>'
    ), row=1, col=2)
    
    fig.update_layout(
        title=dict(text=f'MM-GBSA Energy Components: {system_name}', font=dict(size=16)),
        template='plotly_white',
        height=500,
        showlegend=False
    )
    
    fig.update_xaxes(title_text='Energy (kcal/mol)', row=1, col=1)
    
    if save_path:
        fig.write_html(str(save_path))
        logger.info(f"Saved interactive energy components to {save_path}")
    
    return fig


def create_mmgbsa_dashboard_html(
    frames: List[int],
    delta_g_series: List[float],
    per_residue: Dict[str, Dict[str, float]],
    vdw_mean: float,
    elec_mean: float,
    gb_mean: float,
    delta_g_mean: float,
    delta_g_std: float,
    system_name: str = "System",
    save_path: Optional[Union[str, Path]] = None,
    top_n: int = 30
) -> Optional[go.Figure]:
    """
    Create comprehensive MM-GBSA dashboard with all visualizations.
    
    Parameters
    ----------
    frames : List[int]
        Frame numbers
    delta_g_series : List[float]
        ΔG values for each frame
    per_residue : Dict[str, Dict[str, float]]
        Per-residue energy decomposition
    vdw_mean, elec_mean, gb_mean : float
        Mean energy components
    delta_g_mean, delta_g_std : float
        Mean and std of total ΔG
    system_name : str
        System name for title
    save_path : Path, optional
        Path to save HTML file
    top_n : int
        Number of top residues to show
        
    Returns
    -------
    go.Figure or None
    """
    _check_plotly()
    
    fig = make_subplots(
        rows=3, cols=2,
        specs=[
            [{'colspan': 2}, None],
            [{}, {}],
            [{'colspan': 2}, None]
        ],
        subplot_titles=[
            'Energy Evolution Over Trajectory',
            'Energy Components', 'Top Contributing Residues',
            'Per-Residue Decomposition'
        ],
        vertical_spacing=0.1,
        horizontal_spacing=0.1
    )
    
    frames = np.array(frames)
    delta_g = np.array(delta_g_series)
    
    # Row 1: Energy evolution
    fig.add_trace(go.Scatter(
        x=frames, y=delta_g,
        mode='markers',
        name='Per-frame ΔG',
        marker=dict(color=MMGBSA_COLORS['total'], size=3, opacity=0.4),
        hovertemplate='Frame: %{x}<br>ΔG: %{y:.2f}<extra></extra>'
    ), row=1, col=1)
    
    # Rolling average
    if len(delta_g) >= 50:
        rolling = pd.Series(delta_g).rolling(window=50, center=True).mean()
        fig.add_trace(go.Scatter(
            x=frames, y=rolling,
            mode='lines',
            name='Rolling Avg',
            line=dict(color='#1A535C', width=2)
        ), row=1, col=1)
    
    # Mean line
    fig.add_trace(go.Scatter(
        x=[frames[0], frames[-1]], y=[delta_g_mean, delta_g_mean],
        mode='lines',
        name=f'Mean: {delta_g_mean:.2f}±{delta_g_std:.2f}',
        line=dict(color='red', dash='dash')
    ), row=1, col=1)
    
    # Row 2 Left: Energy components bar
    components = ['VDW', 'Elec', 'GB', 'Total']
    values = [vdw_mean, elec_mean, gb_mean, delta_g_mean]
    colors = [MMGBSA_COLORS['vdw'], MMGBSA_COLORS['elec'], 
              MMGBSA_COLORS['gb'], MMGBSA_COLORS['total']]
    
    fig.add_trace(go.Bar(
        x=components, y=values,
        marker_color=colors,
        text=[f'{v:.1f}' for v in values],
        textposition='outside',
        showlegend=False,
        hovertemplate='%{x}: %{y:.2f} kcal/mol<extra></extra>'
    ), row=2, col=1)
    
    # Row 2 Right: Top residues
    if per_residue:
        sorted_res = sorted(per_residue.items(), 
                           key=lambda x: x[1].get('total', 0))[:15]
        res_names = [r[0] for r in sorted_res]
        res_totals = [r[1].get('total', 0) for r in sorted_res]
        
        fig.add_trace(go.Bar(
            y=res_names, x=res_totals,
            orientation='h',
            marker_color=MMGBSA_COLORS['favorable'],
            showlegend=False,
            hovertemplate='%{y}: %{x:.2f} kcal/mol<extra></extra>'
        ), row=2, col=2)
    
    # Row 3: Per-residue decomposition
    if per_residue:
        sorted_all = sorted(per_residue.items(), 
                           key=lambda x: abs(x[1].get('total', 0)), 
                           reverse=True)[:top_n]
        
        res_names = [r[0] for r in sorted_all]
        totals = [r[1].get('total', 0) for r in sorted_all]
        bar_colors = [MMGBSA_COLORS['favorable'] if t < 0 else MMGBSA_COLORS['unfavorable'] 
                      for t in totals]
        
        fig.add_trace(go.Bar(
            x=res_names, y=totals,
            marker_color=bar_colors,
            showlegend=False,
            hovertemplate='%{x}<br>ΔG: %{y:.2f} kcal/mol<extra></extra>'
        ), row=3, col=1)
    
    fig.update_layout(
        title=dict(
            text=f'MM-GBSA Analysis Dashboard: {system_name}',
            font=dict(size=18)
        ),
        template='plotly_white',
        height=1000,
        showlegend=True,
        legend=dict(x=1.02, y=0.95)
    )
    
    fig.update_yaxes(title_text='ΔG (kcal/mol)', row=1, col=1)
    fig.update_xaxes(title_text='Frame', row=1, col=1)
    fig.update_yaxes(title_text='Energy (kcal/mol)', row=2, col=1)
    fig.update_xaxes(title_text='ΔG (kcal/mol)', row=2, col=2)
    fig.update_xaxes(tickangle=45, tickfont=dict(size=8), row=3, col=1)
    fig.update_yaxes(title_text='ΔG (kcal/mol)', row=3, col=1)
    
    if save_path:
        fig.write_html(str(save_path))
        logger.info(f"Saved MM-GBSA dashboard to {save_path}")
    
    return fig


def generate_all_mmgbsa_plots(
    frames: List[int],
    delta_g_series: List[float],
    per_residue: Dict[str, Dict[str, float]],
    vdw_mean: float,
    elec_mean: float,
    gb_mean: float,
    delta_g_mean: float,
    delta_g_std: float,
    system_name: str = "System",
    output_dir: Union[str, Path] = "mmgbsa_figures",
    prefix: str = "mmgbsa"
) -> Dict[str, Path]:
    """
    Generate all MM-GBSA visualizations (both static and interactive).
    
    Parameters
    ----------
    frames : List[int]
        Frame numbers
    delta_g_series : List[float]
        ΔG values for each frame
    per_residue : Dict[str, Dict[str, float]]
        Per-residue energy decomposition
    vdw_mean, elec_mean, gb_mean : float
        Mean energy components
    delta_g_mean, delta_g_std : float
        Mean and std of total ΔG
    system_name : str
        System name for titles
    output_dir : Path
        Output directory for figures
    prefix : str
        Filename prefix
        
    Returns
    -------
    dict
        {plot_name: saved_path}
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    results = {}
    
    # Static matplotlib plots
    if MATPLOTLIB_AVAILABLE:
        try:
            path = output_dir / f"{prefix}_energy_evolution.png"
            plot_mmgbsa_energy_evolution(frames, delta_g_series, system_name, save_path=path)
            results['energy_evolution_png'] = path
        except Exception as e:
            logger.warning(f"Failed to create energy evolution plot: {e}")
        
        try:
            path = output_dir / f"{prefix}_per_residue.png"
            plot_mmgbsa_per_residue(per_residue, system_name, save_path=path)
            results['per_residue_png'] = path
        except Exception as e:
            logger.warning(f"Failed to create per-residue plot: {e}")
        
        try:
            path = output_dir / f"{prefix}_components.png"
            plot_mmgbsa_components_summary(vdw_mean, elec_mean, gb_mean, delta_g_mean, 
                                          system_name, save_path=path)
            results['components_png'] = path
        except Exception as e:
            logger.warning(f"Failed to create components plot: {e}")
        
        try:
            path = output_dir / f"{prefix}_hotspots.png"
            plot_mmgbsa_hotspots(per_residue, system_name=system_name, save_path=path)
            results['hotspots_png'] = path
        except Exception as e:
            logger.warning(f"Failed to create hotspots plot: {e}")
    
    # Interactive Plotly plots
    if PLOTLY_AVAILABLE:
        try:
            path = output_dir / f"{prefix}_energy_evolution.html"
            plot_mmgbsa_energy_evolution_html(frames, delta_g_series, system_name, save_path=path)
            results['energy_evolution_html'] = path
        except Exception as e:
            logger.warning(f"Failed to create energy evolution HTML: {e}")
        
        try:
            path = output_dir / f"{prefix}_per_residue.html"
            plot_mmgbsa_per_residue_html(per_residue, system_name, save_path=path)
            results['per_residue_html'] = path
        except Exception as e:
            logger.warning(f"Failed to create per-residue HTML: {e}")
        
        try:
            path = output_dir / f"{prefix}_components.html"
            plot_mmgbsa_components_html(vdw_mean, elec_mean, gb_mean, delta_g_mean, 
                                        system_name, save_path=path)
            results['components_html'] = path
        except Exception as e:
            logger.warning(f"Failed to create components HTML: {e}")
        
        try:
            path = output_dir / f"{prefix}_dashboard.html"
            create_mmgbsa_dashboard_html(
                frames, delta_g_series, per_residue,
                vdw_mean, elec_mean, gb_mean, delta_g_mean, delta_g_std,
                system_name, save_path=path
            )
            results['dashboard_html'] = path
        except Exception as e:
            logger.warning(f"Failed to create dashboard HTML: {e}")
    
    logger.info(f"Generated {len(results)} MM-GBSA plots in {output_dir}")
    return results
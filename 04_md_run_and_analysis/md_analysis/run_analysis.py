"""
CCT-OSR1 Binding Analysis - Master Analysis Runner
===================================================
Publication-Grade MD Trajectory Analysis Pipeline

Analyses included:
- RMSD/RMSF structural analysis
- Contact maps and Q native (Best-Hummer-Eaton)
- Binding distance monitoring
- Water interface and residence time
- Conformational clustering with PCA
- Convergence analysis (block averaging, autocorrelation)
"""

import sys
from pathlib import Path
import time
from datetime import datetime
import traceback

if __package__ is None or __package__ == "":
    sys.path.insert(0, str(Path(__file__).parent))

try:
    from .config import (
        PROJECT_ROOT, OUTPUT_DIR, HTML_DIR, FIGURES_DIR, DATA_OUTPUT_DIR,
        TOPOLOGY_PDB, TRAJECTORIES, SAVE_FREQ_PS
    )
except ImportError:  # pragma: no cover
    from config import (
        PROJECT_ROOT, OUTPUT_DIR, HTML_DIR, FIGURES_DIR, DATA_OUTPUT_DIR,
        TOPOLOGY_PDB, TRAJECTORIES, SAVE_FREQ_PS
    )

try:
    from .trajectory_loader import load_and_concatenate_trajectories, save_trajectory_info
    from .analysis_rmsd import run_rmsd_analysis
    from .analysis_contacts import run_contact_analysis
    from .analysis_binding import run_binding_analysis
    from .analysis_ifp import run_ifp_analysis
    from .analysis_water import run_water_analysis
    from .analysis_clustering import run_clustering_analysis
    from .analysis_convergence import run_convergence_analysis
    from .provenance import (
        make_analysis_provenance,
        make_run_provenance,
        write_provenance_json,
    )
except ImportError:  # pragma: no cover
    from trajectory_loader import load_and_concatenate_trajectories, save_trajectory_info
    from analysis_rmsd import run_rmsd_analysis
    from analysis_contacts import run_contact_analysis
    from analysis_binding import run_binding_analysis
    from analysis_ifp import run_ifp_analysis
    from analysis_water import run_water_analysis
    from analysis_clustering import run_clustering_analysis
    from analysis_convergence import run_convergence_analysis
    from provenance import (
        make_analysis_provenance,
        make_run_provenance,
        write_provenance_json,
    )


def create_dashboard_html(output_path: Path = None):
    """Create a master HTML dashboard linking all analyses."""
    if output_path is None:
        output_path = HTML_DIR / "00_dashboard.html"
    
    html_files = sorted(HTML_DIR.glob("*.html"))
    html_files = [f for f in html_files if f.name != "00_dashboard.html"]
    
    html_content = """
<!DOCTYPE html>
<html>
<head>
    <title>CCT-OSR1 Binding Analysis Dashboard</title>
    <style>
        body {
            font-family: 'Segoe UI', Arial, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f5f5f5;
        }
        .header {
            background: linear-gradient(135deg, #1f77b4, #ff7f0e);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
        }
        .header h1 {
            margin: 0;
            font-size: 2.5em;
        }
        .header p {
            margin: 10px 0 0 0;
            opacity: 0.9;
        }
        .grid {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(350px, 1fr));
            gap: 20px;
        }
        .card {
            background: white;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            overflow: hidden;
            transition: transform 0.2s, box-shadow 0.2s;
        }
        .card:hover {
            transform: translateY(-5px);
            box-shadow: 0 5px 20px rgba(0,0,0,0.15);
        }
        .card a {
            text-decoration: none;
            color: inherit;
            display: block;
        }
        .card-header {
            background: #1f77b4;
            color: white;
            padding: 15px 20px;
        }
        .card-header h3 {
            margin: 0;
            font-size: 1.1em;
        }
        .card-body {
            padding: 15px 20px;
        }
        .card-body p {
            color: #666;
            margin: 0;
            font-size: 0.9em;
        }
        .summary {
            background: white;
            border-radius: 10px;
            padding: 20px;
            margin-bottom: 30px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        .summary h2 {
            margin-top: 0;
            color: #1f77b4;
        }
        .summary-grid {
            display: grid;
            grid-template-columns: repeat(4, 1fr);
            gap: 15px;
        }
        .summary-item {
            text-align: center;
            padding: 15px;
            background: #f8f9fa;
            border-radius: 8px;
        }
        .summary-item .value {
            font-size: 1.8em;
            font-weight: bold;
            color: #1f77b4;
        }
        .summary-item .label {
            color: #666;
            font-size: 0.85em;
        }
        .footer {
            text-align: center;
            color: #666;
            margin-top: 30px;
            padding: 20px;
        }
        .category {
            margin-bottom: 30px;
        }
        .category h2 {
            color: #333;
            border-bottom: 2px solid #1f77b4;
            padding-bottom: 10px;
        }
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 CCT-OSR1 / WNK1 Binding Analysis</h1>
        <p>Molecular Dynamics Analysis Dashboard | Generated: TIMESTAMP</p>
    </div>
    
    <div class="summary">
        <h2>📊 Analysis Summary</h2>
        <div class="summary-grid">
            <div class="summary-item">
                <div class="value">96</div>
                <div class="label">Protein Residues</div>
            </div>
            <div class="summary-item">
                <div class="value">17</div>
                <div class="label">Peptide Residues</div>
            </div>
            <div class="summary-item">
                <div class="value">RFxV</div>
                <div class="label">Binding Motif</div>
            </div>
            <div class="summary-item">
                <div class="value">N_FRAMES</div>
                <div class="label">Total Frames</div>
            </div>
        </div>
    </div>
    
    <div class="category">
        <h2>📈 Structural Analysis</h2>
        <div class="grid">
            <div class="card">
                <a href="01_rmsd_timeseries.html" target="_blank">
                    <div class="card-header"><h3>RMSD Time Series</h3></div>
                    <div class="card-body"><p>Root Mean Square Deviation of protein and peptide over simulation time.</p></div>
                </a>
            </div>
            <div class="card">
                <a href="02_rmsf.html" target="_blank">
                    <div class="card-header"><h3>RMSF by Residue</h3></div>
                    <div class="card-body"><p>Root Mean Square Fluctuation per residue, identifying flexible regions.</p></div>
                </a>
            </div>
        </div>
    </div>

    <div class="category">
        <h2>🧷 IFP (MD-IFP)</h2>
        <div class="grid">
            <div class="card">
                <a href="ifp_heatmap.html" target="_blank">
                    <div class="card-header"><h3>IFP Heatmap</h3></div>
                    <div class="card-body"><p>Binary interaction fingerprints over time (MD-IFP style).</p></div>
                </a>
            </div>
            <div class="card">
                <a href="ifp_timeseries.html" target="_blank">
                    <div class="card-header"><h3>IFP Time Series</h3></div>
                    <div class="card-body"><p>Counts of interaction types per frame.</p></div>
                </a>
            </div>
            <div class="card">
                <a href="ifp_occupancy.html" target="_blank">
                    <div class="card-header"><h3>IFP Occupancy</h3></div>
                    <div class="card-body"><p>Top interaction occupancies across the trajectory.</p></div>
                </a>
            </div>
            <div class="card">
                <a href="ifp_residue_summary.html" target="_blank">
                    <div class="card-header"><h3>IFP Residue Summary</h3></div>
                    <div class="card-body"><p>Residue-level breakdown of interaction types.</p></div>
                </a>
            </div>
        </div>
    </div>
    
    <div class="category">
        <h2>🔗 Contact Analysis</h2>
        <div class="grid">
            <div class="card">
                <a href="03_contact_map.html" target="_blank">
                    <div class="card-header"><h3>Contact Frequency Map</h3></div>
                    <div class="card-body"><p>Residue-residue contact frequencies between protein and peptide.</p></div>
                </a>
            </div>
            <div class="card">
                <a href="04_contact_evolution.html" target="_blank">
                    <div class="card-header"><h3>Contact Evolution</h3></div>
                    <div class="card-body"><p>Temporal evolution of peptide contacts with protein.</p></div>
                </a>
            </div>
            <div class="card">
                <a href="05_total_contacts.html" target="_blank">
                    <div class="card-header"><h3>Total Contacts</h3></div>
                    <div class="card-body"><p>Total number of protein-peptide contacts over time.</p></div>
                </a>
            </div>
        </div>
    </div>
    
    <div class="category">
        <h2>🎯 Binding Analysis</h2>
        <div class="grid">
            <div class="card">
                <a href="06_binding_distances.html" target="_blank">
                    <div class="card-header"><h3>Binding Distances</h3></div>
                    <div class="card-body"><p>COM-COM, RFxV, and Phe insertion distances tracking binding events.</p></div>
                </a>
            </div>
            <div class="card">
                <a href="07_distance_distributions.html" target="_blank">
                    <div class="card-header"><h3>Distance Distributions</h3></div>
                    <div class="card-body"><p>Distribution of key binding distances.</p></div>
                </a>
            </div>
        </div>
    </div>
    
    <div class="category">
        <h2>💧 Solvent Analysis</h2>
        <div class="grid">
            <div class="card">
                <a href="08_water_interface.html" target="_blank">
                    <div class="card-header"><h3>Interface Waters</h3></div>
                    <div class="card-body"><p>Water molecules at the protein-peptide interface and around RFxV motif.</p></div>
                </a>
            </div>
            <div class="card">
                <a href="09_water_residence.html" target="_blank">
                    <div class="card-header"><h3>Water Residence Time</h3></div>
                    <div class="card-body"><p>Residence time of water at the binding interface (intermittent correlation).</p></div>
                </a>
            </div>
        </div>
    </div>
    
    <div class="category">
        <h2>📉 Q Native Contacts</h2>
        <div class="grid">
            <div class="card">
                <a href="11_q_native.html" target="_blank">
                    <div class="card-header"><h3>Q Native Time Series</h3></div>
                    <div class="card-body"><p>Fraction of native contacts (Best-Hummer-Eaton switching function, β=5Å, λ=1.8).</p></div>
                </a>
            </div>
        </div>
    </div>
    
    <div class="category">
        <h2>🔮 Conformational Analysis</h2>
        <div class="grid">
            <div class="card">
                <a href="12_clustering.html" target="_blank">
                    <div class="card-header"><h3>Conformational Clustering</h3></div>
                    <div class="card-body"><p>RMSD-based clustering of peptide conformations.</p></div>
                </a>
            </div>
            <div class="card">
                <a href="13_pca.html" target="_blank">
                    <div class="card-header"><h3>PCA Projection</h3></div>
                    <div class="card-body"><p>Principal Component Analysis of peptide conformations.</p></div>
                </a>
            </div>
        </div>
    </div>
    
    <div class="category">
        <h2>📊 Convergence Analysis</h2>
        <div class="grid">
            <div class="card">
                <a href="14_convergence_dashboard.html" target="_blank">
                    <div class="card-header"><h3>Convergence Dashboard</h3></div>
                    <div class="card-body"><p>Block averaging, autocorrelation time, and effective sample size.</p></div>
                </a>
            </div>
            <div class="card">
                <a href="15_block_averaging.html" target="_blank">
                    <div class="card-header"><h3>Block Averaging</h3></div>
                    <div class="card-body"><p>Block averages with SEM for RMSD, contacts, and binding distances.</p></div>
                </a>
            </div>
        </div>
    </div>
    
    <div class="footer">
        <p>CCT-OSR1 Binding Analysis | WNK1-OSR1 Signaling Project</p>
        <p>Analysis framework for MD trajectory exploration</p>
    </div>
</body>
</html>
""".replace("TIMESTAMP", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"Dashboard saved: {output_path}")


def run_all_analyses(stride: int = 1, skip_water: bool = False, skip_convergence: bool = False):
    """
    Run all analyses sequentially.
    
    Parameters
    ----------
    stride : int
        Frame stride for trajectory loading
    skip_water : bool
        Skip water analysis (computationally expensive)
    skip_convergence : bool
        Skip convergence analysis
    """
    start_time = time.time()
    provenance_dir = OUTPUT_DIR / "provenance"
    
    print("=" * 70)
    print("CCT-OSR1 BINDING ANALYSIS - PUBLICATION-GRADE PIPELINE")
    print("=" * 70)
    print(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Stride: {stride}")
    print(f"Output directory: {PROJECT_ROOT}")
    print()
    print("Analyses enabled:")
    print("  [OK] RMSD/RMSF structural analysis")
    print("  [OK] IFP (MD-IFP) interaction fingerprints")
    print("  [OK] Contact maps + Q native (Best-Hummer-Eaton)")
    print("  [OK] Binding distance monitoring")
    print("  [OK] Conformational clustering + PCA")
    print(f"  {'[OK]' if not skip_water else '[SKIP]'} Water interface + residence time")
    print(f"  {'[OK]' if not skip_convergence else '[SKIP]'} Convergence analysis")
    print()

    # Run-level provenance (best-effort, never fails the run)
    run_prov = make_run_provenance(
        name="run_all_analyses",
        repo_root=PROJECT_ROOT,
        inputs=[TOPOLOGY_PDB, *TRAJECTORIES],
        params={
            "stride": stride,
            "skip_water": skip_water,
            "skip_convergence": skip_convergence,
        },
    )
    write_provenance_json(provenance_dir / "run_all_analyses.provenance.json", run_prov)
    
    # Load trajectory once
    print("Loading trajectory...")
    traj = load_and_concatenate_trajectories(stride=stride)
    save_trajectory_info(traj)
    
    # Core analyses (always run)
    analyses = [
        ("RMSD/RMSF", run_rmsd_analysis),
        (
            "IFP (MD-IFP)",
            lambda _t: run_ifp_analysis(
                topology_path=TOPOLOGY_PDB,
                trajectory_paths=TRAJECTORIES,
                output_dir=OUTPUT_DIR,
                generate_plots=True,
                save_results=True,
                verbose=True,
            ),
        ),
        ("Contacts + Q Native", run_contact_analysis),
        ("Binding Distances", run_binding_analysis),
        ("Clustering + PCA", run_clustering_analysis),
    ]
    
    # Optional water analysis
    if not skip_water:
        analyses.append(("Water + Residence Time (stride=10)", 
                        lambda t: run_water_analysis(stride=10, compute_residence=True)))
    
    # Optional convergence analysis
    if not skip_convergence:
        analyses.append(("Convergence Analysis", run_convergence_analysis))
    
    for name, func in analyses:
        print()
        try:
            if name.startswith("Water"):
                func(None)  # Water analysis loads its own trajectory with larger stride
            else:
                func(traj)

            analysis_prov = make_analysis_provenance(
                name=name,
                repo_root=PROJECT_ROOT,
                inputs=[TOPOLOGY_PDB, *TRAJECTORIES],
                outputs=[HTML_DIR, FIGURES_DIR, DATA_OUTPUT_DIR],
                params={
                    "stride": stride,
                    "skip_water": skip_water,
                    "skip_convergence": skip_convergence,
                },
                code_paths=[Path(__file__), Path(__file__).parent / "config.py"],
                status="ok",
            )
            safe_name = "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in name)
            write_provenance_json(provenance_dir / f"analysis.{safe_name}.provenance.json", analysis_prov)
        except Exception as e:
            print(f"ERROR in {name}: {e}")
            traceback.print_exc()

            analysis_prov = make_analysis_provenance(
                name=name,
                repo_root=PROJECT_ROOT,
                inputs=[TOPOLOGY_PDB, *TRAJECTORIES],
                outputs=[HTML_DIR, FIGURES_DIR, DATA_OUTPUT_DIR],
                params={
                    "stride": stride,
                    "skip_water": skip_water,
                    "skip_convergence": skip_convergence,
                },
                code_paths=[Path(__file__), Path(__file__).parent / "config.py"],
                status="error",
                error={
                    "type": type(e).__name__,
                    "message": str(e),
                    "traceback": traceback.format_exc(),
                },
            )
            safe_name = "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "_" for ch in name)
            write_provenance_json(provenance_dir / f"analysis.{safe_name}.provenance.json", analysis_prov)
    
    # Create dashboard
    print()
    print("=" * 70)
    print("CREATING DASHBOARD")
    print("=" * 70)
    create_dashboard_html()
    
    # Summary
    elapsed = time.time() - start_time
    print()
    print("=" * 70)
    print("[OK] ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"Total time: {elapsed/60:.1f} minutes")
    print(f"Output directory: {PROJECT_ROOT}")
    print()
    print("Generated files:")
    print(f"  HTML reports: {HTML_DIR}")
    print(f"  Figures: {FIGURES_DIR}")
    print(f"  Data: {DATA_OUTPUT_DIR}")
    print()
    print(f"Open dashboard: {HTML_DIR / '00_dashboard.html'}")
    print()
    print("Key outputs:")
    print("  - RMSD time series with equilibration detection")
    print("  - Q native contact fraction (differentiable)")
    print("  - Water residence time at interface")
    print("  - Convergence metrics (τ_auto, N_eff, SEM)")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="CCT-OSR1 Binding Analysis - Publication-Grade MD Pipeline"
    )
    parser.add_argument('--stride', type=int, default=1, 
                        help='Frame stride for loading (default: 1)')
    parser.add_argument('--skip-water', action='store_true', 
                        help='Skip water analysis (faster)')
    parser.add_argument('--skip-convergence', action='store_true',
                        help='Skip convergence analysis')
    parser.add_argument('--quick', action='store_true',
                        help='Quick mode: stride=10, skip water and convergence')
    
    args = parser.parse_args()
    
    if args.quick:
        run_all_analyses(stride=10, skip_water=True, skip_convergence=True)
    else:
        run_all_analyses(
            stride=args.stride, 
            skip_water=args.skip_water,
            skip_convergence=args.skip_convergence
        )

# 04 — MD execution & canonical analysis code

Scripts that actually produced the simulations and the trajectory metrics
(RMSD, RMSF, SASA) underlying `03_data` and the figures in `01_figures`.
Contents copied verbatim from the original project archive (paths inside the
docs refer to the original machine layout; see mapping below).

## md_protocol_amber/ — run the simulations (Paper Methods → "Molecular Dynamics")

- `runcomplex_paper_amber.py` — publication-grade OpenMM runner used for all
  CCT–βϕ complexes: ff14SB + TIP3P, 1.2 nm padding, 0.15 M NaCl, two-stage
  minimization → NVT heating 0–300 K (restrained) → NPT equilibration with
  gradual restraint release → 200 ns production × 2 replicas, HMR 4 amu
  (4 fs timestep), Langevin + MC barostat, PME 10 Å cutoff, CUDA/RTX 5080.
- `REPRODUCTION_COMMANDS.md` — exact per-system commands for the 6 complexes
  (OSR1 CCT–GRFQVT 2V3S, WNK2 CCT1 6FBK, SPAK CCT, WNK4 CCT2–PAKKKSCFQITSVTT,
  NRBP1 CCT–RWTC, FERRY3 CCT-like–βϕC).
- `requirements_MD.txt` — MD environment dependencies.

## smic_analysis/ — canonical trajectory analysis (Paper Methods → MD analysis)

The SMIC suite that generated the RMSD/RMSF and SASA numbers later replotted
into the paper figures:

- `run_analysis.py` — orchestrator (RMSD/RMSF, contacts, binding, water,
  clustering, convergence blocks).
- `analysis_rmsd.py` — RMSD + Cα RMSF from self-aligned trajectories
  (source of the corrected `rmsf.csv` and `rmsd.csv` in
  `03_data/organized_by_system_bundle_20260327/<system>/<R>/csv/`).
- `analysis_water.py` — per-residue SASA (Shrake–Rupley, MDTraj) + water
  interface metrics (source of `csv/sasa.csv` in the same bundle).
- `trajectory_loader.py`, `config.py`, `provenance.py`, `visualization.py` —
  support modules (topology/trajectory loading, per-system paths, run metadata,
  report figures).

Note: the RING/network modules of the suite (`analysis_contacts.py`,
`analysis_binding.py`, `analysis_ifp.py`, `analysis_interactions*.py`,
`analysis_network.py`, `ring_import.py`) are not included in this package;
they are available in the original SMIC project folder on demand.

## Not included (removed from package; available on demand)

- RING interaction-network protocol (`RING_PROTOCOL.md`): 50–150 ns window per
  replica, 2 windows × 3 ensembles of 100 conformations = 6 networks/replica,
  RING (Begue et al. 2025) geometric criteria, ≥20% frequency threshold,
  Cytoscape v3.10.4, interactions conserved in both replicas. This produced
  Fig. 1C / Fig. 5E interaction diagrams and Supplementary networks/tables.
- System definitions (`SYSTEMS_MANIFEST.csv`, example complex PDB) referenced
  by `REPRODUCTION_COMMANDS.md`.

## Original-path mapping

| In docs                        | Original location                                              |
|--------------------------------|----------------------------------------------------------------|
| `00_MD_protocol/…`             | this package `04_md_run_and_analysis/md_protocol_amber/`       |
| `01_systems/…`                 | original project archive (not in package; on demand)           |
| `02_analysis_scripts/smic/…`   | this package `04_md_run_and_analysis/smic_analysis/`           |
| `04_RING_validation/…`         | original project archive (not in package; on demand)           |
| `INVENTARIO/…` (trajectories)  | `d:\DESCARGAS_PESADAS\controles_novel_ccts\INVENTARIO\…`       |

Flow: `md_protocol_amber` → trajectories (DCD) → `smic_analysis` → CSVs in
`03_data` → `02figure_scripts/_build_md_per_system_panels.py` → figures in `01_figures`.

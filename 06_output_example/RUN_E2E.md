# RUN_E2E - reproducible runbook (local commands only)

Order: input model -> preparation + MD -> trajectory analysis (CSV) -> final frame -> plots -> verification.
Every command and parameter below comes from the run logs (`simulation.log`) and the scripts in this
bundle; anything that could not be recovered is marked **NOT RECOVERED**.

## 0. Environments

```bash
# MD (Linux + NVIDIA GPU)
conda env create -f 06_output_example/environment/environment_md_run.yml
conda activate cct_md_openmm82
# Post-processing (any OS, CPU)
python -m venv .venv && .venv/bin/pip install -r 06_output_example/environment/requirements_postprocessing.txt
```

## 1. Inputs

| System | Input model (in `05_run_artifacts/<system>/inputs/`) | Source |
|---|---|---|
| osr1_cct | `osr1_cct_input_2V3S_rcsb.pdb` | PDB 2V3S, unmodified RCSB file (`curl -O https://files.rcsb.org/download/2V3S.pdb`) |
| wnk2_cctl1 | `wnk2_cctl1_input_6FBK_rcsb.pdb` | PDB 6FBK, unmodified RCSB file (`curl -O https://files.rcsb.org/download/6FBK.pdb`) |
| spak_cct | `spak_cct_input_stk39-wnk4-relax.pdb` | AF3 model + Rosetta relax (STK39 CCT + WNK4 peptide) |
| wnk4_cctl2 | `wnk4_cctl2_input_WNK4CCT2_RpA.pdb` | AF3 model + Rosetta relax (WNK4 CCT2 + TSC22D2 peptide) |
| nrbp1_cct | **NOT RECOVERED** (`NRBP1PF2-RpC relax.pdb` per log) | only the prepared / equilibrated solvated systems survive |
| ferry_cctl | `ferry_cctl_input_c12orf4_rfx.pdb` | FERRY3 (= C12orf4) CCT-like + peptide model |

## 2. Preparation + MD (one command per replica)

Runner: `06_output_example/scripts/01_run_md__runcomplex_paper_amber.py`
(`04_md_run_and_analysis/md_protocol_amber/runcomplex_paper_amber.py` is the canonical implementation with
`amber14/protein.ff14SB.xml + amber14/tip3p.xml`).

Common protocol (from the runner): PDBFixer (missing heavy atoms, remove waters/heterogens)
when installed, hydrogens at pH 7.4, cubic box 1.2 nm padding, 0.15 M NaCl; HBond constraints, HMR 4 amu,
PME 1.0 nm, CMMotionRemover every 100 steps; minimisation tol 100 then 10 kJ/mol/nm; NVT 0->300 K in
500 ps (2 fs, backbone CA/C/N restraints 500 kJ/mol/nm^2); NPT 1 bar 2 ns with gradual restraint release
(MonteCarloBarostat every 25 steps); production LangevinMiddle 300 K, 1/ps, 4 fs; DCD + energy CSV every
10 ps (2,500 steps); checkpoint + State XML every 10 ns; `<run>_final.xml` at the end.

```bash
python -u runcomplex_paper_amber.py --pdb <INPUT.pdb> --output_dir runs/replica_<N> --replica <N> \
       --ns 200 --seed <SEED> --platform CUDA --gpu_id 0 --prepare  > logs/replica_<N>.log 2>&1
```

| System | Force field | Prep | Seeds R1 / R2 | Atoms R1 / R2 | Notes |
|---|---|---|---|---|---|
| osr1_cct | Amber14SB + TIP3P | PDBFixer | 456 / 789 | 56,181 / 54,428 | initial Nov-2025 run (state at 22.5 ns, energy CSV to ~32 ns); R1 production restarted from `R1_equilibrated_solvated.pdb` (seed random) to 202.5 ns; R2 continuation to 176.44 ns |
| wnk2_cctl1 | Amber14SB + TIP3P | PDBFixer | 42 / 123 | 23,779 / 23,746 | initial run (energy CSV to ~60 ns), then continuation from the 52.5-ns state (target 300 ns, R1 to 184.18 ns, R2 to 182.68 ns) |
| spak_cct | Amber14SB + TIP3P | PDBFixer | 1001 / 1002 | 27,549 / 27,549 | complete, 202.5 ns; see `06_output_example` |
| nrbp1_cct | Amber14SB + TIP3P | basic OpenMM (no PDBFixer) | random / random | 23,189 / 23,189 | R1 = replica1/NRBP1_rep1, R2 = replica2/results_nrbp1pf2 |
| wnk4_cctl2 | Amber14SB + TIP3P | basic OpenMM (no PDBFixer) | random / random | 47,805 / 47,805 | complete, 202.5 ns |
| ferry_cctl | Amber14SB + TIP3P | - | - | 47,143 / 47,167 | complete, 202.5 ns |

## 3. Trajectory analysis -> CSV (needs the DCDs, not shipped)

```bash
python 06_output_example/scripts/02_analysis_from_dcd.py --top <R>_equilibrated_solvated.pdb --traj replica_<N>.dcd --stride 5 --out <R>/analysis_rerun
```
Calls `compute_all_rmsd` / `compute_all_rmsf` from `04_md_run_and_analysis/md_analysis/analysis_rmsd.py`
(self-aligned, frame 0 reference) and MDTraj Shrake-Rupley per-residue SASA (as `analysis_water.py`),
stride 5 = 50 ps as in the original analysis logs. The curated CSVs (`csv/R*_rmsd|rmsf|sasa.csv`) were
produced by the SMIC suite with this code; UniProt renumbering for the figures is applied by
`02figure_scripts` from `03_data/uniprot_mapping.json`. DSSP generation script: **NOT RECOVERED** (output only).

## 4. Final frame

```bash
python 06_output_example/scripts/03_extract_final_frame.py --state <run>_final.xml --topology <R>_equilibrated_solvated.pdb --out <R>_final_frame_202.5ns_solute.pdb
```
Available for osr1_cct R1, spak_cct R1/R2, wnk4_cctl2 R1/R2 (states at 202.5 ns). Others: no production-end
State/checkpoint locally (see 05 README).

## 5. Plots

```bash
python 06_output_example/scripts/05_plot_outputs.py --replica_dir 06_output_example/R1 --label "SPAK CCT R1"   # per replica
py -3.13 02figure_scripts/_build_md_per_system_panels.py                                                   # paper panels
```

## 6. Verification

```powershell
powershell -ExecutionPolicy Bypass -File 06_output_example/run_e2e.ps1 -Python <python>
```
Writes `VERIFICATION_R1/R2.json` and `SHA256SUMS.txt` (results in `06_output_example/README.md`).


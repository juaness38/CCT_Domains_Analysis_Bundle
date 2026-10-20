# 06 - Output example: SPAK CCT, one complete end-to-end run

One full, real run of the paper's MD pipeline for **SPAK CCT (STK39, UniProt Q9UEW8) + WNK4
peptide KPQLVGRFQVTSSKE**, 2 replicas x 200 ns (March 2026), with every output that exists
locally. All files are original run outputs except those marked *regenerated* (made on
2026-10-19 by the scripts in `scripts/` and checked against the originals).

Paper: Magana-Avila et al., *Determinants of CCT-motif specificity in WNK signaling and
expansion of CCT-like domains*, bioRxiv 2026, https://doi.org/10.64898/2026.04.15.718729

## Pipeline

| Step | What | Script | Output here |
|---|---|---|---|
| 0 | Input model (Rosetta-relaxed SPAK CCT + WNK4 peptide complex) | - | `inputs/spak_cct_input_stk39-wnk4-relax.pdb` |
| 1 | Prep (PDBFixer, pH 7.4, 1.2 nm TIP3P box, 0.15 M NaCl) -> 2-stage minimisation -> NVT 0-300 K 500 ps (backbone restraints 500 kJ/mol/nm^2) -> NPT 1 bar 2 ns (restraint release) -> 200 ns production (Amber14SB / TIP3P, HMR 4 amu, 4 fs, LangevinMiddle 300 K 1/ps, PME 1.0 nm, MC barostat) | `scripts/01_run_md__runcomplex_paper_amber.py` | `R*/md/` |
| 2 | Trajectory analysis: Calpha/backbone RMSD, Calpha RMSF (self-aligned), per-residue SASA (Shrake-Rupley), stride 5 (50 ps) | `scripts/02_analysis_from_dcd.py` (SMIC functions from `04_md_run_and_analysis/md_analysis`) | `R*/analysis_csv/` (original) |
| 3 | Final frame (t = 202.5 ns) from the final OpenMM State, solute only | `scripts/03_extract_final_frame.py` | `R*/final_frame/` (regenerated) |
| 4 | Smoke test: the step-2 code runs on a 2-frame (equilibrated, final) trajectory and gives the same final-frame RMSD as step 6 | `scripts/02_analysis_from_dcd.py --state` | `R*/analysis_smoketest/` (regenerated) |
| 5 | Plots | `scripts/05_plot_outputs.py`; paper panels: `02figure_scripts/_build_md_per_system_panels.py` | `R*/figures/` |
| 6 | Verification + SHA256 | `scripts/06_verify_example.py`, `run_e2e.ps1` | `VERIFICATION_R*.json`, `SHA256SUMS.txt` |

Run everything that is reproducible on a normal PC (steps 3-6):

```powershell
powershell -ExecutionPolicy Bypass -File run_e2e.ps1 -Python <path\to\python.exe>
# add -RunMD to repeat step 1 on a CUDA GPU, and -Dcd1/-Dcd2 <file.dcd> to repeat step 2
```

Exact commands, parameters and the other five systems: `RUN_E2E.md`.
Environments: `environment/environment_md_run.yml` (MD) and `environment/requirements_postprocessing.txt` (steps 2-6).

## Files

```
inputs/spak_cct_input_stk39-wnk4-relax.pdb   original input model (identical to the file used in the run)
scripts/                                    01 MD runner, 02 analysis, 03 final frame, 05 plots, 06 verify
R1/, R2/
  md/prepared_solvated.pdb                  solvated system after preparation (R1 only; R2 copy not recovered)
  md/equilibrated_solvated.pdb              system after NVT+NPT = topology of the DCD / start of production
  md/equilibrated_state.xml                 OpenMM State after equilibration (R2 only; R1 copy not recovered)
  md/final_state.xml                        OpenMM State at the end of production (step 51,250,000; 202.5 ns)
  md/last_checkpoint.chk                    last OpenMM checkpoint (binary, platform-specific)
  md/energy.csv                             OpenMM StateDataReporter: step, time, energies, T, volume, density, speed
  analysis_csv/rmsd.csv rmsf.csv sasa.csv   original SMIC analysis outputs (PDB numbering, residue 1 = UniProt 448)
  analysis_csv/dssp_per_residue.csv         secondary structure track used in the paper panels
  final_frame/spak_cct_R*_final_frame_202.5ns_solute.pdb   regenerated, 1,741 atoms (protein + peptide)
  figures/original_*.png                    original per-replica plots
  figures/regenerated_*.png                 regenerated from the CSVs by step 5
  analysis_smoketest/                       step-4 output
VERIFICATION_R1.json, VERIFICATION_R2.json
SHA256SUMS.txt
```

## Verification (2026-10-19)

| Check | R1 | R2 |
|---|---|---|
| analysis CSVs byte-identical to the recovered originals (`05_run_artifacts/spak_cct/csv`) | yes (4/4) | yes (4/4) |
| `energy.csv` ends at the time stored in `final_state.xml` (202.5 ns) | yes | yes |
| final frame CA RMSD vs equilibrated structure / last value of `rmsd.csv` | 2.36 / 2.52 A | 1.82 / 1.90 A |
| per-residue SASA of the final frame vs trajectory-mean SASA (Pearson r, 113 residues) | 0.92 | 0.92 |

The RMSD difference is expected: `rmsd.csv` stops at the last analysed frame (stride 5, 199.95 ns)
while the final frame is at 202.5 ns. A full byte-level re-run of step 2 needs the DCD.

## Not included / caveats

- **Trajectories** `replica_1.dcd`, `replica_2.dcd` (20,000 frames x 27,549 atoms, about 6.2 GiB each)
  are not in this folder or repository; they remain in the authors' storage.
- `rmsd.csv` `time_ps`/`time_ns` columns are x5 mis-scaled (4,000 rows end at "999.75 ns"); the true
  spacing is 50 ps (0-200 ns). The paper figure script corrects this through `03_data/replot_manifest.csv`.
- The MD runner in `scripts/` is the closest surviving version: its force field, preparation and log
  messages match `simulation.log` line for line, but the exact file used in March 2026 was not recovered.
- The DSSP generation script is not in the bundle (only its output).



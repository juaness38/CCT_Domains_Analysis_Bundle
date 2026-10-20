# 05 - Run artifacts (6 systems x 2 replicas)

Per system exactly three folders:

```
<system>/
  inputs/       original input model + the solvated systems written by the MD runner
                (<sys>_input_*.pdb, R*_prepared_solvated.pdb, R*_equilibrated_solvated.pdb)
  csv/          real outputs only: R*_rmsd.csv, R*_rmsf.csv, R*_sasa.csv (SMIC analysis, PDB numbering),
                R*_dssp_per_residue.csv (UniProt numbering, from 03_data), R*_energy*.csv (OpenMM StateDataReporter)
  final_frame/  R*_final_frame_202.5ns_solute.pdb (protein + peptide at the end of production)
```

Paths for the MD run, analysis and figures: `../06_output_example/RUN_E2E.md`. The complete SPAK run
(logs, states, checkpoints, plots, verification) is in `../06_output_example/`.

## Systems

| Folder | Paper system | Input model | Inputs | CSV | Final frame |
|---|---|---|---|---|---|
| `osr1_cct` | OSR1 CCT + GRFQVT (WNK4) | PDB 2V3S (RCSB, unmodified) | yes | R1, R2 | R1 (202.50 ns), R2 (176.44 ns)* |
| `spak_cct` | SPAK CCT + KPQLVGRFQVTSSKE (WNK4) | `stk39-wnk4-relax.pdb` | yes | R1, R2 | R1 (202.50 ns), R2 (202.50 ns) |
| `wnk2_cctl1` | WNK2 CCT-like 1 + WNK1 RFxV peptide | PDB 6FBK (RCSB, unmodified) | yes | R1, R2 | R1 (184.18 ns)*, R2 (182.68 ns)* |
| `nrbp1_cct` | NRBP1 CCT + EPYRRGRWTCMEYYE (TSC22D2) | `NRBP1PF2-RpC relax.pdb` (AF3 + Rosetta) | yes | R1, R2 | R1 (202.50 ns), R2 (202.50 ns) |
| `wnk4_cctl2` | WNK4 CCT-like 2 + PAKKKSCFQITSVTT (TSC22D2) | `WNK4CCT2_RpA.pdb` | yes | R1, R2 | R1 (202.50 ns), R2 (202.50 ns) |
| `ferry_cctl` | FERRY3 (= C12orf4) CCT-like + TSC22D2 peptide | `c12orf4_rfx.pdb` | yes | R1, R2 | R1 (202.50 ns), R2 (202.50 ns) |

\* Shorter production (< 200 ns) due to continuation runs: OSR1 R2 ended at 176.44 ns; WNK2 R1 ended at 184.18 ns and R2 ended at 182.68 ns.

## Notes per system

- **osr1_cct** - `csv/R*_energy_production.csv` is the analysed production run (to 202.5 ns for R1);
  `R*_energy_initial_segment.csv` is the first Nov-2025 run (~32 ns) - do not concatenate.
  R2 final frame was extracted from the production continuation DCD ending at 176.44 ns (`R2_final_frame_176.44ns_solute.pdb`).
  `R2_rmsd.csv` here is the recovered original (120,366 B); `03_data/.../osr1_cct/R2/csv/rmsd.csv` (97,743 B,
  used by the figures) is a corrected version and `rmsd_CORRUPTED_BACKUP.csv` there is a third one - they differ.
- **wnk2_cctl1** - energy in two parts (initial run + continuation from 52.5 ns). Final frames were extracted
  from the production continuation DCDs (`R1_final_frame_184.18ns_solute.pdb` at 184.18 ns and `R2_final_frame_182.68ns_solute.pdb` at 182.68 ns).
- **nrbp1_cct** - original input model `NRBP1PF2-RpC relax.pdb` added. Trajectory analysis and DCD coordinates
  confirm the convention of paper figures and 03_data (R1 = `NRBP1_replica1/NRBP1_rep1`, R2 = `NRBP1_replica2/results_nrbp1pf2`).
  Final frames extracted from DCD at 202.50 ns.
- **ferry_cctl** - FERRY3 is the current name of C12orf4; original run folder name `C12OSR4_R1/R2`,
  energy files originally `replica_1_energy_C12OR3.csv` / `replica_2_energy_C12OR4.csv`. Final frames extracted
  from DCD at 202.50 ns.
- **wnk4_cctl2** - run with the `protein.ff14SB + TIP3P` runner variant (see RUN_E2E.md).
- Not included anywhere in the repository: trajectories (`*.dcd`, 6-45 GB each), solvated multi-frame
  dumps, MMGBSA, interactive HTML. SHA256 of every file: `SHA256SUMS.txt`.


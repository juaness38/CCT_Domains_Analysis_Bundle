# Computational Methods & Software Inventory

This document details the exact software tools, versions, simulation parameters, repository locations, official documentation links, and formal citations (DOIs) used in the molecular dynamics (MD) simulations and post-processing of the six CCT / CCT-like domain systems reported in:

> Magaña-Ávila et al., *Determinants of CCT-motif specificity in WNK signaling and expansion of CCT-like domains*, bioRxiv 2026, [DOI: 10.64898/2026.04.15.718729](https://doi.org/10.64898/2026.04.15.718729).

Only tools and libraries that appear in the simulation scripts, execution logs, OpenMM state/checkpoint XMLs, or post-processing environments are included. **FreeSASA was NOT used**; per-residue solvent accessible surface area (SASA) was computed using the Shrake-Rupley algorithm as implemented in MDTraj.

---

## 1. Software Inventory & Repository Locations

| Software / Library | Version | Role in Study | Primary Location in Repository | Documentation | Citation / DOI |
|---|---|---|---|---|---|
| **OpenMM** | `8.2.0` (simulations)<br>`8.3.1` / `8.5.1` (analysis) | MD simulation engine (CUDA platform) & State XML deserializer | `06_output_example/scripts/01_run_md__runcomplex_paper_amber.py`<br>`04_md_run_and_analysis/md_protocol_amber/runcomplex_paper_amber.py`<br>`06_output_example/scripts/03_extract_final_frame.py` | [openmm.org](https://openmm.org) | [10.1021/acs.jpcb.3c06639](https://doi.org/10.1021/acs.jpcb.3c06639)<br>[10.1371/journal.pcbi.1005630](https://doi.org/10.1371/journal.pcbi.1005630) |
| **PDBFixer** | `1.12.0` | System preparation: protonation, missing atoms, solvation, ion placement | `06_output_example/scripts/01_run_md__runcomplex_paper_amber.py` | [github.com/openmm/pdbfixer](https://github.com/openmm/pdbfixer) | Part of OpenMM ecosystem |
| **MDTraj** | `1.11.0` | Trajectory loading, PBC molecular imaging, solute stripping, Shrake-Rupley SASA | `06_output_example/scripts/02_analysis_from_dcd.py`<br>`06_output_example/scripts/03_extract_final_frame.py`<br>`04_md_run_and_analysis/md_analysis/analysis_water.py` | [mdtraj.org](https://mdtraj.org) | [10.1016/j.bpj.2015.08.015](https://doi.org/10.1016/j.bpj.2015.08.015) |
| **MDAnalysis** | `2.9.0` / `2.10.0` | Backbone and Cα alignment, superposition, RMSD and RMSF calculations | `04_md_run_and_analysis/md_analysis/analysis_rmsd.py`<br>`06_output_example/scripts/runtime_config.py` | [mdanalysis.org](https://www.mdanalysis.org) | [10.1002/jcc.21787](https://doi.org/10.1002/jcc.21787)<br>[10.25080/Majora-629e541a-00e](https://doi.org/10.25080/Majora-629e541a-00e) |
| **DSSP** | Algorithm (Kabsch & Sander) | Secondary structure assignment over production trajectory (`dssp_per_residue.csv`) | `03_data/organized_by_system_bundle_20260327/<sys>/<R>/dssp/`<br>`05_run_artifacts/<sys>/csv/R*_dssp_per_residue.csv` | [swift.cmbi.umcn.nl](https://swift.cmbi.umcn.nl/gv/dssp/) | [10.1002/bip.360221211](https://doi.org/10.1002/bip.360221211) |
| **AlphaFold 3** | Initial models | Complex structure prediction for STK39-WNK4, WNK4-TSC22D2, C12orf4-TSC22D2, NRBP1-TSC22D2 | `05_run_artifacts/<sys>/inputs/` | [alphafoldserver.com](https://alphafoldserver.com) | [10.1038/s41586-024-07487-w](https://doi.org/10.1038/s41586-024-07487-w) |
| **Rosetta** | FastRelax protocol | Pre-equilibration all-atom relaxation of AF3 complexes (`-relax.pdb`) | `05_run_artifacts/<sys>/inputs/` | [rosettacommons.org](https://www.rosettacommons.org) | [10.1002/pro.2389](https://doi.org/10.1002/pro.2389) |
| **NumPy** | `2.4.2` / `2.5.2` | Numerical computing, vector algebra, array manipulations | Core Python scripts | [numpy.org](https://numpy.org) | [10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2) |
| **Pandas** | `2.3.3` | Tabular data processing for RMSD, RMSF, SASA, energy logs, and UniProt mappings | All analysis and figure scripts | [pandas.pydata.org](https://pandas.pydata.org) | [10.5281/zenodo.3509134](https://doi.org/10.5281/zenodo.3509134) |
| **SciPy** | `1.17.1` | Statistical analysis, correlation calculations (Pearson *r*), curve fits | Post-processing & verification scripts | [scipy.org](https://scipy.org) | [10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2) |
| **Matplotlib** | `3.11.1` | Publication-ready 600-dpi panels and validation figures | `02figure_scripts/_build_md_per_system_panels.py`<br>`06_output_example/scripts/05_plot_outputs.py` | [matplotlib.org](https://matplotlib.org) | [10.1109/MCSE.2007.55](https://doi.org/10.1109/MCSE.2007.55) |
| **BioPython** | `1.88` | PDB handling and residue metadata verification | Environment dependencies | [biopython.org](https://biopython.org) | [10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163) |
| **RING** | `3.0` | Residue interaction network generation across structural ensembles | `07_networks/` pipeline | [ring.biocomputingup.it](https://ring.biocomputingup.it) | [10.1093/nar/gkaa338](https://doi.org/10.1093/nar/gkaa338) |
| **Cytoscape** | `3.10.4` | Network visualization and graph rendering | Downstream analysis of `07_networks/` JSONs | [cytoscape.org](https://cytoscape.org) | [10.1101/gr.1239303](https://doi.org/10.1101/gr.1239303) |
| **openpyxl** | `3.1.5` | Excel reference interaction tables reader | `07_networks/select_representative_state.py` | [openpyxl.readthedocs.io](https://openpyxl.readthedocs.io) | [10.5281/zenodo.592534](https://doi.org/10.5281/zenodo.592534) |

---

## 2. Experimental Structures & Model Inputs

1. **Experimental Crystallographic Structures (RCSB PDB)**:
   - **OSR1 CCT**: PDB ID **2V3S** (*Crystal structure of the C-terminal domain of the human kinase OSR1 in complex with a WNK4 peptide*, Villa et al., J. Biol. Chem. 2007, [DOI: 10.1074/jbc.M706698200](https://doi.org/10.1074/jbc.M706698200)). Used exactly as published, without structural modifications.
   - **WNK2 CCT-like 1**: PDB ID **6FBK** (*Crystal structure of the CCT1 domain of human WNK2 in complex with WNK1 peptide*, Pinkas et al., Structure 2018, [DOI: 10.1016/j.str.2018.06.012](https://doi.org/10.1016/j.str.2018.06.012)). Used exactly as published, without structural modifications.

2. **Computational Complex Predictions**:
   - **SPAK CCT**: Generated via AlphaFold 3 (STK39 construct residues 448–545 with WNK4 peptide `KPQLVGRFQVTSSKE`, residues 1010–1024) followed by Rosetta FastRelax (`stk39-wnk4-relax.pdb`, dated 2026-03-13).
   - **WNK4 CCT-like 2**: Generated via AlphaFold 3 (WNK4 CCT2 construct residues 667–751 with TSC22D2 peptide `PAKKKSCFQITSVTT`, residues 5–19) followed by Rosetta FastRelax (`WNK4CCT2_RpA.pdb`, dated 2025-11-29).
   - **FERRY3 (C12orf4) CCT-like**: Generated via AlphaFold 3 with TSC22D2 peptide (`c12orf4_rfx.pdb`).
   - **NRBP1 CCT**: Generated via AlphaFold 3 with TSC22D2 peptide (`EPYRRGRWTCMEYYE`, residues 173–187) followed by Rosetta relax (`NRBP1PF2-RpC relax.pdb`, SHA256 `4fc6f4e3d308d3f6044233c8d2ee2715377687e3520487ac60665e1274a1b200`).

---

## 3. Molecular Dynamics Simulation Protocol (OpenMM)

Simulations were conducted using the GPU-accelerated **OpenMM 8.2.0** engine on NVIDIA hardware with the `CUDA` platform.

### 3.1. System Preparation & Solvation
1. **PDBFixer 1.12.0**:
   - Missing heavy atoms identified and added (`findMissingAtoms`, `addMissingAtoms`).
   - Crystallographic water molecules and heterogens removed (`removeHeterogens(keepWater=False)`).
   - Hydrogen atoms added according to neutral/physiological protonation states at pH 7.4 (`addMissingHydrogens(pH=7.4)`).
2. **Solvation & Ion Neutralization**:
   - Explicit solvent box: cubic box with 1.2 nm minimum padding from solute to box boundary (`addSolvent(padding=1.2*nanometer)`).
   - Ionic strength: 0.15 M sodium chloride (`0.15 * molar`), neutralizing overall net charge with appropriate excess of Na⁺ or Cl⁻ ions.

### 3.2. Force Field Parameters
- **Protein & Water/Ions**:
  - All systems: **Amber14SB (`amber14/protein.ff14SB.xml`)** force field with explicit **TIP3P (`amber14/tip3p.xml`)** water model (`04_md_run_and_analysis/md_protocol_amber/runcomplex_paper_amber.py`).
- **Hydrogen Mass Repartitioning (HMR)**:
  - Non-water hydrogen mass increased to **4.0 amu** (`hydrogenMass=4.0*amu`), redistributing mass from connected heavy atoms to permit long integration timesteps.
- **Constraints & Electrostatics**:
  - Constraints: `HBonds` (rigid bonds for all bonds involving hydrogen atoms).
  - Electrostatics: Particle Mesh Ewald (**PME**) with a real-space cutoff distance of **1.0 nm** (`nonbondedCutoff=1.0*nanometer`).
  - Center of mass motion removal: `CMMotionRemover` applied every 100 simulation steps.

### 3.3. Simulation Stages & Ensembles
1. **Phase 1: Energy Minimization**:
   - Two-stage steepest descent / L-BFGS minimization in OpenMM:
     - Coarse stage: tolerance $100\text{ kJ/mol/nm}$.
     - Fine stage: tolerance $10.0\text{ kJ/mol/nm}$ (maximum 5,000 steps).
2. **Phase 2: NVT Heating**:
   - Ensemble: Canonical ($NVT$), temperature ramp from $0\text{ K} \rightarrow 300\text{ K}$ over $500\text{ ps}$ (timestep $2\text{ fs}$).
   - Restraints: Harmonic positional restraints on protein and peptide backbone atoms ($\text{C}\alpha, \text{C}, \text{N}$) with force constant $k = 500\text{ kJ/mol/nm}^2$.
3. **Phase 3: NPT Equilibration**:
   - Ensemble: Isothermal-isobaric ($NPT$) at $300\text{ K}$ and $1\text{ bar}$ for $2,000\text{ ps}$ ($2.0\text{ ns}$).
   - Barostat: `MonteCarloBarostat` (volume evaluation every 25 steps).
   - Thermostat: `LangevinMiddleIntegrator` with collision frequency $\gamma = 1.0\text{ ps}^{-1}$.
   - Stepwise restraint release over 5 equal blocks (400 ps each):
     - Block 1: 100% ($k = 500\text{ kJ/mol/nm}^2$)
     - Block 2: 50% ($k = 250\text{ kJ/mol/nm}^2$)
     - Block 3: 25% ($k = 125\text{ kJ/mol/nm}^2$)
     - Block 4: 10% ($k = 50\text{ kJ/mol/nm}^2$)
     - Block 5: 0% (unrestrained, $k = 0$)
   - Cumulative equilibration time: $2.5\text{ ns}$ ($0.5\text{ ns } NVT + 2.0\text{ ns } NPT$).
4. **Phase 4: Production MD**:
   - Ensemble: $NPT$ at $300\text{ K}$ and $1\text{ bar}$.
   - Timestep: **4.0 fs** enabled by HMR and HBond constraints.
   - Target production duration: **200.0 ns** (50,000,000 integration steps, bringing simulation clock to $t = 202.5\text{ ns}$).
   - Trajectory recording: DCD reporter every 2,500 steps (**10 ps** intervals).
   - Energy logging: potential, kinetic, total energy and temperature logged every 2,500 steps (**10 ps** intervals).
   - State checkpoints: OpenMM checkpoint (`.chk`) and serialized State XML (`.xml`) every $10\text{ ns}$ and at final completion.

### 3.4. Production Run Durations & Continuation Segments
The OpenMM clock records time including the 2.5 ns equilibration phase ($t = 2.51\text{ ns}$ to $t = 202.50\text{ ns}$ for standard 200 ns production):
- **Full 200 ns Production ($t = 202.50\text{ ns}$, 20,000 DCD frames)**:
  - `spak_cct` (R1, R2)
  - `wnk4_cctl2` (R1, R2)
  - `nrbp1_cct` (R1, R2)
  - `ferry_cctl` (R1, R2)
  - `osr1_cct` (R1)
- **Continuations Shorter than 200 ns (< 200 ns)**:
  - **`osr1_cct` R2**: Total duration **176.44 ns** ($15,394$ production DCD frames). Initial run performed Nov 2025 (State saved at 22.5 ns); production continuation ran from $22.51\text{ ns} \rightarrow 176.44\text{ ns}$.
  - **`wnk2_cctl1` R1**: Total duration **184.18 ns** ($13,168$ production DCD frames). Production continuation restarted from the $52.5\text{ ns}$ state ($52.51\text{ ns} \rightarrow 184.18\text{ ns}$).
  - **`wnk2_cctl1` R2**: Total duration **182.68 ns** ($13,018$ production DCD frames). Production continuation restarted from the $52.5\text{ ns}$ state ($52.51\text{ ns} \rightarrow 182.68\text{ ns}$).

---

## 4. Trajectory Post-Processing & Analysis

All post-processing was executed using Python 3.11 / 3.13 with MDTraj, MDAnalysis, NumPy, and Pandas.

### 4.1. PBC Imaging & Solute Extraction
Prior to analysis, trajectories and final frames were processed to reconstruct intact molecular complexes across periodic boundaries:
1. `MDTraj.Trajectory.image_molecules`:
   - Anchor molecule designated as the protein chain (`residue.is_protein` largest component).
   - Peptide and adjacent molecules mapped into the primary periodic unit cell nearest to the anchor.
2. Solvent and counterions stripped:
   - Selection: `not water and not resname NA and not resname CL` (`atom_slice`).

### 4.2. RMSD & RMSF Calculation
Calculated via `04_md_run_and_analysis/md_analysis/analysis_rmsd.py` (SMIC suite) using **MDAnalysis**:
- **Trajectory alignment**: Coordinate frames aligned to the reference structure (production frame 0, corresponding to $t = 2.51\text{ ns}$) based on backbone atoms ($\text{N}, \text{C}\alpha, \text{C}, \text{O}$).
- **RMSD**: Computed along the trajectory sampled every 5 DCD frames (**stride 5 = 50 ps**, 4,000 points for full 200 ns runs):
  - Kinase/domain Cα RMSD: `protein_ca` selection.
  - Peptide Cα RMSD: `peptide_ca` selection.
  - Complex Cα RMSD: `complex_ca` selection.
  - Pocket-specific Cα RMSD: binding site residues (`pocket_ca` selection).
- **RMSF**: Root Mean Square Fluctuation computed per Cα atom over the production trajectory after aligning frames to the average trajectory structure.

### 4.3. Solvent Accessible Surface Area (SASA)
- **Algorithm**: **Shrake-Rupley numerical method** (MDTraj `mdtraj.shrake_rupley(traj, mode="residue", probe_radius=0.14)`).
- **Parameters**: Probe sphere radius $r_{\text{probe}} = 1.4\text{ \AA}$ ($0.14\text{ nm}$), 960 test points per atom sphere.
- **Output metrics**:
  - Residue-level mean SASA ($\text{nm}^2$) averaged across all sampled trajectory frames (stride 5 = 50 ps).
  - Residue-level standard deviation ($\text{nm}^2$) quantifying pocket exposure breathing.
- **Clarification**: **FreeSASA was NOT utilized** in any stage of this study.

### 4.4. Secondary Structure (DSSP)
- Secondary structure per residue across the trajectory was determined using the **Dictionary of Secondary Structure of Proteins (DSSP)** algorithm (Kabsch & Sander, 1983).
- Outputs are preserved in `03_data/.../dssp/dssp_per_residue.csv` and `05_run_artifacts/.../csv/R*_dssp_per_residue.csv`. The generation wrapper script was not recovered; output data files are preserved verbatim with verified canonical UniProt mapping.

### 4.5. Canonical UniProt Numbering
To ensure consistency across clinical mutations, homologs, and structural constructs, all residue positions in paper figures (`01_figures/panels/`) are converted from raw PDB numbering to canonical UniProt sequence numbering via `03_data/uniprot_mapping.json`:
- **OSR1** (UniProt [O95747](https://www.uniprot.org/uniprotkb/O95747)): Construct residues $433–527$; offset **+431** (PDB residue 1 corresponds to Gly tag, canonical residue 432).
- **SPAK** (UniProt [Q9UEW8](https://www.uniprot.org/uniprotkb/Q9UEW8)): Construct residues $448–545$; offset **+447**.
- **NRBP1** (UniProt [Q9UHY1](https://www.uniprot.org/uniprotkb/Q9UHY1)): Construct residues $441–535$; offset **+440**.
- **WNK2** (UniProt [Q9Y3S1](https://www.uniprot.org/uniprotkb/Q9Y3S1)): Block-mapped due to internal deletion in construct:
  - PDB residues 1–28 $\rightarrow$ canonical residues 455–482.
  - Internal deletion in construct: residues 483–489 (`KLKGKPK`).
  - PDB residues 29–83 $\rightarrow$ canonical residues 490–544.
- **WNK4 CCT-like 2** (UniProt [Q96J92](https://www.uniprot.org/uniprotkb/Q96J92)): Construct residues $667–751$; offset **+666**.
- **FERRY3 / C12orf4** (UniProt [Q9NQ89](https://www.uniprot.org/uniprotkb/Q9NQ89)): Direct mapping (offset **+0**).

---

## 5. Residue Interaction Networks & Ensemble Analysis (`07_networks/`)

Residue interaction networks were constructed across molecular dynamics trajectories using custom Python workflows and the RING 3 engine:

1. **Trajectory-Frame Sampling (`sample_ensembles.py`)**:
   - Production trajectory frames are sampled and split into multiple ensembles using an interleaved sampling scheme rather than block-wise partitioning. This ensures that every ensemble independently samples the entire temporal progression of the trajectory.
2. **Residue Interaction Network Calculation**:
   - Each ensemble is analyzed with the RING 3 algorithm to identify non-covalent interactions (hydrogen bonds, salt bridges, π–π stacking, cation–π, and van der Waals contacts) formatted as Cytoscape JSON graphs.
3. **Consensus Network Construction (`build_consensus.py`)**:
   - Individual ensemble networks are merged into a single consensus graph. Each interaction edge is assigned a reproducibility score ($f_{\text{consensus}}$) reflecting the proportion of independent ensembles in which the contact occurs, along with summary statistics (mean, min, max frequencies and distances).
4. **Interface Contact Filtering & Canonical Renumbering (`renumber_network.py`)**:
   - Consensus networks are filtered to retain inter-chain contacts (between the CCT domain and the peptide ligand) meeting specified frequency thresholds, while excluding intra-chain background interactions.
   - Residue numbering of the CCT domain is mapped against full-length FASTA sequences to match canonical UniProt numbering, ensuring direct alignment with biochemical literature.
5. **Representative Ensemble State Selection (`select_representative_state.py`)**:
   - The single frame whose interface contact distances minimize the sum of squared deviations from the consensus ensemble geometry is objectively identified and extracted as the representative structural model for visualization and structural figures.

---

## 6. Formal References (DOIs)

1. **OpenMM**:
   - Eastman, P., Galvelis, R., Peláez, R. P., et al. (2024). OpenMM 8: Molecular Dynamics Simulation with Machine Learning Potentials. *The Journal of Physical Chemistry B*, 128(1), 109–116. [DOI: 10.1021/acs.jpcb.3c06639](https://doi.org/10.1021/acs.jpcb.3c06639)
   - Eastman, P., Swails, J., Chodera, J. D., et al. (2017). OpenMM 7: Rapid development of high performance algorithms for molecular dynamics. *PLOS Computational Biology*, 13(7), e1005630. [DOI: 10.1371/journal.pcbi.1005630](https://doi.org/10.1371/journal.pcbi.1005630)
2. **MDTraj**:
   - McGibbon, R. T., Beauchamp, K. A., Harrigan, M. P., et al. (2015). MDTraj: A Modern Open Library for the Analysis of Molecular Dynamics Trajectories. *Biophysical Journal*, 109(8), 1528–1532. [DOI: 10.1016/j.bpj.2015.08.015](https://doi.org/10.1016/j.bpj.2015.08.015)
3. **MDAnalysis**:
   - Michaud-Agrawal, N., Denning, E. J., Woolf, T. B., & Beckstein, O. (2011). MDAnalysis: A toolkit for the analysis of molecular dynamics simulations. *Journal of Computational Chemistry*, 32(10), 2319–2327. [DOI: 10.1002/jcc.21787](https://doi.org/10.1002/jcc.21787)
   - Gowers, R. J., Linke, M., Barnoud, J., et al. (2016). MDAnalysis: A Python package for the rapid analysis of molecular dynamics simulations. *Proceedings of the 15th Python in Science Conference*, 98–105. [DOI: 10.25080/Majora-629e541a-00e](https://doi.org/10.25080/Majora-629e541a-00e)
4. **Shrake-Rupley SASA**:
   - Shrake, A., & Rupley, J. A. (1973). Environment and Exposure to Solvent of Protein Atoms. Lysozyme and Insulin. *Journal of Molecular Biology*, 79(2), 351–371. [DOI: 10.1016/0022-2836(73)90011-9](https://doi.org/10.1016/0022-2836(73)90011-9)
5. **DSSP Algorithm**:
   - Kabsch, W., & Sander, C. (1983). Dictionary of protein secondary structure: Pattern recognition of hydrogen-bonded and geometrical features. *Biopolymers*, 22(12), 2577–2637. [DOI: 10.1002/bip.360221211](https://doi.org/10.1002/bip.360221211)
6. **AlphaFold 3**:
   - Abramson, J., Adler, J., Dunger, J., et al. (2024). Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature*, 630, 493–500. [DOI: 10.1038/s41586-024-07487-w](https://doi.org/10.1038/s41586-024-07487-w)
7. **Rosetta FastRelax**:
   - Conway, P., O'Meara, M. J., & Kuhlman, B. (2014). Relaxation of backbone bond geometry improves protein energy landscape modeling. *Protein Science*, 23(1), 47–55. [DOI: 10.1002/pro.2389](https://doi.org/10.1002/pro.2389)
   - Nivón, L. G., Moretti, R., & Baker, D. (2013). A Smooth Small-Molecule Force Field Generates Clustered High-Quality Loop Predictions. *PLoS ONE*, 8(6), e67051. [DOI: 10.1371/journal.pone.0067051](https://doi.org/10.1371/journal.pone.0067051)
8. **Crystal Structures (RCSB PDB)**:
   - Villa, F., Goebel, J., Rafiqi, F. H., et al. (2007). Structural insights into the recognition of substrates and activators by the OSR1 kinase. *Journal of Biological Chemistry*, 282(44), 32549–32560. [DOI: 10.1074/jbc.M706698200](https://doi.org/10.1074/jbc.M706698200) [PDB: 2V3S]
   - Pinkas, D. M., Sanvitale, C. E., Bufton, J. C., et al. (2018). Structural and Biochemical Characterization of WNK Kinase C-Terminal Domains. *Structure*, 26(8), 1089–1098.e4. [DOI: 10.1016/j.str.2018.06.012](https://doi.org/10.1016/j.str.2018.06.012) [PDB: 6FBK]
9. **Scientific Python Stack**:
   - Harris, C. R., Millman, K. J., van der Walt, S. J., et al. (2020). Array programming with NumPy. *Nature*, 585, 357–362. [DOI: 10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2)
   - Virtanen, P., Gommers, R., Oliphant, T. E., et al. (2020). SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. *Nature Methods*, 17, 261–272. [DOI: 10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2)
   - Hunter, J. D. (2007). Matplotlib: A 2D Graphics Environment. *Computing in Science & Engineering*, 9(3), 90–95. [DOI: 10.1109/MCSE.2007.55](https://doi.org/10.1109/MCSE.2007.55)
   - Cock, P. J. A., Antao, T., Chang, J. T., et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*, 25(11), 1422–1423. [DOI: 10.1093/bioinformatics/btp163](https://doi.org/10.1093/bioinformatics/btp163)
10. **Residue Interaction Networks & Graph Analysis**:
    - Clementel, D., Del Conte, A., Monzon, A. M., et al. (2022). RING 3.0: fast calculation of non-covalent interactions in molecular structures. *Nucleic Acids Research*, 50(W1), W651–W656. [DOI: 10.1093/nar/gkaa338](https://doi.org/10.1093/nar/gkaa338)
    - Shannon, P., Markiel, A., Ozier, O., et al. (2003). Cytoscape: a software environment for integrated models of biomolecular interaction networks. *Genome Research*, 13(11), 2498–2504. [DOI: 10.1101/gr.1239303](https://doi.org/10.1101/gr.1239303)


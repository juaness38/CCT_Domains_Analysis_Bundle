# Reproducibility Package for: Determinants of CCT–motif specificity in WNK signaling and expansion ofCCT-like domains* — Magaña-Ávila et al., bioRxiv 2026
https://doi.org/10.64898/2026.04.15.718729

Canonical package for the Control-MD CCT/CCTL paper figures. It contains **exactly**
the scripts, data and final figures that produced the paper panels
(palette: protein/CCT green + peptide orange, UniProt residue numbering). Nothing
here is illustrative; the panels are regenerable from `02figure_scripts` + `03_data`.


## Structure

```
nature_validation_CCT/
├── 01_figures/
│   └── panels/                      6 systems × (PNG 600dpi + SVG) + all-systems overlay + DOCX
├── 02figure_scripts/                ← the exact figure code
│   ├── _build_md_per_system_panels.py      (main figure generator — THIS made the panels)
│   ├── _replot_md_controls_rmsd_rmsf.py    (base module: replica specs, manifest, time axis)
│   ├── _replot_md_controls_rmsf_only_v2.py (historical: produced the corrected RMSF extracts)
│   └── _compute_uniprot_offsets.py         (UniProt mapping via pairwise alignment)
└── 03_data/                         ← canonical input data (single source: organized_by_system)
    ├── replot_manifest.csv                 RMSD/RMSF source-of-truth manifest (authoritative
    │                                       runtime axes; RMSD CSVs resolved from the manifest)
    ├── uniprot_mapping.json                canonical UniProt mapping (protein blocks + ligand
    │                                       mappings, with evidence) — supersedes the original
    │                                       offset map (deleted; it was computed against
    │                                       non-canonical sequences and was off by one for
    │                                       osr1/wnk4 and wrong for wnk2)
    └── organized_by_system_bundle_20260327/
        <system>/<R1|R2>/csv/               rmsd.csv, rmsf.csv, sasa.csv
        <system>/<R1|R2>/dssp/              dssp_per_residue.csv (secondary-structure track)

└── 04_md_run_and_analysis/          ← MD execution + canonical analysis code
    ├── md_protocol_amber/                  OpenMM/ff14SB runner + reproduction commands (Methods: MD)
    └── md_analysis/                        SMIC suite: RMSD/RMSF + SASA(Shrake–Rupley) → CSVs in 03_data
                                            (RING/network modules & system manifests: on demand,
                                            original SMIC project folder)

├── 05_run_artifacts/                 6 systems: <system>/{inputs, csv, final_frame} (see 05_run_artifacts/README.md)
├── 06_output_example/                one complete SPAK CCT run end to end + RUN_E2E.md runbook,
│                                     run_e2e.ps1, environments, verification, SHA256SUMS
└── 07_networks/                      Residue interaction network pipeline (interleaved ensemble
                                      sampling, consensus building, interface contact filtering,
                                      canonical renumbering, representative state selection)
```

## Systems

| Label       | Slug       | UniProt            | Offset (PDB → UniProt) |
|-------------|------------|--------------------|------------------------|
| OSR1 CCT    | osr1_cct   | O95747 (OXSR1)     | +431                   |
| SPAK CCT    | spak_cct   | Q9UEW8 (STK39)     | +447                   |
| WNK2 CCTL1  | wnk2_cctl1 | Q9Y3S1 (WNK2)      | block map (455–482, 490–544) |
| NRBP1 CCT   | nrbp1_cct  | Q9UHY1 (NRBP1)     | +440                   |
| FERRY CCTL  | ferry_cctl | Q9NQ89 (FERRY3)    | 0                      |
| WNK4 CCTL2  | wnk4_cctl2 | Q96J92 (WNK4)      | +666                   |

Offsets were computed by local pairwise alignment (Biopython `PairwiseAligner`,
match +2 / mismatch −1 / open −10 / extend −0.5) of each PDB chain-A CA sequence
against the canonical UniProt sequence; see `_compute_uniprot_offsets.py`.
All alignments cover ≥99% of the construct.

## Data consolidation (why there is a single bundle)

ALL input data lives in `organized_by_system_bundle_20260327` (per system/replica:
`csv/rmsd.csv`, `csv/rmsf.csv`, `csv/sasa.csv`, `dssp/dssp_per_residue.csv`) —
one canonical group, no loose extracts. Former duplicates were removed after
verification:

- `sasa_only_v2/data/*_residue_summary.csv` — was byte-identical to the bundle's
  `csv/sasa.csv` (pure duplicate) → deleted.
- `rmsf_only_v2/*_rmsf_ca_corrected.csv` — value-identical to the bundle for
  5/6 systems; for OSR1 it held the corrected self-aligned Cα peptide values.
  Those corrected values were **folded into the bundle** (`osr1_cct/<R>/csv/rmsf.csv`
  now stores `protein_ca` + `peptide_ca` chains, same schema as every other
  system) and the loose "corrected" files were deleted.
- `uniprot_offset_map.json` / `replot_manifest.json` — superseded (mapping) /
  duplicate format (manifest) → deleted. `replot_manifest.csv` remains: the
  figure script resolves RMSD through it because 3 replicas (OSR1 R2, WNK2
  R1/R2) need the authoritative runtime time axis recorded there (the bundle's
  native axis differs for those runs).

Verification: regenerating all panels from this consolidated package reproduces
the final figures (7/7 PNG, 600 dpi; systems without numbering changes are
pixel-identical).

## Figure conventions

- **Colors**: Protein/CCT chain = `#90D090` (green); peptide/RFxV ligand = `#F28E2B` (orange).
  DSSP track: α-helix `#E07070`, β-strand `#5090D0`, coil `#D0D0D0`.
- **Numbering**: all per-residue data and figure axes use **canonical UniProt numbering**,
  applied to the data by `_apply_uniprot_numbering.py` from `uniprot_mapping.json`
  (sequence-verified against rest.uniprot.org, Sep 2026):
  - OSR1 CCT → O95747 (OXSR1) +431 — construct = UniProt 433–527; pdb residue 1 is an
    N-terminal Gly expression-tag residue numbered by continuity (432). *(was +432)*
  - SPAK CCT → Q9UEW8 (STK39) +447 (UniProt 448–545).
  - WNK2 CCTL1 → Q9Y3S1 (WNK2) block mapping: pdb 1–28 → 455–482 and pdb 29–83 → 490–544
    (the crystallized construct carries an internal deletion of UniProt 483–489, KLKGKPK).
    *(was +455)*
  - NRBP1 CCT → Q9UHY1 (NRBP1) +440 (UniProt 441–535).
  - FERRY CCTL → Q9NQ89 (FERRY3) +0.
  - WNK4 CCTL2 → Q96J92 (WNK4) +666 (UniProt 667–751). *(was +665)*
- **Ligand numbering** (real UniProt of the parent protein, evidence in `uniprot_mapping.json`):
  - OSR1 ligand GRFQVT (2 copies) → WNK4 Q96J92 1015–1020 (PDBe SIFTS for 2V3S + exact
    canonical match). The paired-chain panel shows copy D = 1015–1019.
  - SPAK ligand KPQLVGRFQVTSSKE → WNK4 Q96J92 1010–1024 (exact unique canonical match).
  - WNK2 ligand LTQVVHSAGRRFIVSPV → WNK1 Q9H4A3 1247–1263 (PDBe SIFTS for 6FBK + exact match).
  - WNK4 CCT2 ligand PAKKKSCFQITSVTT (βϕA / KKKSxFQITSV motif) → TSC22D2 O75157 5–19
    (AF3 model WNK4+TSC22D2 of Magaña-Ávila et al. 2026 + exact unique canonical match).
  - NRBP1 and FERRY ligands EPYRRGRWTCMEYYE (RWTC / βϕB motif region, mediates NRBP1
    binding) → TSC22D2 O75157 173–187 (exact unique canonical match).
- **Panels**: per system, 2 rows (R1, R2) × 3 columns — A: Cα RMSD (0–200 ns),
  B: RMSF + DSSP track, C: mean SASA + DSSP track. Protein raw trace α=0.22 under a
  31/7/5-point moving-average bold line; value-mapped gradient fill under RMSF/SASA.
- **Quality**: PNG at 600 dpi (20×12 in) + SVG vector masters + one DOCX with all panels.

## Reproduce

```powershell
py -3.13 02figure_scripts\_build_md_per_system_panels.py
# → rewrites 01_figures/panels in place (PNG 600 dpi + SVG + overlay + DOCX)
```

Dependencies: `matplotlib`, `numpy`, `pandas`, `python-docx` (plus `mdtraj` and
`plotly` imported by the base module). Figure inputs resolve inside the package
(`Path(__file__).parents[1]`). One machine dependency remains: RMSD time axes are
resolved via `03_data/replot_manifest.csv`, whose `rmsd_csv` entries point at the
original INVENTARIO run folders (`d:\DESCARGAS_PESADAS\...`); if absent, the
script falls back to the bundle's `csv/rmsd.csv` (native axis — differs only for
OSR1 R2 / WNK2 R1–R2).

## Provenance chain

1. MD trajectories (OpenMM) → RMSD/RMSF/SASA CSVs per replica (authoritative runtime
   axes recorded in `replot_manifest.csv`).
2. Corrected self-aligned Cα RMSF → folded into the bundle `csv/rmsf.csv`
   (OSR1 included: `protein_ca` + `peptide_ca` chains, canonical numbering).
3. DSSP per-residue assignments → bundle `dssp/dssp_per_residue.csv`.
4. `_compute_uniprot_offsets.py` → `uniprot_offset_map.json` (historical).
5. `_apply_uniprot_numbering.py` + `uniprot_mapping.json` → canonical UniProt numbering
   in all per-residue data (protein + mapped ligands; originals kept in `*_pdb` columns).
6. `_build_md_per_system_panels.py` → `01_figures/panels` (final paper figures).

## 05 / 06 - run artifacts and the end-to-end example (added 2026-10-19)

- `05_run_artifacts/<system>/` - for each of the 6 systems: `inputs/` (original input model + the
  solvated prepared/equilibrated systems written by the MD runner), `csv/` (real RMSD/RMSF/SASA/DSSP and
  OpenMM energy outputs, R1/R2) and `final_frame/` (solute at 202.5 ns where the final State survives).
- `06_output_example/` - SPAK CCT, 2 x 200 ns, every recovered output (logs, States, checkpoints,
  energies, analysis CSVs, figures, final frames) with the scripts that regenerate the final frames and
  plots and verify them against the originals.
- **Runbook**: `06_output_example/RUN_E2E.md` (input -> preparation/MD -> analysis -> CSV -> final frame
  -> plots), one command per step; `06_output_example/run_e2e.ps1` runs the CPU steps.

## 07 - Residue interaction networks (`07_networks/`)

Residue interaction network analysis for CCT domain / peptide complexes sampled by molecular dynamics.
The pipeline measures which residue contacts hold the domain–peptide interface together and how reproducible
those contacts are across independent trajectory samples:

1. `sample_ensembles.py`: Splits trajectory frames into interleaved ensembles covering the full simulation.
2. Residue interaction network generation via RING 3 (Cytoscape JSON output).
3. `build_consensus.py`: Merges per-ensemble networks into a consensus network with edge frequencies and contact distances.
4. `renumber_network.py`: Filters interface contacts above frequency thresholds and maps CCT domain numbering to canonical UniProt.
5. `select_representative_state.py`: Identifies the individual ensemble member whose interface geometry minimally deviates from the consensus state.

Documentation, software versions, and dependencies: see `07_networks/README.md` and `07_networks/software_versions.tsv`.

### Large files not in this repository

| File | Size | Where |
|---|---|---|
| Production trajectories `replica_*.dcd`, 6 systems x 2 replicas | ~6 GiB (SPAK: 20,000 frames x 27,549 atoms = 6.2 GiB) up to tens of GiB each | authors' storage, on request |

No file in the repository is larger than 100 MB; SHA256 sums: `05_run_artifacts/SHA256SUMS.txt`,
`06_output_example/SHA256SUMS.txt`.

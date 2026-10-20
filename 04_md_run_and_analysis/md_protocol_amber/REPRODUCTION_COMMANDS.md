# Reproduction commands — 5 paper complexes
Canonical script: `00_MD_protocol/runcomplex_paper_amber.py` (ff14SB + TIP3P).
Paper platform: CUDA / RTX 5080. DCD report every 10ps.

## 1. OSR1 CCT + GRFQVT (PDB 2V3S / 2SV3 in text)
Input: OSR1 CCT–GRFQVT crystal structure (2V3S) or `Control_1_bound_OSR1(cct)-WNK4(RfXV)-PDB_2v3s/replica_*/ topology`.
```bash
python 00_MD_protocol/runcomplex_paper_amber.py --pdb 2v3s_complex.pdb --prepare --ns 200 --replica 1 --output_dir runs/OSR1_GRFQVT_R1 --seed 1
python 00_MD_protocol/runcomplex_paper_amber.py --pdb 2v3s_complex.pdb --prepare --ns 200 --replica 2 --output_dir runs/OSR1_GRFQVT_R2 --seed 2
```

## 2. WNK2 CCT1 + LTQVVHSAGRRFIVSPVPESRLR (PDB 6FBK)
Input: `Control_bound_2_WNK2(cctl1)-WNK1(17aa-RFxV)-PDB_6FBK/replica_*/` or 6FBK.
```bash
python 00_MD_protocol/runcomplex_paper_amber.py --pdb 6FBK_complex.pdb --prepare --ns 200 --replica 1 --output_dir runs/WNK2_6FBK_R1 --seed 1
python 00_MD_protocol/runcomplex_paper_amber.py --pdb 6FBK_complex.pdb --prepare --ns 200 --replica 2 --output_dir runs/WNK2_6FBK_R2 --seed 2
```

## 3. SPAK CCT + KPQLVGRFQVTSSKE (AF3 model)
Input: SPAK AF3 model + peptide (see `control_4_SPAK-WNK1PEPTIDE/SPAK_R1-*/`).
```bash
python 00_MD_protocol/runcomplex_paper_amber.py --pdb SPAK_KPQLVGRFQVTSSKE.pdb --prepare --ns 200 --replica 1 --output_dir runs/SPAK_R1 --seed 1
python 00_MD_protocol/runcomplex_paper_amber.py --pdb SPAK_KPQLVGRFQVTSSKE.pdb --prepare --ns 200 --replica 2 --output_dir runs/SPAK_R2 --seed 2
```

## 4. WNK4 CCT2 + PAKKKSCFQITSVTT (novel beta-phi-A motif, AF3 model)
Input: `01_systems/WNK4CCT2_RpA.example.pdb` (copy of `wnk4ROUND2_run/WNK4CCT2_RpA.pdb`).
```bash
python 00_MD_protocol/runcomplex_paper_amber.py --pdb 01_systems/WNK4CCT2_RpA.example.pdb --prepare --ns 200 --replica 1 --output_dir runs/WNK4CCT2_R1 --seed 1
python 00_MD_protocol/runcomplex_paper_amber.py --pdb 01_systems/WNK4CCT2_RpA.example.pdb --prepare --ns 200 --replica 2 --output_dir runs/WNK4CCT2_R2 --seed 2
```
Actually executed reference: `INVENTARIO/wnk4ROUND2_run/runs/replica_1|2/` + `outputs/wnk4round2_r*_sasa_water/`.

## 5. NRBP CCT + EPYRRGRWTCMEYYE (beta-phi-C motif, AF3 model)
Input: `NRBP1_replica1/NRBP1_rep1/` or `NRBP1_replica2/results_nrbp1pf2/`.
```bash
python 00_MD_protocol/runcomplex_paper_amber.py --pdb NRBP_EPYRRGRWTCMEYYE.pdb --prepare --ns 200 --replica 1 --output_dir runs/NRBP_R1 --seed 1
python 00_MD_protocol/runcomplex_paper_amber.py --pdb NRBP_EPYRRGRWTCMEYYE.pdb --prepare --ns 200 --replica 2 --output_dir runs/NRBP_R2 --seed 2
```

## 6. FERRY3 CCT-like + beta-phi-C peptide (Fig.5E, simulated + co-IP only)
Same protocol, system: putative FERRY3 domain + beta-phi-C peptide. Inputs in `novel_ccts_and_WNK4(erendira)/novel_ccts/` (check `bound/` vs `unbound/`).

## Downstream analysis (paper window)
Only 50–150ns per replica: 2 windows (50–100, 100–150) x 3 ensembles of 100 conformations = 6 ensembles/replica → RING → consensus → keep only interactions in both replicas, 20% threshold. See `04_RING_validation/RING_PROTOCOL.md`.
```bash
python 02_analysis_scripts/smic/run_analysis.py --top runs/<SYS>_R1/replica_1_equilibrated.pdb --traj runs/<SYS>_R1/replica_1.dcd --out outputs/<SYS>_R1
```

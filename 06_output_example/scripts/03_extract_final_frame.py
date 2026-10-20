"""Step 3 - Final frame: OpenMM final State XML (+ the run's solvated topology PDB) -> solute-only PDB.

Usage:
  python extract_final_frame.py --state replica_1_final.xml --topology replica_1_equilibrated.pdb --out R1_final_frame_solute.pdb

The state XML is the one written by runcomplex_paper_amber.py at the end of production
(simulation.saveState(<run>_final.xml)). Positions are imaged so the protein is whole and the
peptide is placed in the periodic image closest to it; water and Na+/Cl- ions are removed.
"""
import argparse, os, numpy as np, openmm as mm, mdtraj as md

ap = argparse.ArgumentParser()
ap.add_argument("--state", required=True)
ap.add_argument("--topology", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--keep_solvent", action="store_true")
a = ap.parse_args()

with open(a.state) as f:
    st = mm.XmlSerializer.deserialize(f.read())
pos = st.getPositions(asNumpy=True).value_in_unit(mm.unit.nanometer)
t = md.load(a.topology)
if t.n_atoms != pos.shape[0]:
    raise SystemExit(f"atom count mismatch: topology {t.n_atoms} vs state {pos.shape[0]}")
t.xyz = np.asarray(pos, dtype=np.float32)[None]
box = st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(mm.unit.nanometer)
t.unitcell_vectors = np.asarray(box, dtype=np.float32)[None]
t.time = np.array([st.getTime().value_in_unit(mm.unit.picosecond)])
mols = t.topology.find_molecules()
prot = [m for m in mols if any(at.residue.is_protein for at in m)]
anchor = max(prot, key=len)
t = t.image_molecules(anchor_molecules=[set(anchor)], other_molecules=[set(m) for m in mols if m is not anchor], inplace=False)
if not a.keep_solvent:
    t = t.atom_slice(t.topology.select("not water and not resname NA and not resname CL"))
os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
t.save_pdb(a.out)
print(f"time_ns={t.time[0]/1000:.4f} atoms_out={t.n_atoms} box_nm={np.diag(box).round(4).tolist()} -> {a.out}")


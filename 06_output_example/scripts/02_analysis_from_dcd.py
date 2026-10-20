"""Step 2 - Trajectory analysis (RMSD / RMSF / per-residue SASA) with the paper's SMIC functions.

Uses compute_all_rmsd / compute_all_rmsf from 04_md_run_and_analysis/md_analysis/analysis_rmsd.py and
MDTraj Shrake-Rupley (as analysis_water.py), stride 5 over the 10-ps DCD (= 50 ps, 4000 frames / 200 ns),
which is the stride recorded in the original analysis logs.

Usage:
  python 02_analysis_from_dcd.py --top R1/md/equilibrated_solvated.pdb --traj replica_1.dcd --out R1/analysis_rerun --stride 5
  python 02_analysis_from_dcd.py --top R1/md/equilibrated_solvated.pdb --state R1/md/equilibrated_state_or_final.xml ...   (smoke test without DCD)

NOTE: the 200-ns DCDs (~6.2 GiB each) are not shipped (see README). The column layout of this script's CSVs
is a plain dump of the SMIC arrays and is not byte-identical to the curated analysis_csv/*.csv.
"""
import argparse, os, sys, numpy as np, pandas as pd, mdtraj as md
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "_smic_support"))  # runtime_config.py (missing from 04, copied from the SMIC suite)
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "04_md_run_and_analysis", "md_analysis")))
from analysis_rmsd import compute_all_rmsd, compute_all_rmsf  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--top", required=True); ap.add_argument("--traj", nargs="*", default=[])
ap.add_argument("--state", nargs="*", default=[], help="OpenMM State XMLs appended as frames (smoke test)")
ap.add_argument("--stride", type=int, default=5); ap.add_argument("--out", required=True)
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
top = md.load(a.top)
frames = [top]
for t in a.traj:
    frames.append(md.load(t, top=a.top, stride=a.stride))
if a.state:
    import openmm as mm
    for s in a.state:
        st = mm.XmlSerializer.deserialize(open(s).read())
        f = top[:]; f.xyz = np.asarray(st.getPositions(asNumpy=True).value_in_unit(mm.unit.nanometer), dtype=np.float32)[None]
        f.unitcell_vectors = np.asarray(st.getPeriodicBoxVectors(asNumpy=True).value_in_unit(mm.unit.nanometer), dtype=np.float32)[None]
        frames.append(f)
traj = md.join(frames[1:]) if a.traj else md.join(frames)
traj.image_molecules(inplace=True)
traj = traj.atom_slice(traj.topology.select("not water and not resname NA and not resname CL"))
r = compute_all_rmsd(traj)
pd.DataFrame({k: np.asarray(v) for k, v in r.items() if np.ndim(v) == 1 and len(v) == traj.n_frames}).to_csv(os.path.join(a.out, "rmsd_rerun.csv"), index=False)
rows = []
for key, val in compute_all_rmsf(traj).items():
    arrs = [np.asarray(x) for x in val]
    for i in range(len(arrs[0])):
        rows.append([key] + [x[i] if i < len(x) else None for x in arrs])
pd.DataFrame(rows).to_csv(os.path.join(a.out, "rmsf_rerun.csv"), index=False, header=False)
sasa = md.shrake_rupley(traj, mode="residue")
res = list(traj.topology.residues)
pd.DataFrame({"resSeq": [x.resSeq for x in res], "resname": [x.name for x in res], "chain": [x.chain.index for x in res],
              "mean_sasa_nm2": sasa.mean(0), "std_sasa_nm2": sasa.std(0)}).to_csv(os.path.join(a.out, "sasa_rerun.csv"), index=False)
print("frames:", traj.n_frames, "->", a.out, "| rmsd keys:", list(r.keys()))


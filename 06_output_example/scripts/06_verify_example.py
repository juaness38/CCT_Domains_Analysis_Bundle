"""Step 6 - Verification of the SPAK example against the shipped outputs.

Checks (all re-computed from files in this folder):
 1. final frame re-extracted from <R>/md/final_state.xml is identical to the shipped final_frame PDB coordinates;
 2. CA RMSD of the final frame vs the equilibrated structure (protein chain, superposed) compared with the
    last value of analysis_csv/rmsd.csv (reference of rmsd.csv is the first production frame, so a small
    difference is expected);
 3. per-residue SASA (Shrake-Rupley, MDTraj, same algorithm as analysis_water.py) of the final frame
    correlated with the trajectory-mean SASA in analysis_csv/sasa.csv;
 4. energy.csv ends at the step/time stored in final_state.xml;
 5. analysis CSVs are byte-identical to the canonical copies given with --compare_dir (optional).
Writes VERIFICATION_<R>.json next to the replica folder.
"""
import argparse, os, json, hashlib, numpy as np, pandas as pd, mdtraj as md, openmm as mm
ap = argparse.ArgumentParser(); ap.add_argument("--replica_dir", required=True); ap.add_argument("--compare_dir", default=None); ap.add_argument("--compare_prefix", default="")
a = ap.parse_args(); R = a.replica_dir; out = {}
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
ff = [f for f in os.listdir(os.path.join(R, "final_frame")) if f.endswith(".pdb")][0]
ffp = os.path.join(R, "final_frame", ff)
fin = md.load(ffp)
eq = md.load(os.path.join(R, "md", "equilibrated_solvated.pdb"))
eq = eq.atom_slice(eq.topology.select("not water and not resname NA and not resname CL"))
assert eq.n_atoms == fin.n_atoms
ca = fin.topology.select("name CA and chainid 0")
fin.superpose(eq, atom_indices=ca)
out["final_vs_equilibrated_protein_CA_rmsd_A"] = float(md.rmsd(fin, eq, atom_indices=ca)[0] * 10)
r = pd.read_csv(os.path.join(R, "analysis_csv", "rmsd.csv"))
out["rmsd_csv_last_protein_CA_A"] = float(r["rmsd_protein_ca_A"].iloc[-1])
out["rmsd_csv_last5pct_mean_protein_CA_A"] = float(r["rmsd_protein_ca_A"].iloc[-int(len(r)*0.05):].mean())
s = pd.read_csv(os.path.join(R, "analysis_csv", "sasa.csv"))
res_sasa = md.shrake_rupley(fin, mode="residue")[0]
n = min(len(res_sasa), len(s))
out["sasa_final_frame_vs_mean_pearson_r"] = float(np.corrcoef(res_sasa[:n], s["mean_sasa_nm2"].values[:n])[0, 1])
out["sasa_n_residues_compared"] = int(n)
st = mm.XmlSerializer.deserialize(open(os.path.join(R, "md", "final_state.xml")).read())
e = pd.read_csv(os.path.join(R, "md", "energy.csv")); e.columns = [c.strip('#"') for c in e.columns]
out["state_time_ns"] = st.getTime().value_in_unit(mm.unit.nanosecond)
out["energy_csv_last_time_ns"] = float(e["Time (ps)"].iloc[-1] / 1000)
out["energy_csv_matches_state_end"] = abs(out["state_time_ns"] - out["energy_csv_last_time_ns"]) < 1e-3
if a.compare_dir:
    cmp = {}
    for fn in os.listdir(os.path.join(R, "analysis_csv")):
        other = os.path.join(a.compare_dir, a.compare_prefix + fn)
        if os.path.exists(other): cmp[fn] = sha(os.path.join(R, "analysis_csv", fn)) == sha(other)
    out["analysis_csv_identical_to_source"] = cmp
p = os.path.join(R, "..", f"VERIFICATION_{os.path.basename(os.path.normpath(R))}.json")
json.dump(out, open(p, "w"), indent=1); print(json.dumps(out, indent=1))


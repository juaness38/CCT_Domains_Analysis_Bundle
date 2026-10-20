"""Step 5 - Plots: regenerate per-replica figures from the shipped CSVs (no trajectory needed).

Usage: python plot_outputs.py --replica_dir ../R1 --label "SPAK CCT R1"
Reads  <replica_dir>/analysis_csv/{rmsd,rmsf,sasa}.csv and <replica_dir>/md/energy.csv
Writes <replica_dir>/figures/regenerated_{rmsd,rmsf,sasa,energy}.png
"""
import argparse, os, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ap = argparse.ArgumentParser(); ap.add_argument("--replica_dir", required=True); ap.add_argument("--label", default="")
a = ap.parse_args()
A = os.path.join(a.replica_dir, "analysis_csv"); F = os.path.join(a.replica_dir, "figures"); os.makedirs(F, exist_ok=True)
r = pd.read_csv(os.path.join(A, "rmsd.csv"))
n = len(r); t = [i * 200.0 / n for i in range(n)]  # 4000 frames = DCD stride 5 x 10 ps = 50 ps; native time column is x5 mis-scaled (see README)
fig, ax = plt.subplots(figsize=(8, 3.5))
ax.plot(t, r["rmsd_protein_ca_A"], lw=0.6, color="#3b8f3b", label="CCT Calpha")
ax.plot(t, r["rmsd_peptide_backbone_A"], lw=0.6, color="#F28E2B", label="peptide backbone")
ax.set_xlabel("time (ns)"); ax.set_ylabel("RMSD (A)"); ax.set_title(f"{a.label} RMSD"); ax.legend(); fig.tight_layout(); fig.savefig(os.path.join(F, "regenerated_rmsd.png"), dpi=200); plt.close(fig)
f = pd.read_csv(os.path.join(A, "rmsf.csv"))
fig, ax = plt.subplots(figsize=(8, 3.5))
for ch, g in f.groupby("chain", sort=False):
    ax.plot(range(len(g)), g["rmsf_A"], marker=".", lw=0.8, label=str(g["chain_display"].iloc[0]))
ax.set_xlabel("residue index (per chain)"); ax.set_ylabel("RMSF (A)"); ax.set_title(f"{a.label} RMSF"); ax.legend(); fig.tight_layout(); fig.savefig(os.path.join(F, "regenerated_rmsf.png"), dpi=200); plt.close(fig)
s = pd.read_csv(os.path.join(A, "sasa.csv"))
fig, ax = plt.subplots(figsize=(8, 3.5))
for role, g in s.groupby("entity_role", sort=False):
    ax.errorbar(g["resSeq"], g["mean_sasa_nm2"], yerr=g["std_sasa_nm2"], fmt=".", ms=3, lw=0.5, label=role)
ax.set_xlabel("residue (resSeq)"); ax.set_ylabel("mean SASA (nm^2)"); ax.set_title(f"{a.label} per-residue SASA"); ax.legend(); fig.tight_layout(); fig.savefig(os.path.join(F, "regenerated_sasa.png"), dpi=200); plt.close(fig)
e = pd.read_csv(os.path.join(a.replica_dir, "md", "energy.csv"))
e.columns = [c.strip('#"') for c in e.columns]
fig, axs = plt.subplots(3, 1, figsize=(8, 7), sharex=True)
x = e["Time (ps)"] / 1000.0
for axx, col in zip(axs, ["Potential Energy (kJ/mole)", "Temperature (K)", "Density (g/mL)"]):
    axx.plot(x, e[col], lw=0.4); axx.set_ylabel(col)
axs[-1].set_xlabel("time (ns)"); axs[0].set_title(f"{a.label} OpenMM StateData"); fig.tight_layout(); fig.savefig(os.path.join(F, "regenerated_energy.png"), dpi=200); plt.close(fig)
print("plots ->", F)


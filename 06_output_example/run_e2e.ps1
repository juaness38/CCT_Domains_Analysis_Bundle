# run_e2e.ps1 - SPAK CCT end-to-end reproduction.
# Steps 3-6 (final frame, analysis smoke test, plots, verification, checksums) run locally on CPU.
# Step 1 (MD, GPU) runs only with -RunMD; step 2 (trajectory analysis) runs only if the DCDs are given.
# Usage (from this folder):
#   powershell -ExecutionPolicy Bypass -File run_e2e.ps1 -Python <python> [-RunMD] [-Dcd1 R1.dcd] [-Dcd2 R2.dcd]
param([string]$Python = "python", [switch]$RunMD, [string]$Dcd1 = "", [string]$Dcd2 = "")
$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot
function Run([string]$Name, [string[]]$PyArgs) {
  Write-Host "== $Name"
  & $Python @PyArgs
  if ($LASTEXITCODE -ne 0) { throw "$Name failed (exit $LASTEXITCODE)" }
}
if ($RunMD) {
  # STEP 1 - MD. Arguments recovered from R1/R2 md/simulation.log (seed 1001/1002, 200 ns, --prepare).
  foreach ($n in 1, 2) {
    Run "01 MD replica $n" @("-u", "scripts/01_run_md__runcomplex_paper_amber.py", "--pdb", "inputs/spak_cct_input_stk39-wnk4-relax.pdb",
      "--output_dir", "rerun_md/replica_$n", "--replica", "$n", "--ns", "200", "--seed", "100$n", "--platform", "CUDA", "--gpu_id", "0", "--prepare")
  }
}
foreach ($R in "R1", "R2") {
  $dcd = if ($R -eq "R1") { $Dcd1 } else { $Dcd2 }
  if ($dcd) {
    Run "02 analysis $R" @("scripts/02_analysis_from_dcd.py", "--top", "$R/md/equilibrated_solvated.pdb", "--traj", $dcd, "--stride", "5", "--out", "$R/analysis_rerun")
  }
  Run "03 final frame $R" @("scripts/03_extract_final_frame.py", "--state", "$R/md/final_state.xml", "--topology", "$R/md/equilibrated_solvated.pdb",
    "--out", "$R/final_frame/spak_cct_${R}_final_frame_202.5ns_solute.pdb")
  Run "04 analysis smoke test $R" @("scripts/02_analysis_from_dcd.py", "--top", "$R/md/equilibrated_solvated.pdb", "--state", "$R/md/final_state.xml", "--out", "$R/analysis_smoketest")
  Run "05 plots $R" @("scripts/05_plot_outputs.py", "--replica_dir", $R, "--label", "SPAK CCT $R")
  Run "06 verify $R" @("scripts/06_verify_example.py", "--replica_dir", $R, "--compare_dir", "../05_run_artifacts/spak_cct/csv", "--compare_prefix", "${R}_")
}
Write-Host "== SHA256SUMS.txt"
Get-ChildItem -Recurse -File | Where-Object { $_.Name -ne "SHA256SUMS.txt" } | Sort-Object FullName | ForEach-Object {
  $rel = $_.FullName.Substring($PSScriptRoot.Length + 1).Replace("\", "/")
  "{0}  {1}" -f (Get-FileHash $_.FullName -Algorithm SHA256).Hash.ToLower(), $rel } | Set-Content SHA256SUMS.txt -Encoding ascii
Write-Host "E2E done."

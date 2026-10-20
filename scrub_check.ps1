# scrub_check.ps1 - run before every commit/push of this public repository.
# Fails (exit 1) on secrets, infrastructure/orchestration strings, or files > 100 MB.
param([string]$Root = $PSScriptRoot)
$secret = 'gh[opsu]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}|sk-[A-Za-z0-9]{20,}|AKIA[0-9A-Z]{16}|-----BEGIN [A-Z ]*PRIVATE KEY|xox[abp]-[A-Za-z0-9-]{10,}|(password|passwd|api[_-]?key|secret[_-]?key)\s*[:=]\s*\S+'
$infra  = '\b(?:\d{1,3}\.){3}\d{1,3}:\d+\b|\bssh\s+-|\bscp\s+-|root@|/workspace/|/root/|vast\.ai|vast_key|runpod|salad\s?cloud|lambda\s?labs|nvidia-mps|CUDA_MPS_PIPE|\.ssh[\\/]'
$bad = @()
Get-ChildItem $Root -Recurse -File -Force | Where-Object { $_.FullName -notmatch '\\\.git\\' -and $_.Name -ne 'scrub_check.ps1' } | ForEach-Object {
  if ($_.Length -gt 100MB) { $bad += "TOO LARGE ($([math]::Round($_.Length/1MB)) MB): $($_.FullName)" }
  if ($_.Name -match '^\.env|\.pem$|\.key$|id_rsa|id_ed25519') { $bad += "SENSITIVE FILE NAME: $($_.FullName)" }
  if ($_.Extension -in '.md','.py','.ps1','.txt','.csv','.json','.log','.yml','.yaml','.pdb','.xml','.cfg','.ini','.toml','.sh') {
    foreach ($m in (Select-String -LiteralPath $_.FullName -Pattern $secret, $infra -AllMatches)) { $bad += "$($_.FullName):$($m.LineNumber): $($m.Line.Trim().Substring(0,[math]::Min(120,$m.Line.Trim().Length)))" }
  }
}
if ($bad) { $bad | Select-Object -First 200; Write-Host "SCRUB FAILED: $($bad.Count) finding(s)"; exit 1 } else { Write-Host "SCRUB OK: no secrets, infrastructure strings or >100 MB files" }

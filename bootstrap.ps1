$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $root
$expectedPython = "3.12.13"
$python = Get-Command python -ErrorAction Stop
$version = & $python.Source -c "import platform; print(platform.python_version())"
if ($version -ne $expectedPython) { throw "Python $expectedPython required; found $version at $($python.Source)" }
if (-not (Test-Path .venv)) { & $python.Source -m venv .venv }
$venvPython = Join-Path $root ".venv\Scripts\python.exe"
& $venvPython -m pip install --upgrade pip setuptools wheel
& $venvPython -m pip install -r environment\requirements-core.txt
& $venvPython -m pip install -r environment\requirements-science.txt
$cambDir = Join-Path $root "vendor\CAMB"
if (-not (Test-Path $cambDir)) { New-Item -ItemType Directory -Force -Path (Split-Path -Parent $cambDir) | Out-Null; git clone https://github.com/cmbant/CAMB.git $cambDir }
Push-Location $cambDir
try {
  git fetch --all --tags
  git checkout --detach 3ef0272d6f7ba1231128872e56e6d4c12af8267b
  $head = (git rev-parse HEAD).Trim()
  if ($head -ne "3ef0272d6f7ba1231128872e56e6d4c12af8267b") { throw "CAMB commit mismatch: $head" }
  & $venvPython -m pip install .
} finally { Pop-Location }
& $venvPython scripts\verify_artifact.py
Write-Host "Bootstrap complete. Scientific R1 acceptance is NOT implied; run doctor + validation on this Windows host."

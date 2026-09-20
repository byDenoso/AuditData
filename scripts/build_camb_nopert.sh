#!/usr/bin/env bash
set -euo pipefail

ROOT=${1:-.cache/camb-1.6.6-nopert}
VENV=${2:-.venv-nopert}
PROVENANCE=${3:-results/n3_spectra/camb_nopert_provenance.json}

rm -rf "$ROOT" "$VENV"
git clone --recursive --branch 1.6.6 --depth 1 https://github.com/cmbant/CAMB.git "$ROOT"
SOURCE_COMMIT=$(git -C "$ROOT" rev-parse HEAD)

python - "$ROOT/fortran/DarkEnergyQuintessence.f90" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
text = path.read_text(encoding="utf-8")
start = text.index("subroutine TQuintessence_Init")
end = text.index("end subroutine  TQuintessence_Init", start)
block = text[start:end]
needle = "this%num_perturb_equations = 2"
if block.count(needle) != 1:
    raise SystemExit(f"expected one perturbation-equation assignment, found {block.count(needle)}")
patched = block.replace(needle, "this%num_perturb_equations = 0 ! exact-background perturbation ablation")
path.write_text(text[:start] + patched + text[end:], encoding="utf-8")
PY

grep -n "exact-background perturbation ablation" "$ROOT/fortran/DarkEnergyQuintessence.f90"
git -C "$ROOT" diff --check
PATCH_SHA=$(sha256sum "$ROOT/fortran/DarkEnergyQuintessence.f90" | awk '{print $1}')

python -m venv "$VENV"
"$VENV/bin/python" -m pip install --upgrade pip wheel
"$VENV/bin/python" -m pip install 'numpy==2.1.3' 'scipy==1.14.1' 'matplotlib==3.9.2'
"$VENV/bin/python" -m pip install -e "$ROOT"
"$VENV/bin/python" - <<'PY'
import camb
from pathlib import Path
assert camb.__version__ == "1.6.6"
assert "camb-1.6.6-nopert" in str(Path(camb.__file__).resolve())
print("NOPERT_CAMB_IMPORT_PASS", camb.__version__, Path(camb.__file__).resolve())
PY

mkdir -p "$(dirname "$PROVENANCE")"
cat > "$PROVENANCE" <<JSON
{
  "upstream": "https://github.com/cmbant/CAMB",
  "tag": "1.6.6",
  "source_commit": "$SOURCE_COMMIT",
  "patched_file": "fortran/DarkEnergyQuintessence.f90",
  "patched_file_sha256": "$PATCH_SHA",
  "change": "TQuintessence_Init num_perturb_equations 2 -> 0",
  "scientific_scope": "Exact-background scalar perturbation ablation only"
}
JSON

from pathlib import Path
import json,sys
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'src'))
from native_windows_camb.contract import load_json,validate_manifest,validate_reference_values

def test_manifest_rejects_absolute_windows_paths(tmp_path):
    manifest={"python":{"version":"3.12.13","executable":r"C:\\Users\\Dener\\python.exe"}}
    path=tmp_path/'manifest.json'; path.write_text(json.dumps(manifest),encoding='utf-8'); errors=validate_manifest(load_json(path))
    assert any('absolute/local path' in error for error in errors)

def test_manifest_accepts_relocatable_tokens():
    assert validate_manifest({"schema":"ascom-00323-environment/v1","python":{"version":"3.12.13","executable":"${PYTHON}"},"runtime_root":"${RUNTIME_ROOT}"})==[]

def test_reference_values_require_provenance_and_status():
    values={"m0_loglike":{"value":-395.4831594065799,"status":"historical_pass_needs_r1_replay","provenance":"PEER_CLOSURE_EXECUTION_20260729_v2.md"}}
    assert validate_reference_values(values)==[]

def test_reference_values_reject_unqualified_verified_claim():
    errors=validate_reference_values({"x":{"value":1.0,"status":"verified"}}); assert any('provenance' in e for e in errors)

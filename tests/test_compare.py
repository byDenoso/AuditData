from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'src'))
from native_windows_camb.compare import compare_scalar,compare_vector

def test_scalar_absolute_tolerance_passes_at_boundary(): assert compare_scalar(1.000000001,1.0,abs_tol=1e-9,rel_tol=0.0).passed
def test_scalar_fails_outside_tolerance(): assert not compare_scalar(1.00001,1.0,abs_tol=1e-9,rel_tol=0.0).passed
def test_vector_reports_max_abs_and_rms():
    result=compare_vector([1.0,2.0],[1.0,2.000000001],abs_tol=2e-9,rel_tol=0.0); assert result.passed and result.max_abs>0 and result.rms>0

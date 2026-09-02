from __future__ import annotations
from dataclasses import dataclass
import math
from typing import Iterable

@dataclass(frozen=True)
class ScalarComparison:
    passed: bool; observed: float; reference: float; abs_error: float; rel_error: float; abs_tol: float; rel_tol: float
@dataclass(frozen=True)
class VectorComparison:
    passed: bool; max_abs: float; max_rel: float; rms: float; n: int; abs_tol: float; rel_tol: float

def _threshold(reference: float, abs_tol: float, rel_tol: float) -> float:
    return abs_tol + rel_tol * abs(reference)

def compare_scalar(observed: float, reference: float, *, abs_tol: float, rel_tol: float) -> ScalarComparison:
    observed=float(observed); reference=float(reference); abs_error=abs(observed-reference)
    rel_error=abs_error/abs(reference) if reference!=0 else (0.0 if abs_error==0 else math.inf)
    passed=abs_error <= _threshold(reference,abs_tol,rel_tol)+1e-15
    return ScalarComparison(passed,observed,reference,abs_error,rel_error,abs_tol,rel_tol)

def compare_vector(observed: Iterable[float], reference: Iterable[float], *, abs_tol: float, rel_tol: float) -> VectorComparison:
    obs=[float(x) for x in observed]; ref=[float(x) for x in reference]
    if len(obs)!=len(ref) or not obs: raise ValueError("vectors must be non-empty and have equal length")
    diffs=[abs(a-b) for a,b in zip(obs,ref)]; rels=[d/abs(b) if b!=0 else (0.0 if d==0 else math.inf) for d,b in zip(diffs,ref)]
    passed=all(d<=_threshold(b,abs_tol,rel_tol)+1e-15 for d,b in zip(diffs,ref))
    rms=math.sqrt(math.fsum(d*d for d in diffs)/len(diffs))
    return VectorComparison(passed,max(diffs),max(rels),rms,len(diffs),abs_tol,rel_tol)

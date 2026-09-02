from __future__ import annotations
import hashlib,json,math
from typing import Any

def _normalize(value: Any) -> Any:
    if value is None or isinstance(value,(str,bool,int)): return value
    if isinstance(value,float):
        if not math.isfinite(value): raise ValueError("scientific fingerprints forbid NaN/inf")
        return {"__float_hex__":value.hex()}
    if isinstance(value,dict): return {str(k):_normalize(value[k]) for k in sorted(value,key=str)}
    if isinstance(value,(list,tuple)): return [_normalize(x) for x in value]
    raise TypeError(f"unsupported fingerprint value: {type(value).__name__}")

def scientific_fingerprint(contract: dict[str,Any]) -> str:
    payload=json.dumps(_normalize(contract),sort_keys=True,separators=(",",":"),allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()

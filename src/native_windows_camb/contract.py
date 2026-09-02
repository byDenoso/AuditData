from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

_ABSOLUTE_WINDOWS = re.compile(r"^[A-Za-z]:[\\/]")
_ALLOWED_STATUSES = {"historical_pass_needs_r1_replay","pending_r1","verified","failed","not_applicable"}

def load_json(path: str | Path) -> Any:
    return json.loads(Path(path).read_text(encoding="utf-8"))

def _walk(value: Any, path: str = "$()"):
    if isinstance(value, dict):
        for key, item in value.items(): yield from _walk(item, f"{path}.{key}")
    elif isinstance(value, list):
        for i, item in enumerate(value): yield from _walk(item, f"{path}[{i}]")
    else: yield path, value

def validate_manifest(manifest: dict[str, Any]) -> list[str]:
    errors=[]
    for path,value in _walk(manifest):
        if isinstance(value,str):
            if _ABSOLUTE_WINDOWS.match(value) or value.lower().startswith("file:///"):
                errors.append(f"{path}: absolute/local path is forbidden in public manifests")
            if "\\users\\dener" in value.lower() or "d:\\codex" in value.lower():
                errors.append(f"{path}: workstation-specific path is forbidden")
    return errors

def validate_reference_values(values: dict[str, Any]) -> list[str]:
    errors=[]
    for name,spec in values.items():
        if not isinstance(spec,dict): errors.append(f"{name}: entry must be an object"); continue
        status=spec.get("status")
        if status not in _ALLOWED_STATUSES: errors.append(f"{name}: invalid status {status!r}")
        if "value" not in spec: errors.append(f"{name}: value is required")
        if not spec.get("provenance"): errors.append(f"{name}: provenance is required")
    return errors

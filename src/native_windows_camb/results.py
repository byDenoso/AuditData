from __future__ import annotations
from typing import Iterable, Mapping, Any

def publication_ready(matrix: Iterable[Mapping[str, Any]]) -> tuple[bool,list[str]]:
    missing=[str(row["gate_id"]) for row in matrix if bool(row.get("required_for_paper")) and row.get("status")!="verified"]
    return (not missing,missing)

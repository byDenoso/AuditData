from __future__ import annotations
import hashlib
from pathlib import Path
from typing import Mapping

def sha256_file(path: str|Path) -> str:
    digest=hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""): digest.update(chunk)
    return digest.hexdigest()

def verify_hashes(root: str|Path, expected: Mapping[str,str]) -> tuple[bool,list[str]]:
    root=Path(root); errors=[]
    for relative,digest in expected.items():
        path=root/relative
        if not path.is_file(): errors.append(f"missing: {relative}"); continue
        observed=sha256_file(path)
        if observed.lower()!=str(digest).lower(): errors.append(f"sha256 mismatch: {relative}: expected={digest} observed={observed}")
    return (not errors,errors)

from __future__ import annotations
from typing import Any,Mapping

def build_paper_summary(results: Mapping[str,Mapping[str,Any]]) -> dict[str,Any]:
    verified={k:dict(v) for k,v in results.items() if v.get("status")=="verified"}
    nonverified={k:dict(v) for k,v in results.items() if v.get("status")!="verified"}
    return {"schema":"ascom-00323-paper-summary/v1","publication_ready":not nonverified,"verified_results":verified,"nonverified_results":nonverified}

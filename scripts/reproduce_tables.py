#!/usr/bin/env python3
from __future__ import annotations
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'src'))
from native_windows_camb.replay import build_paper_summary
def main():
    results={}
    for path in sorted((ROOT/'results/r1').glob('G*.json')):
        data=json.loads(path.read_text(encoding='utf-8')); results[data.get('gate_id',path.stem)]=data
    summary=build_paper_summary(results); out=ROOT/'results/generated/paper_summary.json'; out.parent.mkdir(parents=True,exist_ok=True); out.write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    table=ROOT/'results/generated/validation_table.md'; lines=['| Gate | Status | Observed |','|---|---|---|']
    for k,v in sorted(results.items()): lines.append(f"| {k} | {v.get('status','unknown')} | `{json.dumps(v.get('observed'))}` |")
    table.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    if not summary['publication_ready']: print('Paper table generated as audit surface only; publication-ready promotion is false.'); return 2
    print(f'PASS: {table}'); return 0
if __name__=='__main__': raise SystemExit(main())

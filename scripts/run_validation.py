#!/usr/bin/env python3
from __future__ import annotations
import argparse,csv,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'src'))
from native_windows_camb.results import publication_ready
def load_matrix():
    rows=[]
    with (ROOT/'validation/validation_matrix.csv').open(newline='',encoding='utf-8') as f:
        for r in csv.DictReader(f):
            r['required_for_paper']=r['required_for_paper'].lower()=='true'; receipt=ROOT/r['r1_receipt']
            if receipt.is_file():
                try:
                    data=json.loads(receipt.read_text(encoding='utf-8'))
                    if data.get('gate_id')==r['gate_id'] and data.get('status')=='verified': r['status']='verified'
                    elif data.get('status')=='failed': r['status']='failed'
                except Exception: r['status']='failed'
            rows.append(r)
    return rows
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--strict',action='store_true'); a=ap.parse_args(); rows=load_matrix(); ready,missing=publication_ready(rows)
    for r in rows: print(f"{r['gate_id']:>3} {r['status']:<32} {r['description']}")
    print(f"\nPUBLICATION_READY={str(ready).lower()}")
    if missing: print('UNVERIFIED_REQUIRED_GATES='+','.join(missing))
    return 1 if a.strict and not ready else 0
if __name__=='__main__': raise SystemExit(main())

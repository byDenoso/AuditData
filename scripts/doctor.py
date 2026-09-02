#!/usr/bin/env python3
from __future__ import annotations
import argparse,importlib,json,platform
from datetime import datetime,timezone
from pathlib import Path
EXPECTED_PY='3.12.13'; EXPECTED={'cobaya':'3.6.2','mpi4py':'4.1.2','sacc':'2.1.2','getdist':'1.7.7'}
def package_version(name):
    try:
        m=importlib.import_module(name); return getattr(m,'__version__','unknown')
    except Exception as exc: return f'IMPORT_ERROR:{type(exc).__name__}:{exc}'
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--strict',action='store_true'); ap.add_argument('--output',type=Path); a=ap.parse_args()
    receipt={'schema':'ascom-00323-doctor/v1','timestamp_utc':datetime.now(timezone.utc).isoformat(),'platform':{'system':platform.system(),'release':platform.release(),'machine':platform.machine()},'python':{'version':platform.python_version(),'executable':'${PYTHON}'},'packages':{k:package_version(k) for k in EXPECTED}}
    receipt['camb_version']=package_version('camb'); errors=[]
    if platform.system()!='Windows': errors.append('native-Windows acceptance requires Windows')
    if platform.python_version()!=EXPECTED_PY: errors.append(f'Python {EXPECTED_PY} required')
    for k,v in EXPECTED.items():
        if receipt['packages'][k]!=v: errors.append(f"{k} expected {v}; observed {receipt['packages'][k]}")
    if str(receipt['camb_version']).startswith('IMPORT_ERROR'): errors.append('CAMB import failed')
    receipt['errors']=errors; receipt['status']='pass' if not errors else 'fail'; text=json.dumps(receipt,indent=2,sort_keys=True)+'\n'
    if a.output: a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(text,encoding='utf-8')
    print(text,end=''); return 1 if a.strict and errors else 0
if __name__=='__main__': raise SystemExit(main())

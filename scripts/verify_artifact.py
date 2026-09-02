#!/usr/bin/env python3
from __future__ import annotations
import argparse,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'src'))
from native_windows_camb.contract import load_json,validate_manifest,validate_reference_values
from native_windows_camb.integrity import sha256_file
REQUIRED=['README.md','REPRODUCING.md','LICENSE','CITATION.cff','manifests/environment.json','manifests/content-identities.json','validation/reference_values.json','validation/tolerances.json','validation/validation_matrix.csv']
def main():
    p=argparse.ArgumentParser(); p.add_argument('--external-data',type=Path); args=p.parse_args(); errors=[]
    for rel in REQUIRED:
        if not (ROOT/rel).is_file(): errors.append(f'missing repository file: {rel}')
    for rel in ['manifests/environment.json','manifests/external-components.json']: errors.extend(validate_manifest(load_json(ROOT/rel)))
    errors.extend(validate_reference_values(load_json(ROOT/'validation/reference_values.json')))
    if args.external_data:
        ext=load_json(ROOT/'manifests/external-components.json')
        for item in ext['components']:
            digest=item.get('expected_sha256')
            if not digest: continue
            path=args.external_data/item['id']
            if not path.is_file(): errors.append(f"missing external payload: {item['id']}")
            elif sha256_file(path)!=digest: errors.append(f"external payload hash mismatch: {item['id']}")
    if errors:
        for e in errors: print(f'FAIL: {e}')
        return 1
    print('PASS: repository contract and public manifests are internally valid'); return 0
if __name__=='__main__': raise SystemExit(main())

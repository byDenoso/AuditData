#!/usr/bin/env python3
import argparse, json, os, shutil, subprocess, sys
from pathlib import Path
import pandas as pd

TH_VOID = {
    'NEXO-V145': 15.667543,
    'NEXO-V14': 21.494811,
    'NEXO-V11': 22.019039,
    'NEXO-V102': 16.802955,
}
TH_WALL = {'NEXO-W13': 22.126411}
BASE = 'https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/mocks/EZmock/dark/v1'
FINDER = 'https://raw.githubusercontent.com/byDenoso/AuditData/nexo-desi-voidscan-temp-20260826/tools/desi_blind_scan.py'

def run(cmd, **kw):
    print('+', ' '.join(map(str,cmd)), flush=True)
    subprocess.run(cmd, check=True, **kw)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--candidates')
    ap.add_argument('--v15')
    ap.add_argument('--v12')
    ap.add_argument('--eboss')
    ap.add_argument('--nulls',type=int,default=100)
    ap.add_argument('--outdir',default='battery')
    args=ap.parse_args()
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True)
    work=Path('ezwork'); work.mkdir(exist_ok=True)
    finder=work/'desi_blind_scan.py'
    run(['curl','-fL','--retry','5',FINDER,'-o',str(finder)])
    rows=[]
    for m in range(1,6):
        print(f'=== EZMOCK {m}/5 ===', flush=True)
        md=work/f'mock{m}'; od=work/f'out{m}'
        shutil.rmtree(md,ignore_errors=True); shutil.rmtree(od,ignore_errors=True)
        md.mkdir(); od.mkdir()
        b=f'{BASE}/mock{m}'
        ngc=md/'LRG_NGC.fits'; sgc=md/'LRG_SGC.fits'
        run(['curl','-fL','--retry','5',f'{b}/LRG_ffa_NGC_clustering.dat.fits','-o',str(ngc)])
        run(['curl','-fL','--retry','5',f'{b}/LRG_ffa_SGC_clustering.dat.fits','-o',str(sgc)])
        log=out/f'ezmock_{m}.log'
        with open(log,'w') as f:
            p=subprocess.run([sys.executable,str(finder),'--lrg-ngc',str(ngc),'--lrg-sgc',str(sgc),'--queries','50000','--mocks','5','--outdir',str(od)],stdout=f,stderr=subprocess.STDOUT)
        if p.returncode!=0:
            print(log.read_text()[-10000:], flush=True)
            raise SystemExit(p.returncode)
        v=pd.read_csv(od/'void_candidates.csv'); w=pd.read_csv(od/'wall_candidates.csv')
        row={'mock':m,'void_count':int(len(v)),'wall_count':int(len(w)),
             'max_void_score':float(v.score.max()) if len(v) else None,
             'max_wall_score':float(w.score.max()) if len(w) else None}
        for cid,t in TH_VOID.items(): row[f'n_void_score_ge_{cid}']=int((v.score>=t).sum()) if len(v) else 0
        for cid,t in TH_WALL.items(): row[f'n_wall_score_ge_{cid}']=int((w.score>=t).sum()) if len(w) else 0
        rows.append(row)
        print(json.dumps(row,indent=2),flush=True)
        shutil.rmtree(md,ignore_errors=True); shutil.rmtree(od,ignore_errors=True)
    df=pd.DataFrame(rows)
    df.to_csv(out/'ezmock_summary.csv',index=False)
    aggregate={'n_mocks':len(df),'mean_void_count':float(df.void_count.mean()),'mean_wall_count':float(df.wall_count.mean()),
               'fraction_mock_max_void_ge_V11':float((df.max_void_score>=TH_VOID['NEXO-V11']).mean()),
               'fraction_mock_max_wall_ge_W13':float((df.max_wall_score>=TH_WALL['NEXO-W13']).mean())}
    for cid in TH_VOID: aggregate[f'mean_count_ge_{cid}']=float(df[f'n_void_score_ge_{cid}'].mean())
    aggregate['mean_count_ge_NEXO-W13']=float(df['n_wall_score_ge_NEXO-W13'].mean())
    (out/'ezmock_aggregate.json').write_text(json.dumps(aggregate,indent=2))
    print('=== EZMOCK SUMMARY ==='); print(df.to_string(index=False)); print(json.dumps(aggregate,indent=2))
if __name__=='__main__': main()

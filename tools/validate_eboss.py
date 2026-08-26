#!/usr/bin/env python3
import argparse, math, json
from pathlib import Path
import numpy as np, pandas as pd
from astropy.io import fits
from scipy.integrate import cumulative_trapezoid
from scipy.spatial import cKDTree
H0=67.66; OM=.3111; C=299792.458; SEED=160826
zg=np.linspace(0,1.3,26001);dcg=(C/H0)*cumulative_trapezoid(1/np.sqrt(OM*(1+zg)**3+(1-OM)),zg,initial=0)
def dc(z):return np.interp(z,zg,dcg)
def xyz(ra,de,z):
 r=dc(np.asarray(z,float));a=np.deg2rad(ra);d=np.deg2rad(de);q=np.cos(d);return np.column_stack([r*q*np.cos(a),r*q*np.sin(a),r*np.sin(d)]).astype('f4')
def center(row):return xyz(np.array([row.ra_deg]),np.array([row.dec_deg]),np.array([row.z]))[0]
def load(paths):
 out=[]
 for reg,p in enumerate(paths):
  with fits.open(p,memmap=True) as H:
   d=H[1].data;n={x.upper():x for x in d.names};print(p,len(d),d.names,flush=True)
   ra=np.asarray(d[n['RA']],float);de=np.asarray(d[n['DEC']],float);z=np.asarray(d[n['Z']],float);m=np.isfinite(ra)&np.isfinite(de)&np.isfinite(z)&(z>=.6)&(z<=1.0)
   out.append((ra[m].astype('f4'),de[m].astype('f4'),z[m].astype('f4'),np.full(m.sum(),reg,'i1')))
   print(' kept',m.sum(),flush=True)
 return tuple(np.concatenate([x[k] for x in out]) for k in range(4))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--ngc',required=True);ap.add_argument('--sgc',required=True);ap.add_argument('--candidates',required=True);ap.add_argument('--mocks',type=int,default=200);ap.add_argument('--outdir',default='eboss');a=ap.parse_args();out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True)
 cand=pd.read_csv(a.candidates);ra,de,z,reg=load([a.ngc,a.sgc]);X=xyz(ra,de,z);T=cKDTree(X);rng=np.random.default_rng(SEED);CEN=np.vstack([center(r) for _,r in cand.iterrows()]);real=[]
 for i,r in cand.iterrows():real.append(T.query_ball_point(CEN[i],float(r.R_probe_mpc),return_length=True))
 real=np.asarray(real,float);null=np.empty((a.mocks,len(cand)))
 for m in range(a.mocks):
  zz=z.copy()
  for g in np.unique(reg):
   ii=np.where(reg==g)[0];zz[ii]=z[rng.permutation(ii)]
  TT=cKDTree(xyz(ra,de,zz))
  for i,r in cand.iterrows():null[m,i]=TT.query_ball_point(CEN[i],float(r.R_probe_mpc),return_length=True)
  if (m+1)%40==0:print('NULL',m+1,'/',a.mocks,flush=True)
 rows=[]
 for i,r in cand.iterrows():
  mu=null[:,i].mean();sd=null[:,i].std(ddof=1);kind=r.kind;lo=kind=='VOID';ext=((null[:,i]<=real[i]).sum() if lo else (null[:,i]>=real[i]).sum());p=(ext+1)/(a.mocks+1);ratio=real[i]/mu if mu>0 else np.nan;zs=(real[i]-mu)/sd if sd>0 else np.nan;support='ok' if mu>=8 else ('low_support' if mu>0 else 'outside_footprint')
  if support!='ok':confirm=False
  else:confirm=(ratio<.85 and zs<-1.8) if lo else (ratio>1.15 and zs>1.8)
  rows.append(dict(candidate_id=r.candidate_id,kind=kind,ra_deg=r.ra_deg,dec_deg=r.dec_deg,z=r.z,R_mpc=r.R_probe_mpc,real_count=int(real[i]),null_mean=float(mu),ratio=float(ratio) if np.isfinite(ratio) else np.nan,zscore=float(zs) if np.isfinite(zs) else np.nan,empirical_p=float(p),support=support,eboss_confirm=bool(confirm)))
 R=pd.DataFrame(rows);R.to_csv(out/'eboss_validation.csv',index=False);summary={'seed':SEED,'mocks':a.mocks,'eboss_objects':int(len(X)),'candidates':len(cand),'with_support':int((R.support=='ok').sum()),'confirmed':int(R.eboss_confirm.sum())};(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));print(R.to_string(index=False))
if __name__=='__main__':main()

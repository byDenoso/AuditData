#!/usr/bin/env python3
import argparse, math, json, gc
from pathlib import Path
import numpy as np, pandas as pd
from astropy.io import fits
from scipy.integrate import cumulative_trapezoid
from scipy.spatial import cKDTree
H0=67.66; OM=.3111; C=299792.458; h=H0/100.; SEED=8262026
RAD=np.array([50.,75.,100.,130.,140.,160.,200.])
zg=np.linspace(0,1.5,30001); dcg=(C/H0)*cumulative_trapezoid(1/np.sqrt(OM*(1+zg)**3+(1-OM)),zg,initial=0)
def dc(z):return np.interp(z,zg,dcg)
def xyz(ra,de,z):
 r=dc(np.asarray(z,float));a=np.deg2rad(ra);d=np.deg2rad(de);q=np.cos(d);return np.column_stack([r*q*np.cos(a),r*q*np.sin(a),r*np.sin(d)]).astype('f4')
def load(paths,zmin=.4,zmax=1.1):
 out=[]
 for reg,p in enumerate(paths):
  with fits.open(p,memmap=True) as H:
   d=H[1].data;n={x.upper():x for x in d.names};ra=np.asarray(d[n['RA']],float);de=np.asarray(d[n['DEC']],float);z=np.asarray(d[n['Z']],float);wt=np.asarray(d[n.get('WEIGHT','WEIGHT_COMP')],float) if ('WEIGHT' in n or 'WEIGHT_COMP' in n) else np.ones(len(d));m=np.isfinite(ra)&np.isfinite(de)&np.isfinite(z)&(z>=zmin)&(z<=zmax);out.append((ra[m].astype('f4'),de[m].astype('f4'),z[m].astype('f4'),wt[m].astype('f4'),np.full(m.sum(),reg,'i1')))
 return tuple(np.concatenate([x[k] for x in out]) for k in range(5))
def centers(c):return xyz(c.ra_deg.values,c.dec_deg.values,c.z.values)
def measure(tree,weights,CEN):
 n=np.zeros((len(CEN),len(RAD)));w=np.zeros_like(n)
 for j,r in enumerate(RAD):
  lists=tree.query_ball_point(CEN,r,workers=-1)
  for i,ids in enumerate(lists):n[i,j]=len(ids);w[i,j]=weights[np.asarray(ids,dtype=int)].sum() if len(ids) else 0
 return n,w
def validate(paths,cands,tracer,mocks,zmin,zmax,rng):
 ra,de,z,wt,reg=load(paths,zmin,zmax);X=xyz(ra,de,z);CEN=centers(cands);realT=cKDTree(X);rn,rw=measure(realT,wt,CEN);NN=np.empty((mocks,len(cands),len(RAD)));WW=np.empty_like(NN)
 for m in range(mocks):
  zz=z.copy()
  for g in np.unique(reg):
   ii=np.where(reg==g)[0];zz[ii]=z[rng.permutation(ii)]
  T=cKDTree(xyz(ra,de,zz));NN[m],WW[m]=measure(T,wt,CEN)
  if (m+1)%20==0:print(tracer,'NULL',m+1,'/',mocks,flush=True)
  del T
 mu=NN.mean(0);sd=NN.std(0,ddof=1);muw=WW.mean(0);sdw=WW.std(0,ddof=1)
 rows=[]
 for i,r in cands.iterrows():
  j=int(np.argmin(abs(RAD-float(r.R_probe_mpc))));kind=r.kind;lower=kind=='VOID';ex=(NN[:,i,j]<=rn[i,j]).sum() if lower else (NN[:,i,j]>=rn[i,j]).sum();exw=(WW[:,i,j]<=rw[i,j]).sum() if lower else (WW[:,i,j]>=rw[i,j]).sum();p=(ex+1)/(mocks+1);pw=(exw+1)/(mocks+1)
  rows.append(dict(candidate_id=r.candidate_id,tracer=tracer,R_test_mpc=float(RAD[j]),real_count=int(rn[i,j]),null_mean=float(mu[i,j]),ratio=float(rn[i,j]/max(mu[i,j],1)),zscore=float((rn[i,j]-mu[i,j])/max(sd[i,j],1e-9)),empirical_p=float(p),real_weight=float(rw[i,j]),null_weight_mean=float(muw[i,j]),weighted_ratio=float(rw[i,j]/max(muw[i,j],1e-9)),weighted_zscore=float((rw[i,j]-muw[i,j])/max(sdw[i,j],1e-9)),weighted_empirical_p=float(pw)))
 del realT,X,ra,de,z,wt,reg,NN,WW;gc.collect();return pd.DataFrame(rows)
def astra_crossmatch(cands,root):
 rows=[];allg=[]
 for p in Path(root).rglob('*.fits*'):
  try:
   with fits.open(p,memmap=False) as H:
    if len(H)<2 or H[1].data is None:continue
    d=H[1].data;n={x.upper():x for x in d.names}
    if not {'RA','DEC','REDSHIFT'}.issubset(n):continue
    ra=np.asarray(d[n['RA']],float);de=np.asarray(d[n['DEC']],float);z=np.asarray(d[n['REDSHIFT']],float);re=np.asarray(d[n['R_EFF']],float) if 'R_EFF' in n else np.full(len(d),np.nan);vid=np.asarray(d[n['VOID_ID']]) if 'VOID_ID' in n else np.arange(len(d));good=np.isfinite(ra)&np.isfinite(de)&np.isfinite(z)
    if good.any():allg.append((xyz(ra[good],de[good],z[good]),re[good],vid[good],np.full(good.sum(),str(p))))
  except Exception as e:print('skip',p,e,flush=True)
 if not allg:return pd.DataFrame()
 GX=np.vstack([x[0] for x in allg]);RE=np.concatenate([x[1] for x in allg]);VID=np.concatenate([x[2] for x in allg]);FN=np.concatenate([x[3] for x in allg]);T=cKDTree(GX);C=centers(cands);dist,idx=T.query(C,k=1)
 for i,r in cands.iterrows():
  k=idx[i];reff=float(RE[k]/h) if np.isfinite(RE[k]) else np.nan;ov=float(dist[i]/max(float(r.R_probe_mpc)+reff,1)) if np.isfinite(reff) else np.nan;rows.append(dict(candidate_id=r.candidate_id,nearest_astra_distance_mpc=float(dist[i]),astra_R_eff_mpc=reff,center_overlap_index=ov,astra_file=str(FN[k]),astra_void_id=str(VID[k]),astra_center_overlap=bool(np.isfinite(ov) and ov<1.0)))
 return pd.DataFrame(rows)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--lrg-ngc',required=True);ap.add_argument('--lrg-sgc',required=True);ap.add_argument('--elg-ngc',required=True);ap.add_argument('--elg-sgc',required=True);ap.add_argument('--candidates',required=True);ap.add_argument('--astra-root',required=True);ap.add_argument('--mocks',type=int,default=100);ap.add_argument('--outdir',default='validation');a=ap.parse_args();out=Path(a.outdir);out.mkdir(parents=True,exist_ok=True);cands=pd.read_csv(a.candidates);rng=np.random.default_rng(SEED)
 L=validate([a.lrg_ngc,a.lrg_sgc],cands,'LRG',a.mocks,.4,1.1,rng);E=validate([a.elg_ngc,a.elg_sgc],cands,'ELG',a.mocks,.8,1.1,rng);A=astra_crossmatch(cands,a.astra_root);V=L.merge(E,on='candidate_id',suffixes=('_lrg','_elg')).merge(A,on='candidate_id',how='left');V['robust_both_tracers']=((V.ratio_lrg<.85)&(V.ratio_elg<.85)&(V.weighted_ratio_lrg<.9)&(V.weighted_ratio_elg<.9))|((V.ratio_lrg>1.15)&(V.ratio_elg>1.15)&(V.weighted_ratio_lrg>1.1)&(V.weighted_ratio_elg>1.1));V.to_csv(out/'frozen_validation.csv',index=False);A.to_csv(out/'astra_crossmatch.csv',index=False);summary={'seed':SEED,'null_shuffles_per_tracer':a.mocks,'frozen_candidates':len(cands),'robust_both_tracers':int(V.robust_both_tracers.sum()),'astra_center_overlaps':int(A.astra_center_overlap.sum()) if len(A) else None};(out/'validation_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));print(V.to_string(index=False))
if __name__=='__main__':main()

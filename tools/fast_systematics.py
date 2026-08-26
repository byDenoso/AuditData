#!/usr/bin/env python3
import argparse,json,math
from pathlib import Path
import numpy as np,pandas as pd
from astropy.io import fits
from scipy.integrate import cumulative_trapezoid
from scipy.spatial import cKDTree
H0=67.66;OM=.3111;OL=1-OM;C=299792.458;RNG=np.random.default_rng(82620261335)
ZG=np.linspace(0,1.6,30000);DC=(C/H0)*cumulative_trapezoid(1/np.sqrt(OM*(1+ZG)**3+OL),ZG,initial=0)
def dc(z):return np.interp(z,ZG,DC)
def xyz(ra,de,z):
 r=dc(np.asarray(z,float));a=np.deg2rad(np.asarray(ra,float));d=np.deg2rad(np.asarray(de,float));q=np.cos(d);return np.c_[r*q*np.cos(a),r*q*np.sin(a),r*np.sin(d)]
def ctr(r):return xyz([r.ra_deg],[r.dec_deg],[r.z])[0]
def load(paths,zlo,zhi):
 A=[]
 for reg,p in enumerate(paths):
  with fits.open(p,memmap=True) as h:
   t=h[1].data;z=np.asarray(t['Z'],float);m=np.isfinite(z)&(z>=zlo)&(z<=zhi);names=set(t.names)
   A.append((np.asarray(t['RA'],float)[m],np.asarray(t['DEC'],float)[m],z[m],np.asarray(t['WEIGHT'],float)[m] if 'WEIGHT'in names else np.ones(m.sum()),np.full(m.sum(),reg,int)))
 return tuple(np.concatenate([x[i] for x in A]) for i in range(5))
def shuffle_z(z,reg,local=False):
 z2=z.copy()
 for g in np.unique(reg):
  ix=np.where(reg==g)[0]
  if not local:z2[ix]=z[ix][RNG.permutation(len(ix))]
  else:
   bins=np.floor(z[ix]/.05).astype(int)
   for b in np.unique(bins):
    j=ix[np.where(bins==b)[0]];z2[j]=z[j][RNG.permutation(len(j))]
 return z2
def query(tree,pts,R,w=None):
 out=[]
 for p in pts:
  ids=tree.query_ball_point(p,R)
  out.append(len(ids) if w is None else float(np.sum(w[ids])))
 return np.asarray(out,float)
def run_tracer(name,paths,zlo,cand,nbase=25,nlocal=15):
 ra,de,z,w,reg=load(paths,zlo,1.1);real=cKDTree(xyz(ra,de,z));rows=[]
 centers=np.array([ctr(r) for _,r in cand.iterrows()]);Rs=np.array([r.R_mpc for _,r in cand.iterrows()],float)
 # centers + six 25-Mpc perturbations per candidate
 allpts=[];owner=[]
 for i,c in enumerate(centers):
  allpts.append(c);owner.append((i,0))
  for ax in range(3):
   for s in (-1,1):q=c.copy();q[ax]+=25*s;allpts.append(q);owner.append((i,len([x for x in owner if x[0]==i])))
 allpts=np.asarray(allpts)
 # observed base, weighted, perturb
 obs=[];obsw=[];pert=[]
 for i,(c,R) in enumerate(zip(centers,Rs)):
  obs.append(query(real,[c],R)[0]);obsw.append(query(real,[c],R,w)[0]);pp=[allpts[k] for k,o in enumerate(owner) if o[0]==i];pert.append(query(real,pp,R))
 obs=np.array(obs);obsw=np.array(obsw)
 nb=[];nbw=[];npert=[[] for _ in range(len(cand))];nloc=[]
 for k in range(nbase):
  t=cKDTree(xyz(ra,de,shuffle_z(z,reg,False)));nb.append([query(t,[c],R)[0] for c,R in zip(centers,Rs)]);nbw.append([query(t,[c],R,w)[0] for c,R in zip(centers,Rs)])
  for i,R in enumerate(Rs):pp=[allpts[j] for j,o in enumerate(owner) if o[0]==i];npert[i].append(query(t,pp,R))
 for k in range(nlocal):
  t=cKDTree(xyz(ra,de,shuffle_z(z,reg,True)));nloc.append([query(t,[c],R)[0] for c,R in zip(centers,Rs)])
 nb=np.asarray(nb);nbw=np.asarray(nbw);nloc=np.asarray(nloc)
 for i,(_,r) in enumerate(cand.iterrows()):
  direction=-1 if r.kind=='VOID' else 1;mu=nb[:,i].mean();sd=max(nb[:,i].std(ddof=1),1);wm=nbw[:,i].mean();ws=max(nbw[:,i].std(ddof=1),1);lm=nloc[:,i].mean();ls=max(nloc[:,i].std(ddof=1),1)
  pr=pert[i]/np.maximum(np.mean(npert[i],axis=0),1)
  p=(1+(nb[:,i]<=obs[i]).sum())/(nbase+1) if direction<0 else (1+(nb[:,i]>=obs[i]).sum())/(nbase+1)
  rows.append(dict(candidate_id=r.candidate_id,tracer=name,ratio=obs[i]/mu,z=(obs[i]-mu)/sd,p=p,wratio=obsw[i]/wm,wz=(obsw[i]-wm)/ws,local_ratio=obs[i]/lm,local_z=(obs[i]-lm)/ls,center_med=float(np.median(pr)),center_worst=float(np.max(pr) if direction<0 else np.min(pr))))
 return pd.DataFrame(rows),real,(ra,de,z,w,reg)
def version_test(paths,cand):
 ra,de,z,w,reg=load(paths,.4,1.1);tr=cKDTree(xyz(ra,de,z));centers=np.array([ctr(r) for _,r in cand.iterrows()]);Rs=np.array(cand.R_mpc,float);obs=np.array([query(tr,[c],R)[0] for c,R in zip(centers,Rs)]);N=[]
 for k in range(15):t=cKDTree(xyz(ra,de,shuffle_z(z,reg)));N.append([query(t,[c],R)[0] for c,R in zip(centers,Rs)])
 N=np.asarray(N);return {cand.iloc[i].candidate_id:obs[i]/max(N[:,i].mean(),1) for i in range(len(cand))}
def scale_profile(tr,cand):
 out={}
 for _,r in cand.iterrows():
  c=ctr(r);out[r.candidate_id]={str(R):int(query(tr,[c],R)[0]) for R in [50,75,100,130,160,200]}
 return out
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--lrg',nargs=2,required=True);ap.add_argument('--elg',nargs=2,required=True);ap.add_argument('--lrg12',nargs=2,required=True);ap.add_argument('--cand',required=True);ap.add_argument('--out',default='fastbattery');a=ap.parse_args();o=Path(a.out);o.mkdir(exist_ok=True);cand=pd.read_csv(a.cand)
 L,trL,_=run_tracer('LRG',a.lrg,.4,cand);E,trE,_=run_tracer('ELG',a.elg,.8,cand);V=version_test(a.lrg12,cand);D=pd.concat([L,E]);D['v12_ratio']=D.candidate_id.map(V);D.to_csv(o/'tests.csv',index=False)
 prof=scale_profile(trL,cand);(o/'scale_profiles.json').write_text(json.dumps(prof,indent=2))
 verdict=[]
 for cid,g in D.groupby('candidate_id'):
  k=cand.set_index('candidate_id').loc[cid,'kind'];vals=[]
  for _,r in g.iterrows():
   if k=='VOID':vals += [r.ratio<.82,r.wratio<.84,r.local_ratio<.84,r.center_worst<.92]
   else:vals += [r.ratio>1.10,r.wratio>1.10,r.local_ratio>1.08,r.center_worst>1.05]
  v12=V[cid];vals.append(v12<.88 if k=='VOID' else v12>1.08);verdict.append((cid,'SURVIVES' if all(vals) else 'FAIL',sum(vals),len(vals)))
 S=pd.DataFrame(verdict,columns=['candidate_id','verdict','passed','total']);S.to_csv(o/'verdict.csv',index=False);print(D.to_string(index=False));print(S.to_string(index=False));print(json.dumps(prof,indent=2))
if __name__=='__main__':main()

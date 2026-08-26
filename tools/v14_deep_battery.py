#!/usr/bin/env python3
import argparse, json, math
from pathlib import Path
import numpy as np, pandas as pd
from astropy.io import fits
from scipy.integrate import cumulative_trapezoid
from scipy.spatial import cKDTree

H0=67.66; OM=.3111; OL=1-OM; C=299792.458
RA0=193.309380; DEC0=2.514639; Z0=.822175
RNG=np.random.default_rng(2608261401)
ZG=np.linspace(0,1.3,30000); DCG=(C/H0)*cumulative_trapezoid(1/np.sqrt(OM*(1+ZG)**3+OL),ZG,initial=0)
def dc(z): return np.interp(np.asarray(z,float),ZG,DCG)
def xyz(ra,de,z):
 r=dc(z); a=np.deg2rad(np.asarray(ra,float)); d=np.deg2rad(np.asarray(de,float)); q=np.cos(d)
 return np.c_[r*q*np.cos(a),r*q*np.sin(a),r*np.sin(d)]
def cxyz(z=Z0): return xyz([RA0],[DEC0],[z])[0]
def load(path,zlo=.65,zhi=1.0):
 with fits.open(path,memmap=True) as h:
  t=h[1].data; names={n.upper():n for n in t.names}; z=np.asarray(t[names['Z']],float)
  m=np.isfinite(z)&(z>=zlo)&(z<=zhi)
  ra=np.asarray(t[names['RA']],float)[m]; de=np.asarray(t[names['DEC']],float)[m]; z=z[m]
  w=np.asarray(t[names['WEIGHT']],float)[m] if 'WEIGHT' in names else np.ones(m.sum())
  ph=np.asarray(t[names['PHOTSYS']])[m] if 'PHOTSYS' in names else np.array(['?']*m.sum())
  return ra,de,z,w,ph

def qcount(tree,c,R,w=None):
 ids=tree.query_ball_point(c,R)
 return (float(len(ids)),float(np.sum(w[ids])) if w is not None else float(len(ids)))

def sphere_profiles(ra,de,z,w, rra,rde,rz,rw, label):
 X=xyz(ra,de,z); RX=xyz(rra,rde,rz); t=cKDTree(X); rt=cKDTree(RX); c=cxyz()
 # random normalization in local redshift slab
 dsel=(z>.75)&(z<.90); rsel=(rz>.75)&(rz<.90); alpha=dsel.sum()/max(rsel.sum(),1); alphaw=np.sum(w[dsel])/max(np.sum(rw[rsel]),1)
 radii=np.array([25,40,50,60,75,90,110,130,160,200],float)
 rows=[]; lastD=lastR=lastDw=lastRw=0
 for R in radii:
  D,Dw=qcount(t,c,R,w); RR,RRw=qcount(rt,c,R,rw); exp=alpha*RR; expw=alphaw*RRw
  rows.append(dict(tracer=label,R_mpc=R,data=D,random=RR,ratio=D/max(exp,1e-9),weighted_ratio=Dw/max(expw,1e-9),shell_ratio=(D-lastD)/max(alpha*(RR-lastR),1e-9),shell_weighted_ratio=(Dw-lastDw)/max(alphaw*(RRw-lastRw),1e-9)))
  lastD,lastR,lastDw,lastRw=D,RR,Dw,RRw
 return pd.DataFrame(rows),t,rt,alpha,alphaw,X,RX

def z_sweep(tree,rt,alpha,w,rw):
 rows=[]
 for zc in np.arange(.76,.886,.005):
  c=cxyz(zc); D,Dw=qcount(tree,c,75,w); RR,RRw=qcount(rt,c,75,rw)
  rows.append(dict(z=zc,data=D,random=RR,ratio=D/max(alpha*RR,1e-9),weighted_ratio=Dw/max(alpha*RRw,1e-9)))
 return pd.DataFrame(rows)

def center_grid(tree,rt,alpha,w,rw):
 c=cxyz(); rows=[]
 for dx in [-50,-25,0,25,50]:
  for dy in [-50,-25,0,25,50]:
   for dz in [-50,-25,0,25,50]:
    q=c+np.array([dx,dy,dz]); D,Dw=qcount(tree,q,75,w); RR,RRw=qcount(rt,q,75,rw)
    rows.append(dict(dx=dx,dy=dy,dz=dz,ratio=D/max(alpha*RR,1e-9),data=D,random=RR))
 return pd.DataFrame(rows)

def angular_mask_test(rra,rde,rz):
 # use random points only: compare random angular support at V14 to valid same-dec centers
 m=(rz>Z0-.0125)&(rz<Z0+.0125); ra=rra[m]; de=rde[m]
 U=xyz(ra,de,np.ones(len(ra))); U/=np.linalg.norm(U,axis=1)[:,None]; tr=cKDTree(U)
 theta=75/dc(Z0); chord=2*np.sin(theta/2); u0=xyz([RA0],[DEC0],[1.0])[0]; u0/=np.linalg.norm(u0)
 obs=tr.query_ball_point(u0,chord,return_length=True)
 if len(U)>1500: idx=RNG.choice(len(U),1500,replace=False)
 else: idx=np.arange(len(U))
 vals=tr.query_ball_point(U[idx],chord,return_length=True)
 return dict(random_slice_objects=int(len(U)),theta_deg=float(np.degrees(theta)),candidate_random_count=int(obs),control_median=float(np.median(vals)),control_p05=float(np.quantile(vals,.05)),control_p95=float(np.quantile(vals,.95)),candidate_to_median=float(obs/max(np.median(vals),1)),low_tail_empirical_p=float((np.sum(vals<=obs)+1)/(len(vals)+1)))

def fine_slice_angular(ra,de,z,label):
 m=(z>Z0-.0125)&(z<Z0+.0125); ra=ra[m];de=de[m]
 U=xyz(ra,de,np.ones(len(ra))); U/=np.linalg.norm(U,axis=1)[:,None]; tr=cKDTree(U)
 theta=75/dc(Z0); chord=2*np.sin(theta/2); u0=xyz([RA0],[DEC0],[1.0])[0];u0/=np.linalg.norm(u0)
 obs=tr.query_ball_point(u0,chord,return_length=True)
 idx=RNG.choice(len(U),min(5000,len(U)),replace=False); vals=tr.query_ball_point(U[idx],chord,return_length=True)
 return dict(tracer=label,objects=int(len(U)),obs=int(obs),median=float(np.median(vals)),ratio=float(obs/max(np.median(vals),1)),empirical_p=float((np.sum(vals<=obs)+1)/(len(vals)+1)),q05=float(np.quantile(vals,.05)),q95=float(np.quantile(vals,.95)))

def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--lrg',required=True);ap.add_argument('--elg',required=True);ap.add_argument('--random',required=True);ap.add_argument('--out',default='v14out');a=ap.parse_args();o=Path(a.out);o.mkdir(parents=True,exist_ok=True)
 lra,lde,lz,lw,lph=load(a.lrg); era,ede,ez,ew,eph=load(a.elg,.8,1.0); rra,rde,rz,rw,rph=load(a.random)
 L,lt,rt,alpha,alphaw,LX,RX=sphere_profiles(lra,lde,lz,lw,rra,rde,rz,rw,'LRG')
 # random catalog is LRG-selected; ELG profile uses structured same-z sky rather than LRG random normalization
 Erows=[]; et=cKDTree(xyz(era,ede,ez)); c=cxyz();
 for R in [25,40,50,60,75,90,110,130,160,200]:
  D,Dw=qcount(et,c,R,ew); Erows.append(dict(tracer='ELG',R_mpc=R,data=D,weighted_data=Dw))
 E=pd.DataFrame(Erows)
 L.to_csv(o/'lrg_profiles.csv',index=False);E.to_csv(o/'elg_profiles.csv',index=False)
 zsw=z_sweep(lt,rt,alpha,lw,rw);zsw.to_csv(o/'redshift_sweep.csv',index=False)
 cg=center_grid(lt,rt,alpha,lw,rw);cg.to_csv(o/'center_grid.csv',index=False)
 mask=angular_mask_test(rra,rde,rz); fine=[fine_slice_angular(lra,lde,lz,'LRG'),fine_slice_angular(era,ede,ez,'ELG')]
 # compensation metric: strongest shell beyond 75, and integrated return toward unity
 outer=L[L.R_mpc>=90]; comp=float(outer.shell_ratio.max()); rcomp=float(outer.loc[outer.shell_ratio.idxmax(),'R_mpc'])
 basin=float((cg.ratio<.8).mean()); minrow=cg.loc[cg.ratio.idxmin()].to_dict()
 summary={'candidate':'NEXO-V14','ra':RA0,'dec':DEC0,'z':Z0,'R_nominal_mpc':75,'lrg_objects':len(lz),'elg_objects':len(ez),'random_objects':len(rz),'random_alpha':alpha,'mask_test':mask,'fine_slice_angular':fine,'compensation_shell_max_ratio':comp,'compensation_shell_R_mpc':rcomp,'center_grid_fraction_ratio_lt_0p8':basin,'center_grid_min':minrow,'lrg_ratio_75':float(L[L.R_mpc==75].ratio.iloc[0]),'lrg_weighted_ratio_75':float(L[L.R_mpc==75].weighted_ratio.iloc[0])}
 (o/'summary.json').write_text(json.dumps(summary,indent=2,default=float));print(json.dumps(summary,indent=2,default=float));print('LRG PROFILE');print(L.to_string(index=False));print('FINE SLICE',fine);print('MASK',mask)
if __name__=='__main__': main()

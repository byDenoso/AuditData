#!/usr/bin/env python3
import argparse,json,math
from pathlib import Path
import numpy as np,pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize_scalar
C=299792.458; OM=.33; HREF=67.5; RNG=np.random.default_rng(2608262016)

def E(z): return np.sqrt(OM*(1+z)**3+1-OM)
def chi_of_z(z,H=HREF):
 z=np.asarray(z,float);zg=np.linspace(0,max(.3,float(np.max(z))*1.02),12000);ch=(C/H)*cumulative_trapezoid(1/E(zg),zg,initial=0);return np.interp(z,zg,ch)
def mu_model(z,H):
 z=np.asarray(z,float); d=(1+z)*chi_of_z(z,H); return 5*np.log10(np.maximum(d,1e-12))+25
def unitvec(ra,dec):
 a=np.deg2rad(np.asarray(ra,float));d=np.deg2rad(np.asarray(dec,float));q=np.cos(d);return np.c_[q*np.cos(a),q*np.sin(a),np.sin(d)]
def fit_H(z,mu,sig):
 def f(H):
  r=mu-mu_model(z,H);return float(np.sum((r/np.maximum(sig,.04))**2))
 r=minimize_scalar(f,bounds=(55,85),method='bounded');return float(r.x),float(r.fun)
def fib_axes(n=500):
 i=np.arange(n);phi=(1+5**.5)/2;th=2*np.pi*i/phi;z=1-2*(i+.5)/n;r=np.sqrt(1-z*z);return np.c_[r*np.cos(th),r*np.sin(th),z]
def hemisphere_scan(df,zcol,zmax):
 m=(df.IS_CALIBRATOR==0)&(df[zcol]>=.01)&(df[zcol]<zmax)&np.isfinite(df.MU_SH0ES)&np.isfinite(df.RA)&np.isfinite(df.DEC)
 q=df.loc[m];z=q[zcol].to_numpy(float);mu=q.MU_SH0ES.to_numpy(float);sig=q.MU_SH0ES_ERR_DIAG.to_numpy(float);U=unitvec(q.RA,q.DEC)
 best=None
 for ax in fib_axes():
  s=U@ax; a=s>=0;b=~a
  if a.sum()<20 or b.sum()<20:continue
  h1,_=fit_H(z[a],mu[a],sig[a]);h2,_=fit_H(z[b],mu[b],sig[b]);dh=abs(h1-h2)
  if best is None or dh>best[0]:best=(dh,h1,h2,ax,a.sum(),b.sum())
 dh,h1,h2,ax,n1,n2=best;ra=np.degrees(np.arctan2(ax[1],ax[0]))%360;dec=np.degrees(np.arcsin(ax[2]))
 return dict(z_definition=zcol,zmax=zmax,N=len(q),H_plus=h1,H_minus=h2,deltaH_abs=dh,deltaH_percent=100*dh/((h1+h2)/2),axis_ra_deg=float(ra),axis_dec_deg=float(dec),N_plus=int(n1),N_minus=int(n2))
def zmin_scan(df,zcol):
 out=[]
 for zmin in [.005,.01,.015,.023,.03,.05]:
  m=(df.IS_CALIBRATOR==0)&(df[zcol]>=zmin)&(df[zcol]<.15)&np.isfinite(df.MU_SH0ES)
  q=df.loc[m];H,ch=fit_H(q[zcol].to_numpy(float),q.MU_SH0ES.to_numpy(float),q.MU_SH0ES_ERR_DIAG.to_numpy(float));out.append(dict(z_definition=zcol,zmin=zmin,zmax=.15,N=len(q),H0=H,chi2=ch))
 return out

def synth_bias(df,R,delta,offset_frac,norient=120):
 # noiseless synthetic Hubble flow seen by an observer embedded in a spherical top-hat underdensity
 m=(df.IS_CALIBRATOR==0)&(df.zHD>=.023)&(df.zHD<.15)&np.isfinite(df.RA)&np.isfinite(df.DEC)
 q=df.loc[m]; zcos=q.zHD.to_numpy(float);U=unitvec(q.RA,q.DEC);r=chi_of_z(zcos,HREF);X=U*r[:,None]
 f=OM**.55; alpha=-(f*HREF*delta/3.0) # km/s/Mpc; outward for delta<0
 vals=[]
 for _ in range(norient):
  v=RNG.normal(size=3);v/=np.linalg.norm(v);obs=v*(offset_frac*R);center=np.zeros(3)
  robs=obs-center; vobs=alpha*robs if np.linalg.norm(robs)<R else np.zeros(3)
  rel=X+obs-center; rr=np.linalg.norm(rel,axis=1);vs=np.zeros_like(rel);inside=rr<R;vs[inside]=alpha*rel[inside]
  vlos=np.sum((vs-vobs)*U,axis=1)
  zapp=zcos+(1+zcos)*vlos/C
  mu=mu_model(zcos,HREF) # true luminosity distances
  Hfit,_=fit_H(zapp,mu,np.full(len(mu),.08));vals.append(100*(Hfit/HREF-1))
 return dict(R_Mpc=R,delta_density=delta,offset_fraction=offset_frac,N=len(q),bias_mean_percent=float(np.mean(vals)),bias_sd_percent=float(np.std(vals,ddof=1)),bias_p05=float(np.quantile(vals,.05)),bias_p95=float(np.quantile(vals,.95)),bias_max=float(np.max(vals)))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--data',required=True);ap.add_argument('--out',default='observerh0');a=ap.parse_args();o=Path(a.out);o.mkdir(parents=True,exist_ok=True);d=pd.read_csv(a.data,sep=r'\s+')
 zs=[]
 for zc in ['zHD','zCMB']:zs+=zmin_scan(d,zc)
 pd.DataFrame(zs).to_csv(o/'zmin_H0.csv',index=False)
 hs=[]
 for zc in ['zHD','zCMB']:
  for zm in [.03,.05,.08,.15]:hs.append(hemisphere_scan(d,zc,zm))
 pd.DataFrame(hs).to_csv(o/'hemisphere_H0.csv',index=False)
 sims=[]
 for R in [40,70,100,150,300]:
  for de in [-.1,-.2,-.3,-.4]:
   for off in [0,.25,.5,.75]:sims.append(synth_bias(d,R,de,off))
 S=pd.DataFrame(sims);S.to_csv(o/'observer_void_grid.csv',index=False)
 # empirical perspective shift directly from same SNe zCMB -> zHD
 m=(d.IS_CALIBRATOR==0)&(d.zHD>=.023)&(d.zHD<.15)&np.isfinite(d.MU_SH0ES)
 q=d.loc[m]; hh,_=fit_H(q.zHD.to_numpy(float),q.MU_SH0ES.to_numpy(float),q.MU_SH0ES_ERR_DIAG.to_numpy(float));hc,_=fit_H(q.zCMB.to_numpy(float),q.MU_SH0ES.to_numpy(float),q.MU_SH0ES_ERR_DIAG.to_numpy(float))
 direct=dict(N=len(q),H0_zHD=hh,H0_zCMB=hc,delta_kmsmpc=hh-hc,delta_percent=100*(hh/hc-1))
 # requirements for target boosts in linear centered top-hat approx delta=-3 dH/(f H)
 f=OM**.55; req=[]
 for pct in [1,3,5,8]:req.append(dict(target_H_bias_percent=pct,required_delta_density=-3*(pct/100)/f))
 summary={'direct_redshift_correction_effect':direct,'max_simulated_bias':S.loc[S.bias_mean_percent.idxmax()].to_dict(),'target_density_requirements':req}
 (o/'summary.json').write_text(json.dumps(summary,indent=2,default=float));pd.DataFrame(req).to_csv(o/'required_density.csv',index=False);print(json.dumps(summary,indent=2,default=float));print(pd.DataFrame(zs).to_string(index=False));print(pd.DataFrame(hs).to_string(index=False));print(S.sort_values('bias_mean_percent',ascending=False).head(15).to_string(index=False))
if __name__=='__main__':main()
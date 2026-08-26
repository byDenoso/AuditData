#!/usr/bin/env python3
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.integrate import cumulative_trapezoid
C=299792.458; OM=.33; RNG=np.random.default_rng(2608262117)
PEER_BACKGROUNDS=[70.391,71.62]
TARGETS=[72.0,72.5,73.0,73.5]
RGRID=[70,100,150,200,250,300]
DGRID=[-.10,-.15,-.20,-.25,-.30,-.35,-.40]
OFFGRID=[0,.25,.50,.75]

def E(z): return np.sqrt(OM*(1+z)**3+1-OM)
def chi_of_z(z,H):
 z=np.asarray(z,float);zg=np.linspace(0,max(.3,float(np.max(z))*1.02),12000);ch=(C/H)*cumulative_trapezoid(1/E(zg),zg,initial=0);return np.interp(z,zg,ch)
def mu_model(z,H):
 z=np.asarray(z,float);d=(1+z)*chi_of_z(z,H);return 5*np.log10(np.maximum(d,1e-12))+25
def unitvec(ra,dec):
 a=np.deg2rad(np.asarray(ra,float));d=np.deg2rad(np.asarray(dec,float));q=np.cos(d);return np.c_[q*np.cos(a),q*np.sin(a),np.sin(d)]
def fit_H_exact(z,mu,sig,Href=73.0):
 # Exact because, at fixed Omega_m, d_L(z,H)=d_L(z,Href)*Href/H.
 ref=mu_model(z,Href);w=1/np.maximum(np.asarray(sig,float),.04)**2;delta=np.asarray(mu,float)-ref
 mean_delta=float(np.sum(w*delta)/np.sum(w));H=float(Href*10**(-mean_delta/5))
 resid=delta+5*np.log10(H/Href);chi=float(np.sum(w*resid*resid));return H,chi
def synth_bias(df,Hbg,R,delta,offset_frac,norient=160):
 m=(df.IS_CALIBRATOR==0)&(df.zHD>=.023)&(df.zHD<.15)&np.isfinite(df.RA)&np.isfinite(df.DEC)
 q=df.loc[m];zcos=q.zHD.to_numpy(float);U=unitvec(q.RA,q.DEC);r=chi_of_z(zcos,Hbg);X=U*r[:,None];mu=mu_model(zcos,Hbg);sig=np.full(len(mu),.08)
 f=OM**.55;alpha=-(f*Hbg*delta/3.0);vals=[]
 for _ in range(norient):
  axis=RNG.normal(size=3);axis/=np.linalg.norm(axis);obs=axis*(offset_frac*R)
  vobs=alpha*obs if np.linalg.norm(obs)<R else np.zeros(3)
  rel=X+obs;rr=np.linalg.norm(rel,axis=1);vs=np.zeros_like(rel);inside=rr<R;vs[inside]=alpha*rel[inside]
  vlos=np.sum((vs-vobs)*U,axis=1);zapp=zcos+(1+zcos)*vlos/C
  Hfit,_=fit_H_exact(zapp,mu,sig,Href=Hbg);vals.append(100*(Hfit/Hbg-1))
 vals=np.asarray(vals)
 return dict(H0_background=Hbg,R_Mpc=R,delta_density=delta,offset_fraction=offset_frac,N=len(q),bias_mean_percent=float(vals.mean()),bias_sd_percent=float(vals.std(ddof=1)),bias_p05=float(np.quantile(vals,.05)),bias_p50=float(np.quantile(vals,.50)),bias_p95=float(np.quantile(vals,.95)),bias_min=float(vals.min()),bias_max=float(vals.max()),H0_app_mean=float(Hbg*(1+vals.mean()/100)),H0_app_p05=float(Hbg*(1+np.quantile(vals,.05)/100)),H0_app_p95=float(Hbg*(1+np.quantile(vals,.95)/100)))
def direct_correction(df):
 m=(df.IS_CALIBRATOR==0)&(df.zHD>=.023)&(df.zHD<.15)&np.isfinite(df.MU_SH0ES);q=df.loc[m]
 hhd,_=fit_H_exact(q.zHD.to_numpy(float),q.MU_SH0ES.to_numpy(float),q.MU_SH0ES_ERR_DIAG.to_numpy(float),73.0)
 hcmb,_=fit_H_exact(q.zCMB.to_numpy(float),q.MU_SH0ES.to_numpy(float),q.MU_SH0ES_ERR_DIAG.to_numpy(float),73.0)
 return dict(N=len(q),H0_zHD=hhd,H0_zCMB=hcmb,delta_kmsmpc=hhd-hcmb,delta_percent=100*(hhd/hcmb-1))
def target_rows(S):
 out=[];f=OM**.55
 for Hbg in PEER_BACKGROUNDS:
  sub=S[S.H0_background==Hbg].copy()
  for target in TARGETS:
   need=100*(target/Hbg-1);centered_delta=-3*(need/100)/f
   for off in OFFGRID:
    s=sub[sub.offset_fraction==off].copy();s['abs_target_error']=np.abs(s.H0_app_mean-target);best=s.sort_values('abs_target_error').iloc[0]
    out.append(dict(H0_background=Hbg,target_H0=target,required_boost_percent=need,linear_all_inside_required_delta=centered_delta,offset_fraction=off,best_R_Mpc=float(best.R_Mpc),best_delta_density=float(best.delta_density),best_H0_app_mean=float(best.H0_app_mean),best_H0_app_p05=float(best.H0_app_p05),best_H0_app_p95=float(best.H0_app_p95),target_error=float(best.abs_target_error)))
 return pd.DataFrame(out)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--data',required=True);ap.add_argument('--out',default='peerh0');a=ap.parse_args();o=Path(a.out);o.mkdir(parents=True,exist_ok=True);d=pd.read_csv(a.data,sep=r'\s+')
 rows=[]
 for Hbg in PEER_BACKGROUNDS:
  for R in RGRID:
   for de in DGRID:
    for off in OFFGRID:rows.append(synth_bias(d,Hbg,R,de,off))
 S=pd.DataFrame(rows);S.to_csv(o/'peer_observer_void_grid.csv',index=False);T=target_rows(S);T.to_csv(o/'peer_target_matches.csv',index=False);direct=direct_correction(d)
 H=PEER_BACKGROUNDS[0];need=100*(73/H-1);exact=T[(T.H0_background==H)&(T.target_H0==73)].sort_values('offset_fraction');fam=S[(S.H0_background==H)&(S.delta_density==-0.20)&(S.R_Mpc.isin([150,200,250,300]))].copy()
 summary={'peer_anchor_free_H0':H,'peer_local_bracket_H0':PEER_BACKGROUNDS[1],'target_H0':73.0,'required_anchor_free_boost_percent':need,'direct_Pantheon_redshift_correction_effect':direct,'best_grid_matches_anchor_free_to_73':exact.to_dict(orient='records'),'delta_minus_0p20_family':fam.to_dict(orient='records')}
 (o/'summary.json').write_text(json.dumps(summary,indent=2,default=float));print(json.dumps(summary,indent=2,default=float));print('\nTARGET MATCHES\n',T.to_string(index=False));print('\nTOP ANCHOR-FREE TO 73\n',S[S.H0_background==H].assign(err=lambda x:abs(x.H0_app_mean-73)).sort_values('err').head(25).to_string(index=False))
if __name__=='__main__':main()

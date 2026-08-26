#!/usr/bin/env python3
import argparse,json,math
from pathlib import Path
import numpy as np
from astropy.io import fits
from scipy.spatial import cKDTree
RA0=193.309380;DEC0=2.514639;Z0=.822175
H0=67.66;OM=.3111;C=299792.458
RNG=np.random.default_rng(2608261417)

def dc(z):
 zz=np.linspace(0,z,10000);e=np.sqrt(OM*(1+zz)**3+(1-OM));return (C/H0)*np.trapezoid(1/e,zz)
def uv(ra,de):
 a=np.deg2rad(np.asarray(ra,float));d=np.deg2rad(np.asarray(de,float));q=np.cos(d);return np.c_[q*np.cos(a),q*np.sin(a),np.sin(d)]
def load(p):
 with fits.open(p,memmap=True) as h:
  t=h[1].data;z=np.asarray(t['Z'],float);m=np.isfinite(z)&(z>Z0-.0125)&(z<Z0+.0125)
  return np.asarray(t['RA'],float)[m],np.asarray(t['DEC'],float)[m],z[m]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--lrg',required=True);ap.add_argument('--random',required=True);ap.add_argument('--out',default='supportout');a=ap.parse_args();o=Path(a.out);o.mkdir(parents=True,exist_ok=True)
 ra,de,z=load(a.lrg);rra,rde,rz=load(a.random);D=uv(ra,de);R=uv(rra,rde);dt=cKDTree(D);rt=cKDTree(R);u0=uv([RA0],[DEC0])[0]
 theta=75/dc(Z0);ch=2*np.sin(theta/2);obs=int(dt.query_ball_point(u0,ch,return_length=True));rs0=int(rt.query_ball_point(u0,ch,return_length=True))
 n=min(20000,len(R));idx=RNG.choice(len(R),n,replace=False);cent=R[idx];sep=np.linalg.norm(cent-u0,axis=1);cent=cent[sep>2*ch]
 rs=np.asarray(rt.query_ball_point(cent,ch,return_length=True));ds=np.asarray(dt.query_ball_point(cent,ch,return_length=True));
 lo=.85*rs0;hi=1.15*rs0;keep=(rs>=lo)&(rs<=hi);vals=ds[keep];supports=rs[keep]
 # tighter 10% diagnostic too
 k10=(rs>=.9*rs0)&(rs<=1.1*rs0);v10=ds[k10]
 def stat(v):
  return {'n':int(len(v)),'mean':float(np.mean(v)),'median':float(np.median(v)),'sd':float(np.std(v,ddof=1)),'q05':float(np.quantile(v,.05)),'q95':float(np.quantile(v,.95)),'zscore':float((obs-np.mean(v))/np.std(v,ddof=1)),'lower_tail_p':float((np.sum(v<=obs)+1)/(len(v)+1))} if len(v)>2 else {'n':int(len(v))}
 out={'candidate':'NEXO-V14','z_slice':[Z0-.0125,Z0+.0125],'aperture_deg':float(np.degrees(theta)),'obs_lrg':obs,'candidate_random_support':rs0,'all_sample_n':int(len(ds)),'support_15pct':stat(vals),'support_15pct_mean_random_support':float(np.mean(supports)) if len(supports) else None,'support_10pct':stat(v10)}
 (o/'support_matched.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':main()

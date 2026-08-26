#!/usr/bin/env python3
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize

C=299792.458
RNG=np.random.default_rng(2608261901)
H0REF=73.0; OM=.33
SHELLS=[(.005,.01),(.01,.02),(.02,.03),(.03,.05),(.05,.08),(.08,.12),(.12,.2)]


def E(z,om=OM): return np.sqrt(om*(1+z)**3+(1-om))
def mu_model(z,H0=H0REF,om=OM):
    z=np.asarray(z,float);zg=np.linspace(0,max(.25,float(z.max())*1.02),10000)
    chi=(C/H0)*cumulative_trapezoid(1/E(zg,om),zg,initial=0)
    dl=(1+z)*np.interp(z,zg,chi)
    return 5*np.log10(dl)+25

def unitvec(ra,dec):
    a=np.deg2rad(np.asarray(ra,float));d=np.deg2rad(np.asarray(dec,float));c=np.cos(d)
    return np.c_[c*np.cos(a),c*np.sin(a),np.sin(d)]

def vec_to_radec(v):
    v=np.asarray(v,float);v=v/np.linalg.norm(v)
    ra=np.degrees(np.arctan2(v[1],v[0]))%360
    dec=np.degrees(np.arcsin(v[2]))
    return float(ra),float(dec)

def wls(X,y,sig):
    w=1/np.maximum(sig,0.03)**2
    XtW=X.T*w
    M=XtW@X
    b=XtW@y
    beta=np.linalg.solve(M,b)
    cov=np.linalg.inv(M)
    r=y-X@beta
    chi=float(np.sum(w*r*r))
    return beta,cov,chi

def fit_shell(df,zcol,lo,hi,nperm=2000):
    m=(df[zcol]>=lo)&(df[zcol]<hi)&(df.IS_CALIBRATOR==0)&np.isfinite(df.MU_SH0ES)&np.isfinite(df.RA)&np.isfinite(df.DEC)
    q=df.loc[m].copy()
    if len(q)<25:return None
    z=q[zcol].to_numpy(float); mu=q.MU_SH0ES.to_numpy(float); sig=q.MU_SH0ES_ERR_DIAG.to_numpy(float)
    U=unitvec(q.RA,q.DEC)
    base=mu_model(z,H0REF,OM)
    resid=mu-base
    # monopole + dipole in distance modulus
    X=np.c_[np.ones(len(q)),U]
    beta,cov,chi=wls(X,resid,sig)
    dip=beta[1:]; amp=float(np.linalg.norm(dip)); amp_err=float(np.sqrt(np.trace(cov[1:,1:])))
    ra,dec=vec_to_radec(dip) if amp>0 else (np.nan,np.nan)
    # quadrupole five independent basis + dipole + monopole
    x,y,zv=U[:,0],U[:,1],U[:,2]
    Q=np.c_[x*x-y*y,2*zv*zv-x*x-y*y,2*x*y,2*x*zv,2*y*zv]
    X2=np.c_[np.ones(len(q)),U,Q]
    b2,c2,chi2=wls(X2,resid,sig)
    dchi2_quad=float(chi-chi2)
    # permutation null for dipole amplitude, sky fixed and residuals shuffled
    amps=[];dqs=[]
    for _ in range(nperm):
        yp=RNG.permutation(resid)
        bp,_,chip=wls(X,yp,sig);amps.append(float(np.linalg.norm(bp[1:])))
        bq,_,chiq=wls(X2,yp,sig);dqs.append(float(chip-chiq))
    amps=np.asarray(amps);dqs=np.asarray(dqs)
    p=(np.sum(amps>=amp)+1)/(len(amps)+1)
    pq=(np.sum(dqs>=dchi2_quad)+1)/(len(dqs)+1)
    # convert dipole magnitude in mu to fractional H dipole approximately: dH/H=-ln(10)/5 dmu
    fracH=math.log(10)/5*amp
    return dict(z_definition=zcol,zlo=lo,zhi=hi,zmid=.5*(lo+hi),N=int(len(q)),monopole_mag=float(beta[0]),dipole_mag=amp,dipole_ra_deg=ra,dipole_dec_deg=dec,dipole_H_fraction=fracH,dipole_H_percent=100*fracH,dipole_empirical_p=float(p),quadrupole_delta_chi2=dchi2_quad,quadrupole_empirical_p=float(pq),chi2_dipole=chi,chi2_dipole_quadrupole=chi2)

def fit_vpec_alignment(df,zcol):
    m=(df.IS_CALIBRATOR==0)&(df[zcol]>.005)&(df[zcol]<.08)&np.isfinite(df.VPEC)&np.isfinite(df.RA)&np.isfinite(df.DEC)&np.isfinite(df.MU_SH0ES)
    q=df.loc[m].copy();U=unitvec(q.RA,q.DEC);v=q.VPEC.to_numpy(float);sig=np.full(len(q),250.)
    X=np.c_[np.ones(len(q)),U]
    beta,cov,chi=wls(X,v,sig);vec=beta[1:];amp=float(np.linalg.norm(vec));ra,dec=vec_to_radec(vec)
    amps=[]
    for _ in range(3000):
        bp,_,_=wls(X,RNG.permutation(v),sig);amps.append(float(np.linalg.norm(bp[1:])))
    amps=np.asarray(amps);p=(np.sum(amps>=amp)+1)/(len(amps)+1)
    return dict(z_definition=zcol,N=int(len(q)),vpec_dipole_kms=amp,ra_deg=ra,dec_deg=dec,empirical_p=float(p),monopole_kms=float(beta[0]))

def angular_sep(ra1,dec1,ra2,dec2):
    u=unitvec([ra1],[dec1])[0];v=unitvec([ra2],[dec2])[0]
    return float(np.degrees(np.arccos(np.clip(u@v,-1,1))))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data',required=True);ap.add_argument('--out',default='anisotropy');a=ap.parse_args();o=Path(a.out);o.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(a.data,sep=r'\s+')
    rows=[]
    for zcol in ['zHD','zCMB']:
        for lo,hi in SHELLS:
            r=fit_shell(d,zcol,lo,hi)
            if r: rows.append(r)
    R=pd.DataFrame(rows);R.to_csv(o/'anisotropy_shells.csv',index=False)
    vp=[fit_vpec_alignment(d,zcol) for zcol in ['zHD','zCMB']]
    pd.DataFrame(vp).to_csv(o/'vpec_dipole.csv',index=False)
    # reference directions: Virgo/M87 J2000, approximate Local Void direction from Galactic center-side obscured region represented by l~20,b~0 converted approximately to ICRS offline constant
    refs={'Virgo_M87':(187.706,12.391),'LocalVoid_proxy':(277.0,-14.0),'CMB_dipole_apex':(167.94,-6.94)}
    comp=[]
    for _,r in R.iterrows():
        for name,(ra,dec) in refs.items(): comp.append(dict(z_definition=r.z_definition,zlo=r.zlo,zhi=r.zhi,reference=name,separation_deg=angular_sep(r.dipole_ra_deg,r.dipole_dec_deg,ra,dec)))
    pd.DataFrame(comp).to_csv(o/'axis_comparisons.csv',index=False)
    # summary select low-z strongest per definition
    summary={}
    for zcol in ['zHD','zCMB']:
        s=R[R.z_definition==zcol]
        if len(s):
            best=s.loc[s.dipole_empirical_p.idxmin()].to_dict();summary[zcol]={'strongest_dipole':best,'vpec_dipole':next(x for x in vp if x['z_definition']==zcol)}
    (o/'summary.json').write_text(json.dumps(summary,indent=2,default=float))
    print(json.dumps(summary,indent=2,default=float));print(R.to_string(index=False));print(pd.DataFrame(vp).to_string(index=False))
if __name__=='__main__':main()

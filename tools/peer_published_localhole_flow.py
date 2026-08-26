#!/usr/bin/env python3
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid

C=299792.458
OM=0.33
RNG=np.random.default_rng(2608262217)
HBACK=[70.391,71.62]
BGAL=[1.0,1.1,1.2]
OFFSETS=[0.,25.,50.,75.,100.]
# Published Local Hole envelope: ~20% at 100 h^-1 Mpc, ~13% cumulative at 150 h^-1 Mpc.
PROFILES=[
 ('low',-.18,-.12),
 ('nominal',-.20,-.13),
 ('high',-.22,-.14),
]
BOUNDARIES=[('uncompensated',None),('comp1p5',1.5),('comp2p0',2.0)]
NORI=320

def E(z): return np.sqrt(OM*(1+z)**3 + 1-OM)
def chi_of_z(z,H):
    z=np.asarray(z,float); zg=np.linspace(0,max(.25,float(np.max(z))*1.02),25000)
    ch=(C/H)*cumulative_trapezoid(1/E(zg),zg,initial=0)
    return np.interp(z,zg,ch)
def mu_model(z,H):
    z=np.asarray(z,float); d=(1+z)*chi_of_z(z,H)
    return 5*np.log10(np.maximum(d,1e-12))+25
def unitvec(ra,dec):
    a=np.deg2rad(np.asarray(ra,float)); d=np.deg2rad(np.asarray(dec,float)); q=np.cos(d)
    return np.c_[q*np.cos(a),q*np.sin(a),np.sin(d)]
def shell_delta(d1,dcum,R1,R2):
    return (dcum*R2**3-d1*R1**3)/(R2**3-R1**3)
def enclosed_delta(r,R1,R2,d1,d2,boundary,compfac):
    r=np.asarray(r,float); out=np.zeros_like(r)
    # enclosed deficit mass in units 4pi/3 rho_bar; Mdelta = delta * R^3 pieces.
    m1=d1*R1**3
    m2=m1+d2*(R2**3-R1**3)
    x=r>0
    rr=r[x]; val=np.zeros_like(rr)
    q1=rr<=R1; val[q1]=d1
    q2=(rr>R1)&(rr<=R2); val[q2]=(m1+d2*(rr[q2]**3-R1**3))/rr[q2]**3
    q3=rr>R2
    if boundary=='uncompensated':
        val[q3]=m2/rr[q3]**3
    else:
        R3=compfac*R2
        d3=-m2/(R3**3-R2**3)
        mid=q3&(rr<=R3); far=rr>R3
        val[mid]=(m2+d3*(rr[mid]**3-R2**3))/rr[mid]**3
        val[far]=0.0
    out[x]=val
    return out
def vfield(pos,H,R1,R2,d1,d2,boundary,compfac):
    r=np.linalg.norm(pos,axis=-1)
    db=enclosed_delta(r,R1,R2,d1,d2,boundary,compfac)
    f=OM**.55
    # vector v = -(H f / 3) * delta_bar(<r) * position
    return -(H*f/3.0)*db[...,None]*pos
def fit_h_analytic(zapp,mutrue,Hbg):
    mub=mu_model(zapp,Hbg)
    c=float(np.mean(mutrue-mub))
    return Hbg*10**(-c/5.0)
def simulate(df,Hbg,bgal,pname,d1g,dcumg,D,bname,compfac):
    h=Hbg/100.; R1=100./h; R2=150./h
    d1=d1g/bgal; dcum=dcumg/bgal; d2=shell_delta(d1,dcum,R1,R2)
    m=(df.IS_CALIBRATOR==0)&(df.zHD>=.023)&(df.zHD<.15)&np.isfinite(df.RA)&np.isfinite(df.DEC)
    q=df.loc[m]; z=q.zHD.to_numpy(float); U=unitvec(q.RA,q.DEC); rr=chi_of_z(z,Hbg); mut=mu_model(z,Hbg)
    vals=[]; dip=[]
    axes=[np.array([1.,0.,0.])] if D==0 else []
    if D>0:
        for _ in range(NORI):
            a=RNG.normal(size=3); a/=np.linalg.norm(a); axes.append(a)
    for a in axes:
        obs=a*D; src=obs[None,:]+U*rr[:,None]
        vo=vfield(obs[None,:],Hbg,R1,R2,d1,d2,bname,compfac)[0]
        vs=vfield(src,Hbg,R1,R2,d1,d2,bname,compfac)
        vlos=np.sum((vs-vo[None,:])*U,axis=1)
        zapp=z+(1+z)*vlos/C
        hf=fit_h_analytic(zapp,mut,Hbg)
        vals.append(hf)
        # diagnostic amplitude of line-of-sight velocity dipole regression (km/s)
        X=np.c_[np.ones(len(U)),U]
        beta=np.linalg.lstsq(X,vlos,rcond=None)[0]
        dip.append(float(np.linalg.norm(beta[1:])))
    vals=np.asarray(vals); dip=np.asarray(dip)
    return dict(H0_background=Hbg,galaxy_bias=bgal,profile=pname,delta_g_inner=d1g,delta_g_cumulative_R2=dcumg,
                R1_Mpc=R1,R2_Mpc=R2,delta_m_inner=d1,delta_m_shell=d2,boundary=bname,compensation_factor=compfac if compfac else np.nan,
                observer_offset_Mpc=D,N_SN=int(len(q)),N_orient=int(len(vals)),
                H0_app_mean=float(vals.mean()),H0_app_p05=float(np.quantile(vals,.05)),H0_app_p50=float(np.quantile(vals,.5)),H0_app_p95=float(np.quantile(vals,.95)),
                boost_mean_percent=float(100*(vals.mean()/Hbg-1)),boost_p05_percent=float(100*(np.quantile(vals,.05)/Hbg-1)),boost_p95_percent=float(100*(np.quantile(vals,.95)/Hbg-1)),
                vlos_dipole_mean_kms=float(dip.mean()),vlos_dipole_p95_kms=float(np.quantile(dip,.95)))
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--data',required=True); ap.add_argument('--out',default='publishedhole'); a=ap.parse_args()
    o=Path(a.out);o.mkdir(parents=True,exist_ok=True); df=pd.read_csv(a.data,sep=r'\s+')
    rows=[]
    for H in HBACK:
      for bg in BGAL:
       for pn,d1,dc in PROFILES:
        for bn,cf in BOUNDARIES:
         for D in OFFSETS: rows.append(simulate(df,H,bg,pn,d1,dc,D,bn,cf))
    R=pd.DataFrame(rows);R.to_csv(o/'published_localhole_flow_grid.csv',index=False)
    # frozen nominal rows
    P=R[(R.galaxy_bias==1.1)&(R.profile=='nominal')].copy();P.to_csv(o/'nominal_b1p1.csv',index=False)
    # closest rows to H0=73 per background, but do not treat as fit. Pure readback.
    clos=[]
    for H in HBACK:
        s=R[R.H0_background==H].copy();s['abs_to_73']=abs(s.H0_app_mean-73.0);clos.append(s.sort_values('abs_to_73').head(12))
    Cc=pd.concat(clos);Cc.to_csv(o/'closest_to_73_diagnostic.csv',index=False)
    # robust envelope of nominal-profile predictions at modest offsets <=50 Mpc, b 1-1.2 across boundary conditions.
    env={}
    for H in HBACK:
        s=R[(R.H0_background==H)&(R.profile=='nominal')&(R.observer_offset_Mpc<=50)]
        env[str(H)]={'H0_min_mean':float(s.H0_app_mean.min()),'H0_max_mean':float(s.H0_app_mean.max()),'H0_median_mean':float(s.H0_app_mean.median()),
                     'boost_min_percent':float(s.boost_mean_percent.min()),'boost_max_percent':float(s.boost_mean_percent.max()),
                     'fraction_mean_ge_72p5':float((s.H0_app_mean>=72.5).mean()),'fraction_mean_ge_73':float((s.H0_app_mean>=73).mean())}
    # Nominal exact center, b=1.0/1.1/1.2 by boundary.
    core=R[(R.profile=='nominal')&(R.observer_offset_Mpc==0)].copy()
    summary={'input_profile':'published Local Hole envelope: delta_g~-0.20 to 100 h^-1 Mpc and cumulative ~-0.13 to 150 h^-1 Mpc',
             'velocity_model':'linear spherical v=-(H f/3) delta_bar(<r) r; includes observer velocity and exterior flow',
             'robust_nominal_envelope_offsets_le_50Mpc':env,
             'nominal_center_rows':core.to_dict(orient='records')}
    (o/'summary.json').write_text(json.dumps(summary,indent=2,default=float))
    print(json.dumps(summary,indent=2,default=float)); print('\nNOMINAL b=1.1\n',P.to_string(index=False)); print('\nCLOSEST DIAGNOSTIC\n',Cc.to_string(index=False))
if __name__=='__main__':main()

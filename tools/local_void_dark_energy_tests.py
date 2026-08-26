#!/usr/bin/env python3
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize
from scipy.linalg import cho_factor, cho_solve

C=299792.458
R_CUTS=np.array([40.,70.,100.,150.,300.])


def E(z,om,w0=-1.,wa=0.):
    z=np.asarray(z,float); a=1/(1+z)
    de=(1-om)*a**(-3*(1+w0+wa))*np.exp(-3*wa*(1-a))
    return np.sqrt(om*(1+z)**3+de)

def mu_model(z,H0,om,w0=-1.,wa=0.):
    z=np.asarray(z,float)
    zg=np.linspace(0,max(2.5,float(np.max(z))*1.01),20000)
    chi=(C/H0)*cumulative_trapezoid(1/E(zg,om,w0,wa),zg,initial=0)
    d=(1+z)*np.interp(z,zg,chi)
    return 5*np.log10(np.maximum(d,1e-9))+25

def comoving_z(R,H0=70,om=.3):
    zg=np.linspace(0,.15,20000); chi=(C/H0)*cumulative_trapezoid(1/E(zg,om),zg,initial=0)
    return float(np.interp(R,chi,zg))

def load_cov(path,n):
    vals=np.loadtxt(path)
    if vals.ndim==1 and len(vals)==n*n+1 and int(vals[0])==n:
        vals=vals[1:]
    elif vals.ndim==1 and len(vals)==n*n:
        pass
    else:
        vals=vals.reshape(-1)
        if len(vals)==n*n+1: vals=vals[1:]
    return vals.reshape(n,n)

def smooth_step(z,zc,width_frac=.12):
    w=max(zc*width_frac,0.002)
    return .5*(1-np.tanh((z-zc)/w))

def hess_err(fun,x,i,step):
    x=np.array(x,float); f0=fun(x)
    xp=x.copy();xm=x.copy();xp[i]+=step;xm[i]-=step
    d2=(fun(xp)-2*f0+fun(xm))/(step*step)
    return math.sqrt(2/d2) if d2>0 else np.nan

def fit_base(z,mu,cf,om0=.31, ztype='zHD'):
    def f(p):
        H0,om=p
        if not (50<H0<90 and .1<om<.6): return 1e50
        r=mu-mu_model(z,H0,om)
        return float(r@cho_solve(cf,r))
    res=minimize(f,[73,.31],method='Nelder-Mead',options={'maxiter':3000,'xatol':1e-8,'fatol':1e-5})
    return res, f

def fit_step(z,mu,cf,zc,base):
    def f(p):
        H0,om,A=p
        if not (50<H0<90 and .1<om<.6 and -.5<A<.5): return 1e50
        r=mu-(mu_model(z,H0,om)+A*smooth_step(z,zc))
        return float(r@cho_solve(cf,r))
    x0=[base.x[0],base.x[1],0]
    res=minimize(f,x0,method='Nelder-Mead',options={'maxiter':4000,'xatol':1e-8,'fatol':1e-5})
    err=hess_err(f,res.x,2,0.002)
    return res,err

def fit_cpl(z,mu,cf,omfix=.3114, corr=None):
    y=mu.copy() if corr is None else mu-corr
    def f(p):
        H0,w0,wa=p
        if not (50<H0<90 and -2.5<w0<0 and -3<wa<2 and w0+wa<0): return 1e50
        r=y-mu_model(z,H0,omfix,w0,wa)
        return float(r@cho_solve(cf,r))
    res=minimize(f,[72.5,-.9,-.2],method='Nelder-Mead',options={'maxiter':6000,'xatol':1e-8,'fatol':1e-5})
    return res

def fit_shell_H0(z,mu,Csub,om,lo,hi):
    m=(z>=lo)&(z<hi)
    if m.sum()<8:return None
    cc=Csub[np.ix_(m,m)]
    try: cf=cho_factor(cc,lower=True,check_finite=False)
    except Exception:return None
    zz=z[m]; yy=mu[m]
    def f(x):
        H=float(x[0]);r=yy-mu_model(zz,H,om);return float(r@cho_solve(cf,r))
    r=minimize(f,[73],method='Nelder-Mead')
    return dict(zlo=lo,zhi=hi,n=int(m.sum()),H0=float(r.x[0]),chi2=float(r.fun))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--data',required=True);ap.add_argument('--cov',required=True);ap.add_argument('--out',default='voidde');a=ap.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    d=pd.read_csv(a.data,delim_whitespace=True)
    Cfull=load_cov(a.cov,len(d))
    rows=[]; summaries={}
    for zcol in ['zHD','zCMB']:
        zall=d[zcol].to_numpy(float); muall=d.MU_SH0ES.to_numpy(float)
        mask=(d.IS_CALIBRATOR.to_numpy(int)==0)&np.isfinite(zall)&np.isfinite(muall)&(zall>0.005)&(zall<2.3)
        idx=np.where(mask)[0]; z=zall[idx]; mu=muall[idx]; cov=Cfull[np.ix_(idx,idx)]
        # tiny jitter only if needed
        for eps in [0,1e-10,1e-8,1e-6]:
            try: cf=cho_factor(cov+np.eye(len(cov))*eps,lower=True,check_finite=False);break
            except Exception: continue
        base,basefun=fit_base(z,mu,cf)
        H0b,omb=map(float,base.x)
        step_rows=[]
        best=None
        for R in R_CUTS:
            zc=comoving_z(R,H0b,omb)
            rr,Aerr=fit_step(z,mu,cf,zc,base)
            H0,om,A=map(float,rr.x); hrat=10**(-A/5); f0=omb**.55; delta=-3*(hrat-1)/f0
            row=dict(z_definition=zcol,Rcut_Mpc=float(R),zcut=zc,N_inside=int((z<zc).sum()),H0_background=H0,Omega_m=om,A_mag=A,A_err=Aerr,A_sig=A/Aerr if np.isfinite(Aerr) and Aerr>0 else np.nan,delta_chi2=float(rr.fun-base.fun),Hlocal_over_Hbg=hrat,deltaH_percent=100*(hrat-1),linear_delta_density=delta)
            rows.append(row);step_rows.append(row)
            if best is None or rr.fun<best[0]:best=(rr.fun,row)
        # shell H0 with global Om
        shells=[]
        for lo,hi in [(0.005,.01),(.01,.023),(.023,.05),(.05,.1),(.1,.15),(.15,.3),(.3,.6)]:
            x=fit_shell_H0(z,mu,cov,omb,lo,hi)
            if x: x['z_definition']=zcol;shells.append(x)
        pd.DataFrame(shells).to_csv(out/f'shell_H0_{zcol}.csv',index=False)
        # CPL fit, before/after best void correction. Use fixed Om=DESI+CMB+Pantheon+ mean.
        cpl0=fit_cpl(z,mu,cf,.3114)
        brow=best[1]; corr=brow['A_mag']*smooth_step(z,brow['zcut'])
        cpl1=fit_cpl(z,mu,cf,.3114,corr=corr)
        # compare shapes of published DESI DR2 CPL vs LCDM, normalized at z=0.5 to remove H0/M offset
        zz=np.array([.01,.02,.03,.05,.07,.1,.2,.3,.5,.8,1.0])
        mu_l=mu_model(zz,68.17,.3027,-1,0); mu_c=mu_model(zz,67.51,.3114,-.838,-.62)
        shape=(mu_c-mu_l); shape-=np.interp(.5,zz,shape)
        void_curve=brow['A_mag']*smooth_step(zz,brow['zcut']); void_curve-=np.interp(.5,zz,void_curve)
        pd.DataFrame({'z':zz,'DESI_CPL_minus_LCDM_mag_norm_z0p5':shape,'best_void_mag_norm_z0p5':void_curve}).to_csv(out/f'cpl_vs_void_{zcol}.csv',index=False)
        summaries[zcol]=dict(N=int(len(z)),baseline=dict(H0=H0b,Omega_m=omb,chi2=float(base.fun)),best_step=brow,cpl_fixed_Om_before=dict(H0=float(cpl0.x[0]),w0=float(cpl0.x[1]),wa=float(cpl0.x[2]),chi2=float(cpl0.fun)),cpl_fixed_Om_after_void_correction=dict(H0=float(cpl1.x[0]),w0=float(cpl1.x[1]),wa=float(cpl1.x[2]),chi2=float(cpl1.fun)),cpl_shift=dict(dH0=float(cpl1.x[0]-cpl0.x[0]),dw0=float(cpl1.x[1]-cpl0.x[1]),dwa=float(cpl1.x[2]-cpl0.x[2])))
    pd.DataFrame(rows).to_csv(out/'step_tests.csv',index=False)
    (out/'summary.json').write_text(json.dumps(summaries,indent=2))
    print(json.dumps(summaries,indent=2));print(pd.DataFrame(rows).to_string(index=False))
if __name__=='__main__':main()

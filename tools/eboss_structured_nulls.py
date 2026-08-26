#!/usr/bin/env python3
import argparse, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.integrate import cumulative_trapezoid
from scipy.spatial import cKDTree

H0=67.66; OM=0.3111; OL=1-OM; C=299792.458
RNG=np.random.default_rng(82620261335)
ZG=np.linspace(0,1.3,25000)
DCG=(C/H0)*cumulative_trapezoid(1/np.sqrt(OM*(1+ZG)**3+OL),ZG,initial=0)
def dc(z): return np.interp(np.asarray(z,float),ZG,DCG)
def xyz(ra,dec,z):
    r=dc(z); a=np.deg2rad(ra); d=np.deg2rad(dec); cd=np.cos(d)
    return np.c_[r*cd*np.cos(a),r*cd*np.sin(a),r*np.sin(d)]
def cx(r): return xyz([r.ra_deg],[r.dec_deg],[r.z])[0]

def getcol(t,*names):
    lut={n.upper():n for n in t.names}
    for n in names:
        if n.upper() in lut: return np.asarray(t[lut[n.upper()]])
    raise KeyError((names,t.names))

def load(paths,zmin=.6,zmax=1.0):
    chunks=[]
    for reg,p in enumerate(paths):
        with fits.open(p,memmap=True) as h:
            t=h[1].data
            ra=getcol(t,'RA','RA_TARGET'); de=getcol(t,'DEC','DEC_TARGET'); z=getcol(t,'Z','Z_REDROCK')
            m=np.isfinite(ra)&np.isfinite(de)&np.isfinite(z)&(z>=zmin)&(z<=zmax)
            chunks.append(pd.DataFrame({'ra':ra[m].astype(float),'dec':de[m].astype(float),'z':z[m].astype(float),'region':reg}))
    return pd.concat(chunks,ignore_index=True)

def onep(vals,obs,kind):
    return (np.sum(vals<=obs)+1)/(len(vals)+1) if kind=='VOID' else (np.sum(vals>=obs)+1)/(len(vals)+1)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--data',nargs=2,required=True); ap.add_argument('--random',nargs=2,required=True); ap.add_argument('--cand',required=True); ap.add_argument('--out',default='ebossnull'); ap.add_argument('--nrep',type=int,default=100)
    a=ap.parse_args(); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    dat=load(a.data); ran=load(a.random); cands=pd.read_csv(a.cand); cands=cands[(cands.z>=.6)&(cands.z<=1.0)].reset_index(drop=True)
    dxyz=xyz(dat.ra.values,dat.dec.values,dat.z.values); rxyz=xyz(ran.ra.values,ran.dec.values,ran.z.values)
    dt=cKDTree(dxyz); rt=cKDTree(rxyz)
    rows=[]
    for _,r in cands.iterrows():
        c=cx(r); R=float(r.R_mpc); obs=float(dt.query_ball_point(c,R,return_length=True))
        # Assign footprint region from nearest random point in each cap.
        nearest=[]
        for reg in [0,1]:
            ii=np.where(ran.region.values==reg)[0]; tr=cKDTree(rxyz[ii]); nearest.append((float(tr.query(c)[0]),reg))
        nd,reg=min(nearest); dr=dat[dat.region==reg]; rr=ran[ran.region==reg]
        # Direct random-catalog expectation with local redshift normalization.
        dz=.10
        dm=(dr.z.values>=r.z-dz)&(dr.z.values<=r.z+dz); rm=(rr.z.values>=r.z-dz)&(rr.z.values<=r.z+dz)
        alpha=dm.sum()/max(rm.sum(),1)
        ridx=np.where(ran.region.values==reg)[0]
        rc=float(cKDTree(rxyz[ridx]).query_ball_point(c,R,return_length=True))
        exp=alpha*rc
        rows.append(dict(candidate_id=r.candidate_id,kind=r.kind,family='EBOSS_OFFICIAL_RANDOM',observed=obs,null_mean=exp,ratio=obs/max(exp,1),zscore=(obs-exp)/math.sqrt(max(exp,1)),empirical_p=np.nan,coverage_random_count=rc,nearest_random_mpc=nd,region=reg))
        # Same-z random-footprint centers: real eBOSS universe, valid selection locations.
        pool=rr[np.abs(rr.z-r.z)<.025]
        if len(pool)>0:
            take=RNG.choice(len(pool),min(1200,len(pool)),replace=False); pp=pool.iloc[take]
            cc=xyz(pp.ra.values,pp.dec.values,np.full(len(pp),r.z)); cc=cc[np.linalg.norm(cc-c,axis=1)>2*R]
            vals=np.array([dt.query_ball_point(q,R,return_length=True) for q in cc],float)
            if len(vals):
                mu=vals.mean(); sd=max(vals.std(ddof=1),1)
                rows.append(dict(candidate_id=r.candidate_id,kind=r.kind,family='EBOSS_SAMEZ_RANDOM_SKY',observed=obs,null_mean=mu,ratio=obs/max(mu,1),zscore=(obs-mu)/sd,empirical_p=onep(vals,obs,r.kind),coverage_random_count=len(vals),nearest_random_mpc=nd,region=reg))
        # Same sky, nearby z: angular mask is identical to the candidate location.
        zpool=rr[(np.abs(rr.z-r.z)>.03)&(np.abs(rr.z-r.z)<.16)].z.values
        if len(zpool):
            zz=RNG.choice(zpool,min(1000,len(zpool)),replace=False); cc=xyz(np.full(len(zz),r.ra_deg),np.full(len(zz),r.dec_deg),zz)
            vals=np.array([dt.query_ball_point(q,R,return_length=True) for q in cc],float)
            mu=vals.mean(); sd=max(vals.std(ddof=1),1)
            rows.append(dict(candidate_id=r.candidate_id,kind=r.kind,family='EBOSS_SAMESKY_NEARBYZ',observed=obs,null_mean=mu,ratio=obs/max(mu,1),zscore=(obs-mu)/sd,empirical_p=onep(vals,obs,r.kind),coverage_random_count=len(vals),nearest_random_mpc=nd,region=reg))
    # Catalogue counterfactual nulls shared across candidates.
    obs=np.array([dt.query_ball_point(cx(r),float(r.R_mpc),return_length=True) for _,r in cands.iterrows()],float)
    for fam in ['REGION_Z_SHUFFLE','ZBIN0025_SKY']:
        arr=[]
        ra=dat.ra.values.copy(); de=dat.dec.values.copy(); z=dat.z.values.copy(); reg=dat.region.values.copy()
        for k in range(a.nrep):
            if fam=='REGION_Z_SHUFFLE':
                z2=z.copy()
                for g in [0,1]:
                    ii=np.where(reg==g)[0]; z2[ii]=z[RNG.permutation(ii)]
                tx=xyz(ra,de,z2)
            else:
                ra2=ra.copy(); de2=de.copy(); keys=reg*1000+np.floor(z/.025).astype(int)
                for key in np.unique(keys):
                    ii=np.where(keys==key)[0]
                    if len(ii)>1:
                        jj=RNG.permutation(ii); ra2[ii]=ra[jj]; de2[ii]=de[jj]
                tx=xyz(ra2,de2,z)
            tr=cKDTree(tx); arr.append([tr.query_ball_point(cx(r),float(r.R_mpc),return_length=True) for _,r in cands.iterrows()])
        arr=np.asarray(arr,float)
        for j,(_,r) in enumerate(cands.iterrows()):
            vals=arr[:,j]; mu=vals.mean(); sd=max(vals.std(ddof=1),1)
            rows.append(dict(candidate_id=r.candidate_id,kind=r.kind,family='EBOSS_'+fam,observed=obs[j],null_mean=mu,ratio=obs[j]/max(mu,1),zscore=(obs[j]-mu)/sd,empirical_p=onep(vals,obs[j],r.kind),coverage_random_count=np.nan,nearest_random_mpc=np.nan,region=np.nan))
    df=pd.DataFrame(rows); df.to_csv(out/'eboss_structured.csv',index=False)
    verdict=[]
    for cid,g in df.groupby('candidate_id'):
        kind=g.kind.iloc[0]; valid=g[(g.family!='EBOSS_OFFICIAL_RANDOM') | (g.coverage_random_count>=20)]
        if kind=='VOID': ok=bool((valid.ratio<.90).all())
        else: ok=bool((valid.ratio>1.10).all())
        verdict.append({'candidate_id':cid,'kind':kind,'families':len(valid),'survives':ok,'worst_ratio':float(valid.ratio.max() if kind=='VOID' else valid.ratio.min()),'max_empirical_p':float(valid.empirical_p.dropna().max()) if valid.empirical_p.notna().any() else None})
    vd=pd.DataFrame(verdict); vd.to_csv(out/'verdict.csv',index=False)
    print('VERDICT'); print(vd.to_string(index=False)); print('DETAIL'); print(df.to_string(index=False))
    (out/'meta.json').write_text(json.dumps({'data_objects':len(dat),'random_objects':len(ran),'nrep':a.nrep,'seed':82620261335},indent=2))
if __name__=='__main__': main()

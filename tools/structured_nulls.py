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
ZG=np.linspace(0,1.6,30000)
DCG=(C/H0)*cumulative_trapezoid(1/np.sqrt(OM*(1+ZG)**3+OL),ZG,initial=0)

def dc(z): return np.interp(np.asarray(z,float),ZG,DCG)
def xyz(ra,dec,z):
    r=dc(z); a=np.deg2rad(ra); d=np.deg2rad(dec); cd=np.cos(d)
    return np.c_[r*cd*np.cos(a), r*cd*np.sin(a), r*np.sin(d)]
def center_xyz(row): return xyz([row.ra_deg],[row.dec_deg],[row.z])[0]

def load(paths,zmin,zmax):
    out=[]
    for region,p in enumerate(paths):
        with fits.open(p,memmap=True) as h:
            t=h[1].data; names=set(t.names)
            ra=np.asarray(t['RA'],float); de=np.asarray(t['DEC'],float); z=np.asarray(t['Z'],float)
            m=np.isfinite(ra)&np.isfinite(de)&np.isfinite(z)&(z>=zmin)&(z<=zmax)
            ph=np.asarray(t['PHOTSYS']).astype(str)[m] if 'PHOTSYS' in names else np.full(m.sum(),'X')
            out.append(pd.DataFrame({'ra':ra[m],'dec':de[m],'z':z[m],'region':region,'photsys':ph}))
    return pd.concat(out,ignore_index=True)

def count_candidates(tree,cands):
    return np.array([tree.query_ball_point(center_xyz(r),float(r.R_mpc),return_length=True) for _,r in cands.iterrows()],float)

def one_sided_p(null,obs,kind):
    if kind=='VOID': return (np.sum(null<=obs)+1)/(len(null)+1)
    return (np.sum(null>=obs)+1)/(len(null)+1)

def summaries(family,arr,obs,cands,tracer):
    rows=[]
    for j,(_,r) in enumerate(cands.iterrows()):
        x=np.asarray(arr)[:,j]
        mu=float(x.mean()); sd=float(max(x.std(ddof=1),1.0))
        rows.append(dict(tracer=tracer,candidate_id=r.candidate_id,kind=r.kind,family=family,
                         observed=float(obs[j]),null_mean=mu,null_sd=sd,ratio=float(obs[j]/mu),
                         zscore=float((obs[j]-mu)/sd),empirical_p=float(one_sided_p(x,obs[j],r.kind)),
                         q05=float(np.quantile(x,.05)),q50=float(np.quantile(x,.5)),q95=float(np.quantile(x,.95))))
    return rows

def permute_by_groups(values,keys,mode='permute'):
    out=values.copy()
    order=np.argsort(keys,kind='stable'); sk=keys[order]
    cuts=np.r_[0,np.flatnonzero(sk[1:]!=sk[:-1])+1,len(order)]
    for a,b in zip(cuts[:-1],cuts[1:]):
        ids=order[a:b]
        if len(ids)<2: continue
        if mode=='permute': out[ids]=values[RNG.permutation(ids)]
        else:
            k=int(RNG.integers(1,len(ids)))
            out[ids]=np.roll(values[ids],k)
    return out

def transformed_counts(cat,cands,family,nrep):
    ra=cat.ra.to_numpy(float); de=cat.dec.to_numpy(float); z=cat.z.to_numpy(float)
    reg=cat.region.to_numpy(int); ph=cat.photsys.astype(str).to_numpy()
    out=[]
    if family=='REGION_PHOTSYS_Z':
        _,phcode=np.unique(ph,return_inverse=True); keys=reg*10+phcode
    elif family=='ANGCELL10x5_Z':
        keys=reg*100000 + (ra//10).astype(int)*100 + ((de+90)//5).astype(int)
    elif family=='ANGCELL5x5_Z':
        keys=reg*1000000 + (ra//5).astype(int)*1000 + ((de+90)//5).astype(int)
    elif family=='CELL10x5_ROLL':
        keys=reg*100000 + (ra//10).astype(int)*100 + ((de+90)//5).astype(int)
    elif family=='ZBIN0025_SKY':
        keys=reg*1000 + np.floor(z/0.025).astype(int)
    else:
        raise ValueError(family)
    for i in range(nrep):
        if family=='ZBIN0025_SKY':
            idx=np.arange(len(z)); pidx=permute_by_groups(idx,keys,'permute').astype(int)
            tx=xyz(ra[pidx],de[pidx],z)
        else:
            mode='roll' if family=='CELL10x5_ROLL' else 'permute'
            z2=permute_by_groups(z,keys,mode)
            tx=xyz(ra,de,z2)
        tr=cKDTree(tx)
        out.append(count_candidates(tr,cands))
        if (i+1)%10==0: print(f'{family} {i+1}/{nrep}',flush=True)
    return np.asarray(out)

def infer_region(cat,cands):
    regs=[]
    for _,r in cands.iterrows():
        c=center_xyz(r); ds=[]
        for reg in sorted(cat.region.unique()):
            sub=cat[cat.region==reg]
            tr=cKDTree(xyz(sub.ra.values,sub.dec.values,sub.z.values))
            ds.append((tr.query(c,k=1)[0],reg))
        regs.append(min(ds)[1])
    return regs

def conditional_controls(cat,cands,obs,nctrl=1000):
    tr=cKDTree(xyz(cat.ra.values,cat.dec.values,cat.z.values))
    regs=infer_region(cat,cands); rows=[]
    for j,(_,r) in enumerate(cands.iterrows()):
        reg=regs[j]; R=float(r.R_mpc); c0=center_xyz(r)
        sub=cat[cat.region==reg]
        # Same-z, different-sky: preserves the observed universe and the redshift slice.
        s=sub[np.abs(sub.z-r.z)<0.025]
        if len(s)>0:
            ids=RNG.choice(len(s),min(nctrl,len(s)),replace=False)
            cc=xyz(s.ra.values[ids],s.dec.values[ids],np.full(len(ids),r.z))
            far=np.linalg.norm(cc-c0,axis=1)>2*R
            cc=cc[far]
            vals=np.array([tr.query_ball_point(q,R,return_length=True) for q in cc],float)
            if len(vals):
                mu=vals.mean(); sd=max(vals.std(ddof=1),1)
                rows.append(dict(tracer='',candidate_id=r.candidate_id,kind=r.kind,family='REAL_UNIVERSE_SAMEZ_SKY',observed=obs[j],null_mean=mu,null_sd=sd,ratio=obs[j]/mu,zscore=(obs[j]-mu)/sd,empirical_p=one_sided_p(vals,obs[j],r.kind),q05=np.quantile(vals,.05),q50=np.quantile(vals,.5),q95=np.quantile(vals,.95)))
        # Same sky, nearby redshift controls: preserves angular selection exactly.
        pool=sub[(np.abs(sub.z-r.z)>0.03)&(np.abs(sub.z-r.z)<0.18)]
        if len(pool)>0:
            zz=RNG.choice(pool.z.values,min(nctrl,len(pool)),replace=False)
            cc=xyz(np.full(len(zz),r.ra_deg),np.full(len(zz),r.dec_deg),zz)
            vals=np.array([tr.query_ball_point(q,R,return_length=True) for q in cc],float)
            mu=vals.mean(); sd=max(vals.std(ddof=1),1)
            rows.append(dict(tracer='',candidate_id=r.candidate_id,kind=r.kind,family='REAL_UNIVERSE_SAMESKY_Z',observed=obs[j],null_mean=mu,null_sd=sd,ratio=obs[j]/mu,zscore=(obs[j]-mu)/sd,empirical_p=one_sided_p(vals,obs[j],r.kind),q05=np.quantile(vals,.05),q50=np.quantile(vals,.5),q95=np.quantile(vals,.95)))
    return rows

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--tracer',choices=['LRG','ELG'],required=True); ap.add_argument('--ngc',required=True); ap.add_argument('--sgc',required=True); ap.add_argument('--candidates',required=True); ap.add_argument('--nrep',type=int,default=50); ap.add_argument('--outdir',default='structured')
    a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    zmin,zmax=(0.4,1.1) if a.tracer=='LRG' else (0.8,1.1)
    cat=load([a.ngc,a.sgc],zmin,zmax); cands=pd.read_csv(a.candidates)
    cands=cands[(cands.z>=zmin)&(cands.z<=zmax)].reset_index(drop=True)
    realtree=cKDTree(xyz(cat.ra.values,cat.dec.values,cat.z.values)); obs=count_candidates(realtree,cands)
    rows=[]
    for x in conditional_controls(cat,cands,obs,1000): x['tracer']=a.tracer; rows.append(x)
    families=['REGION_PHOTSYS_Z','ANGCELL10x5_Z','ANGCELL5x5_Z','CELL10x5_ROLL','ZBIN0025_SKY']
    for fam in families:
        arr=transformed_counts(cat,cands,fam,a.nrep)
        rows.extend(summaries(fam,arr,obs,cands,a.tracer))
    df=pd.DataFrame(rows)
    df.to_csv(out/f'structured_{a.tracer}.csv',index=False)
    # Frozen verdict: candidate must have directionally anomalous ratio in every structured family;
    # p is descriptive because these are conditional nulls, not a global discovery p-value.
    v=[]
    for cid,g in df.groupby('candidate_id'):
        kind=g.kind.iloc[0]
        if kind=='VOID': ok=bool((g.ratio<0.90).all())
        else: ok=bool((g.ratio>1.10).all())
        v.append({'tracer':a.tracer,'candidate_id':cid,'kind':kind,'families':len(g),'survives_all_structured':ok,'max_empirical_p':float(g.empirical_p.max()),'worst_ratio':float(g.ratio.max() if kind=='VOID' else g.ratio.min())})
    vd=pd.DataFrame(v); vd.to_csv(out/f'verdict_{a.tracer}.csv',index=False)
    meta={'tracer':a.tracer,'objects':len(cat),'nrep_per_transformed_family':a.nrep,'families':families+['REAL_UNIVERSE_SAMEZ_SKY','REAL_UNIVERSE_SAMESKY_Z'],'seed':82620261335}
    (out/f'meta_{a.tracer}.json').write_text(json.dumps(meta,indent=2))
    print(vd.to_string(index=False)); print(df.to_string(index=False))
if __name__=='__main__': main()

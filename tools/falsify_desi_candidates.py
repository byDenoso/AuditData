#!/usr/bin/env python3
import argparse, json, math
from pathlib import Path
import numpy as np, pandas as pd
from astropy.io import fits
from scipy.integrate import cumulative_trapezoid
from scipy.spatial import cKDTree

H0=67.66; OM=0.3111; OL=1-OM; C=299792.458; RNG=np.random.default_rng(82620261251)
ZG=np.linspace(0,1.6,30000); DCG=(C/H0)*cumulative_trapezoid(1/np.sqrt(OM*(1+ZG)**3+OL),ZG,initial=0)
def dc(z): return np.interp(z,ZG,DCG)
def xyz(ra,dec,z):
    r=dc(np.asarray(z,float)); a=np.deg2rad(ra); d=np.deg2rad(dec); cd=np.cos(d)
    return np.c_[r*cd*np.cos(a), r*cd*np.sin(a), r*np.sin(d)]
def center(ra,dec,z): return xyz(np.array([ra]),np.array([dec]),np.array([z]))[0]

def load(path,zmin=.4,zmax=1.1):
    with fits.open(path,memmap=True) as h:
        t=h[1].data; names=set(t.names); ra=np.asarray(t['RA'],float); de=np.asarray(t['DEC'],float); z=np.asarray(t['Z'],float)
        m=np.isfinite(z)&(z>=zmin)&(z<=zmax)
        out={'ra':ra[m],'dec':de[m],'z':z[m]}
        for k in ['WEIGHT','FRAC_TLOBS_TILES','PHOTSYS']:
            if k in names: out[k.lower()]=np.asarray(t[k])[m]
        return out

def cat_xyz(cat,mask=None):
    if mask is None: mask=np.ones(len(cat['z']),bool)
    return xyz(cat['ra'][mask],cat['dec'][mask],cat['z'][mask])

def count_ratio(cat,c,R,nnull=40,mask=None,weighted=False,local_bins=False):
    if mask is None: mask=np.ones(len(cat['z']),bool)
    ra,dec,z=cat['ra'][mask],cat['dec'][mask],cat['z'][mask]
    w=np.asarray(cat.get('weight',np.ones(len(cat['z']))))[mask].astype(float)
    pts=xyz(ra,dec,z); tr=cKDTree(pts)
    ids=tr.query_ball_point(c,R)
    obs=w[ids].sum() if weighted else len(ids)
    null=[]
    for _ in range(nnull):
        if local_bins:
            z2=z.copy(); bins=np.floor(z/0.05).astype(int)
            for b in np.unique(bins):
                ii=np.where(bins==b)[0]; z2[ii]=z[RNG.permutation(ii)]
        else: z2=z[RNG.permutation(len(z))]
        mt=cKDTree(xyz(ra,dec,z2)); jj=mt.query_ball_point(c,R)
        null.append(w[jj].sum() if weighted else len(jj))
    null=np.asarray(null,float); mu=null.mean(); sd=max(null.std(ddof=1),1)
    return obs/mu,(obs-mu)/sd,(np.sum(null <= obs)+1)/(len(null)+1),obs,mu

def perturbed_centers(c,offset=25):
    out=[c]
    for ax in range(3):
        for s in (-1,1):
            q=c.copy(); q[ax]+=s*offset; out.append(q)
    return out

def pca_shape(points):
    if len(points)<20:return (np.nan,np.nan,np.nan,np.nan)
    X=points-np.median(points,axis=0); vals,vec=np.linalg.eigh(np.cov(X,rowvar=False)); vec=vec[:,np.argsort(vals)[::-1]]; pr=X@vec
    ext=np.percentile(pr,95,axis=0)-np.percentile(pr,5,axis=0); a,b,c=ext
    return float(a),float(b),float(c),float(b/max(c,1e-9))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--v15-lrg-ngc',required=True); ap.add_argument('--v15-lrg-sgc',required=True); ap.add_argument('--v15-elg-ngc',required=True); ap.add_argument('--v15-elg-sgc',required=True); ap.add_argument('--v12-lrg-ngc'); ap.add_argument('--v12-lrg-sgc'); ap.add_argument('--candidates',required=True); ap.add_argument('--outdir',default='falsification')
    a=ap.parse_args(); out=Path(a.outdir); out.mkdir(exist_ok=True)
    cands=pd.read_csv(a.candidates)
    lrgN=load(a.v15_lrg_ngc); lrgS=load(a.v15_lrg_sgc); elgN=load(a.v15_elg_ngc,.8,1.1); elgS=load(a.v15_elg_sgc,.8,1.1)
    lrg12N=load(a.v12_lrg_ngc) if a.v12_lrg_ngc else None; lrg12S=load(a.v12_lrg_sgc) if a.v12_lrg_sgc else None
    rows=[]
    for _,r in cands.iterrows():
        kind=r['kind']; cc=center(r.ra_deg,r.dec_deg,r.z); R=float(r.R_mpc); hemi='N' if r.dec_deg>-10 else 'S'
        L=lrgN if hemi=='N' else lrgS; E=elgN if hemi=='N' else elgS; L12=lrg12N if hemi=='N' else lrg12S
        res={'candidate_id':r.candidate_id,'kind':kind,'ra_deg':r.ra_deg,'dec_deg':r.dec_deg,'z':r.z,'R_mpc':R}
        direction = -1 if kind=='VOID' else 1
        for label,cat in [('LRG',L),('ELG',E)]:
            if label=='ELG' and not (.8<=r.z<=1.1): continue
            q,zs,p,obs,mu=count_ratio(cat,cc,R,60,weighted=False); res[f'{label}_ratio']=q;res[f'{label}_z']=zs;res[f'{label}_p']=p
            qw,zsw,pw,_,_=count_ratio(cat,cc,R,60,weighted=True); res[f'{label}_wratio']=qw;res[f'{label}_wz']=zsw
            ql,zsl,pl,_,_=count_ratio(cat,cc,R,60,weighted=False,local_bins=True); res[f'{label}_localz_ratio']=ql;res[f'{label}_localz_z']=zsl
            # center perturbation stability, direct ratio versus fresh nulls per offset
            prs=[]
            for pc in perturbed_centers(cc,25): prs.append(count_ratio(cat,pc,R,20)[0])
            res[f'{label}_center_med_ratio']=float(np.median(prs)); res[f'{label}_center_worst_ratio']=float(max(prs) if kind=='VOID' else min(prs))
            # 50% thinning repeated, compare observed local count to full-catalog shuffled expectation via ratio
            thin=[]
            for _ in range(20):
                m=RNG.random(len(cat['z']))<.5; thin.append(count_ratio(cat,cc,R,8,mask=m)[0])
            res[f'{label}_thin_med_ratio']=float(np.median(thin));res[f'{label}_thin_q90_ratio']=float(np.quantile(thin,.9 if kind=='VOID' else .1))
            # completeness cuts if present
            if 'frac_tlobs_tiles' in cat:
                for cut in (.5,.8):
                    m=np.asarray(cat['frac_tlobs_tiles'],float)>=cut
                    if m.sum()>10000:
                        qc,zc,_,_,_=count_ratio(cat,cc,R,30,mask=m);res[f'{label}_frac{int(cut*10)}_ratio']=qc;res[f'{label}_frac{int(cut*10)}_z']=zc
            # photometric-system split
            if 'photsys' in cat:
                ph=np.asarray(cat['photsys']).astype(str)
                for ps in np.unique(ph):
                    m=ph==ps
                    if m.sum()>10000:
                        qp,zp,_,_,_=count_ratio(cat,cc,R,25,mask=m);res[f'{label}_PHOTSYS_{ps}_ratio']=qp;res[f'{label}_PHOTSYS_{ps}_z']=zp
        if L12 is not None:
            q12,z12,p12,_,_=count_ratio(L12,cc,R,60); res['LRG_v12_ratio']=q12;res['LRG_v12_z']=z12;res['LRG_v12_p']=p12
        # radial profile + geometry
        tr=cKDTree(cat_xyz(L)); radii=[50,75,100,130,160,200]
        for RR in radii:
            qR,zR,_,_,_=count_ratio(L,cc,RR,30);res[f'LRG_ratio_R{RR}']=qR;res[f'LRG_z_R{RR}']=zR
        if kind=='WALL':
            ids=tr.query_ball_point(cc,140); aa,bb,cc3,pl=pca_shape(tr.data[np.asarray(ids)]);res['axis1']=aa;res['axis2']=bb;res['axis3']=cc3;res['planarity']=pl
            shapes=[]
            for pc in perturbed_centers(cc,25):
                ids=tr.query_ball_point(pc,140); shapes.append(pca_shape(tr.data[np.asarray(ids)])[3])
            res['planarity_perturb_min']=float(np.nanmin(shapes));res['planarity_perturb_med']=float(np.nanmedian(shapes))
        rows.append(res)
    df=pd.DataFrame(rows);df.to_csv(out/'falsification_results.csv',index=False)
    # verdict rule frozen here: both tracers directional, local-z directional, weighted directional; version if available; center and thinning retain direction; wall planarity >1.4
    verdicts=[]
    for _,r in df.iterrows():
        if r.kind=='VOID':
            tests=[]
            for lab in ['LRG','ELG']:
                if f'{lab}_ratio' in r and pd.notna(r.get(f'{lab}_ratio')):
                    tests += [r[f'{lab}_ratio']<.80,r[f'{lab}_wratio']<.82,r[f'{lab}_localz_ratio']<.82,r[f'{lab}_center_worst_ratio']<.90,r[f'{lab}_thin_q90_ratio']<.90]
            if pd.notna(r.get('LRG_v12_ratio',np.nan)): tests.append(r.LRG_v12_ratio<.85)
            verdicts.append('SURVIVES' if all(tests) else 'FAIL')
        else:
            tests=[r.LRG_ratio>1.15,r.LRG_wratio>1.15,r.LRG_localz_ratio>1.12,r.LRG_center_worst_ratio>1.10,r.LRG_thin_q90_ratio>1.10]
            if pd.notna(r.get('ELG_ratio',np.nan)): tests += [r.ELG_ratio>1.10,r.ELG_wratio>1.10]
            tests += [r.get('planarity',0)>1.4,r.get('planarity_perturb_min',0)>1.25]
            verdicts.append('SURVIVES' if all(tests) else 'FAIL')
    df['verdict']=verdicts;df.to_csv(out/'falsification_results.csv',index=False)
    with open(out/'summary.json','w') as f: json.dump({'n':len(df),'survives':int((df.verdict=='SURVIVES').sum()),'fail':int((df.verdict=='FAIL').sum())},f,indent=2)
    print(df[['candidate_id','kind','verdict']].to_string(index=False)); print((out/'summary.json').read_text())
if __name__=='__main__': main()

#!/usr/bin/env python3
import argparse, json, math, gc
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.integrate import cumulative_trapezoid
from scipy.spatial import cKDTree

H0=67.66; OM=0.3111; OL=1-OM; C=299792.458; SEED=82620261
SCALE_RADII=np.array([50.,75.,100.,125.,140.,160.,200.])

zg=np.linspace(0,1.5,30001); inv=1/np.sqrt(OM*(1+zg)**3+OL)
dcg=(C/H0)*cumulative_trapezoid(inv,zg,initial=0)
def dc(z): return np.interp(z,zg,dcg)
def iz(r): return np.interp(r,dcg,zg)

def xyz(ra,de,z):
    r=dc(np.asarray(z,float)); a=np.deg2rad(ra); d=np.deg2rad(de); cd=np.cos(d)
    return np.column_stack((r*cd*np.cos(a),r*cd*np.sin(a),r*np.sin(d))).astype('f4')

def center_xyz(ra,de,z): return xyz(np.array([ra]),np.array([de]),np.array([z]))[0].astype(float)

def local_basis(q):
    u=q/np.linalg.norm(q)
    ref=np.array([0.,0.,1.])
    if abs(np.dot(u,ref))>.9: ref=np.array([0.,1.,0.])
    e1=np.cross(u,ref); e1/=np.linalg.norm(e1)
    e2=np.cross(u,e1); e2/=np.linalg.norm(e2)
    return u,e1,e2

def load_cat(paths,zmin,zmax):
    chunks=[]
    for region,p in enumerate(paths):
        print('LOAD',p,flush=True)
        with fits.open(p,memmap=True) as h:
            d=h[1].data; names={n.upper():n for n in d.names}
            ra=np.asarray(d[names['RA']],float); de=np.asarray(d[names['DEC']],float); z=np.asarray(d[names['Z']],float)
            good=np.isfinite(ra)&np.isfinite(de)&np.isfinite(z)&(z>=zmin)&(z<=zmax)
            if 'WEIGHT' in names: w=np.asarray(d[names['WEIGHT']],float)
            else:
                ws=[]
                for n in ['WEIGHT_SYSTOT','WEIGHT_CP','WEIGHT_NOZ']:
                    if n in names: ws.append(np.asarray(d[names[n]],float))
                w=np.prod(ws,axis=0) if ws else np.ones(len(d),float)
            w=np.where(np.isfinite(w)&(w>0),w,1.)
            ph=np.asarray(d[names['PHOTSYS']]).astype(str) if 'PHOTSYS' in names else np.full(len(d),'?',dtype='<U1')
            chunks.append((ra[good].astype('f4'),de[good].astype('f4'),z[good].astype('f4'),w[good].astype('f4'),ph[good],np.full(good.sum(),region,'i1')))
            print(' kept',int(good.sum()),flush=True)
    return tuple(np.concatenate([x[k] for x in chunks]) for k in range(6))

def shuffled_xyz(ra,de,z,reg,rng):
    zz=z.copy()
    for g in np.unique(reg):
        ii=np.where(reg==g)[0]; zz[ii]=z[rng.permutation(ii)]
    return xyz(ra,de,zz)

def count(tree,q,r): return int(tree.query_ball_point(q,float(r),return_length=True))
def ids(tree,q,r): return np.asarray(tree.query_ball_point(q,float(r)),dtype=int)

def build_geometries(cands):
    G={}
    for _,r in cands.iterrows():
        q=center_xyz(r.ra_deg,r.dec_deg,r.z); R=float(r.radius_mpc); kind=r.kind
        u,e1,e2=local_basis(q)
        jit=[]
        for a in (-.30,0,.30):
            for b in (-.30,0,.30):
                for c in (-.30,0,.30): jit.append(q+R*(a*u+b*e1+c*e2))
        radial=[q+s*R*u for s in (-1.,-.5,0,.5,1.)]
        ring=[q+1.5*R*(math.cos(t)*e1+math.sin(t)*e2) for t in np.linspace(0,2*math.pi,12,endpoint=False)]
        G[r.candidate_id]=dict(q=q,R=R,kind=kind,u=u,e1=e1,e2=e2,jitter=np.asarray(jit),radial=np.asarray(radial),ring=np.asarray(ring))
    return G

def octant_counts(X,tree,q,R,u,e1,e2):
    ii=ids(tree,q,R)
    if len(ii)==0:return np.zeros(8,int)
    D=X[ii]-q
    b0=(D@u>=0).astype(int); b1=(D@e1>=0).astype(int); b2=(D@e2>=0).astype(int)
    k=4*b0+2*b1+b2
    return np.bincount(k,minlength=8)

def evaluate_desi(cands,paths,zrange,nnull,label,rng):
    ra,de,z,w,ph,reg=load_cat(paths,*zrange); X=xyz(ra,de,z); tree=cKDTree(X); G=build_geometries(cands)
    w95=np.minimum(w,np.nanpercentile(w,95))
    obs={}
    for _,r in cands.iterrows():
        cid=r.candidate_id; g=G[cid]; q,R=g['q'],g['R']; ii=ids(tree,q,R)
        obs[cid]=dict(
            nominal=count(tree,q,R), weight=float(w[ii].sum()), weight95=float(w95[ii].sum()),
            scales=np.array([count(tree,q,rr) for rr in SCALE_RADII],float),
            jitter=np.array([count(tree,qq,R) for qq in g['jitter']],float),
            radial=np.array([count(tree,qq,R) for qq in g['radial']],float),
            ring=np.array([count(tree,qq,R) for qq in g['ring']],float),
            oct=octant_counts(X,tree,q,R,g['u'],g['e1'],g['e2']).astype(float),
            phN=float(np.sum(ph[ii]=='N')), phS=float(np.sum(ph[ii]=='S')))
    null={cid:{k:[] for k in ['nominal','weight','weight95','scales','jitter','radial','ring','oct','phN','phS']} for cid in G}
    for m in range(nnull):
        if (m+1)%20==0: print(label,'NULL',m+1,'/',nnull,flush=True)
        Xm=shuffled_xyz(ra,de,z,reg,rng); tm=cKDTree(Xm)
        for _,r in cands.iterrows():
            cid=r.candidate_id; g=G[cid]; q,R=g['q'],g['R']; ii=ids(tm,q,R)
            null[cid]['nominal'].append(len(ii)); null[cid]['weight'].append(float(w[ii].sum())); null[cid]['weight95'].append(float(w95[ii].sum()))
            null[cid]['scales'].append([count(tm,q,rr) for rr in SCALE_RADII])
            null[cid]['jitter'].append([count(tm,qq,R) for qq in g['jitter']])
            null[cid]['radial'].append([count(tm,qq,R) for qq in g['radial']])
            null[cid]['ring'].append([count(tm,qq,R) for qq in g['ring']])
            null[cid]['oct'].append(octant_counts(Xm,tm,q,R,g['u'],g['e1'],g['e2']))
            null[cid]['phN'].append(float(np.sum(ph[ii]=='N'))); null[cid]['phS'].append(float(np.sum(ph[ii]=='S')))
        del tm,Xm; gc.collect()
    rows=[]
    for _,r in cands.iterrows():
        cid=r.candidate_id; kind=r.kind; O=obs[cid]; N={k:np.asarray(v,float) for k,v in null[cid].items()}
        def rat(o,n): return float(o/max(np.mean(n),1e-9))
        nom=rat(O['nominal'],N['nominal']); wr=rat(O['weight'],N['weight']); w95r=rat(O['weight95'],N['weight95'])
        scale_mu=np.mean(N['scales'],axis=0); scale_ratio=O['scales']/np.maximum(scale_mu,1)
        jit_mu=np.mean(N['jitter'],axis=0); jit_ratio=O['jitter']/np.maximum(jit_mu,1)
        rad_mu=np.mean(N['radial'],axis=0); rad_ratio=O['radial']/np.maximum(rad_mu,1)
        ring_mu=np.mean(N['ring'],axis=0); ring_ratio=O['ring']/np.maximum(ring_mu,1)
        oct_mu=np.mean(N['oct'],axis=0); oct_ratio=O['oct']/np.maximum(oct_mu,1)
        if kind=='VOID':
            scale_n=int(np.sum(scale_ratio<.90)); jit_frac=float(np.mean(jit_ratio<.90)); oct_n=int(np.sum(oct_ratio<.90));
            radial_rank=int(np.argsort(rad_ratio).tolist().index(2)+1); ring_pass=bool(nom <= np.median(ring_ratio)-.05)
            sign_nom=nom<.85; weight_pass=(wr<.85 and w95r<.85); scale_pass=scale_n>=4; jitter_pass=jit_frac>=.70; oct_pass=oct_n>=5; radial_pass=radial_rank<=3
        else:
            scale_n=int(np.sum(scale_ratio>1.10)); jit_frac=float(np.mean(jit_ratio>1.10)); oct_n=int(np.sum(oct_ratio>1.10));
            radial_rank=int(np.argsort(-rad_ratio).tolist().index(2)+1); ring_pass=bool(nom >= np.median(ring_ratio)+.05)
            sign_nom=nom>1.15; weight_pass=(wr>1.15 and w95r>1.15); scale_pass=scale_n>=4; jitter_pass=jit_frac>=.70; oct_pass=oct_n>=5; radial_pass=radial_rank<=3
        ph_results={}
        ph_pass=True; ph_tested=0
        for P in ['N','S']:
            o=O['ph'+P]; n=N['ph'+P]; mu=float(np.mean(n)); ratio=(o/mu if mu>0 else np.nan)
            supported=mu>=10
            ok=(ratio<.90 if kind=='VOID' else ratio>1.10) if supported else None
            ph_results[P]=(mu,ratio,supported,ok)
            if supported: ph_tested+=1; ph_pass=ph_pass and bool(ok)
        if ph_tested==0: ph_pass=None
        rows.append(dict(candidate_id=cid,kind=kind,dataset=label,nominal_ratio=nom,weighted_ratio=wr,weight95_ratio=w95r,
            nominal_pass=sign_nom,weight_stress_pass=weight_pass,scale_persist_n=scale_n,scale_pass=scale_pass,jitter_fraction=jit_frac,jitter_pass=jitter_pass,
            radial_center_rank=radial_rank,radial_pass=radial_pass,ring_median_ratio=float(np.median(ring_ratio)),ring_pass=ring_pass,octants_same_sign=oct_n,octant_pass=oct_pass,
            photsys_tested=ph_tested,photsys_pass=ph_pass,phN_null_mean=ph_results['N'][0],phN_ratio=ph_results['N'][1],phS_null_mean=ph_results['S'][0],phS_ratio=ph_results['S'][1],
            scale_profile=';'.join(f'{R:.0f}:{x:.3f}' for R,x in zip(SCALE_RADII,scale_ratio)),radial_profile=';'.join(f'{s}:{x:.3f}' for s,x in zip(['-1R','-.5R','0','+.5R','+1R'],rad_ratio))))
    return pd.DataFrame(rows)

def evaluate_eboss(cands,files,zrange,label,nnull,rng):
    ra,de,z,w,ph,reg=load_cat(files,*zrange); X=xyz(ra,de,z); t=cKDTree(X); G=build_geometries(cands)
    obs={cid:count(t,g['q'],g['R']) for cid,g in G.items()}; nn={cid:[] for cid in G}
    for m in range(nnull):
        if (m+1)%50==0: print(label,'NULL',m+1,'/',nnull,flush=True)
        Xm=shuffled_xyz(ra,de,z,reg,rng); tm=cKDTree(Xm)
        for cid,g in G.items(): nn[cid].append(count(tm,g['q'],g['R']))
        del Xm,tm; gc.collect()
    rows=[]
    for _,r in cands.iterrows():
        cid=r.candidate_id; arr=np.asarray(nn[cid],float); mu=float(arr.mean()); sd=float(arr.std(ddof=1)); real=float(obs[cid]); ratio=real/mu if mu>0 else np.nan; zz=(real-mu)/max(sd,1.)
        support=mu>=15 and zrange[0]<=r.z<=zrange[1]
        if r.kind=='VOID': passed=bool(support and ratio<.85 and zz<-2.)
        else: passed=bool(support and ratio>1.15 and zz>2.)
        emp=(1+np.sum(arr<=real))/(len(arr)+1) if r.kind=='VOID' else (1+np.sum(arr>=real))/(len(arr)+1)
        rows.append(dict(candidate_id=cid,kind=r.kind,dataset=label,real=real,null_mean=mu,ratio=ratio,zscore=zz,empirical_p=float(emp),support=support,pass_independent=passed))
    return pd.DataFrame(rows)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--candidates',required=True); ap.add_argument('--v15',required=True); ap.add_argument('--v12',required=True); ap.add_argument('--eboss',required=True); ap.add_argument('--outdir',default='battery'); ap.add_argument('--nulls',type=int,default=100); args=ap.parse_args()
    out=Path(args.outdir); out.mkdir(parents=True,exist_ok=True); rng=np.random.default_rng(SEED); cands=pd.read_csv(args.candidates)
    v15=Path(args.v15);v12=Path(args.v12);eb=Path(args.eboss)
    desi15=evaluate_desi(cands,[v15/'LRG_NGC.fits',v15/'LRG_SGC.fits'],(.4,1.1),args.nulls,'DESI_v1.5_LRG',rng)
    desi12=evaluate_desi(cands,[v12/'LRG_NGC.fits',v12/'LRG_SGC.fits'],(.4,1.1),max(50,args.nulls//2),'DESI_v1.2_LRG',rng)
    # Cross-version retention is defined against nominal density contrast from v1.5.
    merged=desi15[['candidate_id','kind','nominal_ratio']].merge(desi12[['candidate_id','nominal_ratio']],on='candidate_id',suffixes=('_v15','_v12'))
    cv=[]
    for _,r in merged.iterrows():
        if r.kind=='VOID': retention=(1-r.nominal_ratio_v12)/max(1-r.nominal_ratio_v15,1e-6); sign=r.nominal_ratio_v12<.85
        else: retention=(r.nominal_ratio_v12-1)/max(r.nominal_ratio_v15-1,1e-6); sign=r.nominal_ratio_v12>1.15
        cv.append(dict(candidate_id=r.candidate_id,catalog_version_ratio_v15=r.nominal_ratio_v15,catalog_version_ratio_v12=r.nominal_ratio_v12,contrast_retention=float(retention),catalog_version_pass=bool(sign and retention>=.70)))
    cv=pd.DataFrame(cv)
    eframes=[]
    specs=[('eBOSS_LRG',('eBOSS_LRG_NGC.fits','eBOSS_LRG_SGC.fits'),(.6,1.0)),('eBOSS_LRGpCMASS',('eBOSS_LRGpCMASS_NGC.fits','eBOSS_LRGpCMASS_SGC.fits'),(.6,1.0)),('eBOSS_ELG',('eBOSS_ELG_NGC.fits','eBOSS_ELG_SGC.fits'),(.6,1.1)),('eBOSS_QSO',('eBOSS_QSO_NGC.fits','eBOSS_QSO_SGC.fits'),(.8,1.1))]
    for label,(a,b),zr in specs:
        eframes.append(evaluate_eboss(cands,[eb/a,eb/b],zr,label,200,rng))
    eboss=pd.concat(eframes,ignore_index=True)
    desi15.to_csv(out/'desi_v15_stability.csv',index=False); desi12.to_csv(out/'desi_v12_stability.csv',index=False); cv.to_csv(out/'catalog_version.csv',index=False); eboss.to_csv(out/'eboss_independent.csv',index=False)
    summary=[]
    for _,c in cands.iterrows():
        cid=c.candidate_id; d=desi15[desi15.candidate_id==cid].iloc[0]; vv=cv[cv.candidate_id==cid].iloc[0]; ee=eboss[(eboss.candidate_id==cid)&(eboss.support)]
        independent_pass=int(ee.pass_independent.sum()); independent_tests=int(len(ee))
        tests={'nominal':bool(d.nominal_pass),'weights':bool(d.weight_stress_pass),'scales':bool(d.scale_pass),'jitter':bool(d.jitter_pass),'radial':bool(d.radial_pass),'transverse_ring':bool(d.ring_pass),'octants':bool(d.octant_pass),'catalog_v1.2':bool(vv.catalog_version_pass)}
        if d.photsys_pass is not None and not pd.isna(d.photsys_pass): tests['photsys']=bool(d.photsys_pass)
        internal_pass=sum(tests.values()); internal_n=len(tests)
        verdict='SURVIVES_INTERNAL'
        if internal_pass<max(6,internal_n-2): verdict='WEAKENED'
        if independent_tests>0 and independent_pass>0: verdict='INDEPENDENTLY_REPLICATED' if verdict=='SURVIVES_INTERNAL' else 'MIXED'
        summary.append(dict(candidate_id=cid,kind=c.kind,internal_pass=internal_pass,internal_tests=internal_n,independent_pass=independent_pass,independent_tests=independent_tests,verdict=verdict,failed_tests=';'.join(k for k,v in tests.items() if not v)))
    summary=pd.DataFrame(summary); summary.to_csv(out/'battery_summary.csv',index=False)
    meta={'seed':SEED,'desi_nulls':args.nulls,'eboss_nulls':200,'candidates':len(cands),'note':'All candidate coordinates frozen before this battery. Independent surveys not used in discovery.'}; (out/'battery_meta.json').write_text(json.dumps(meta,indent=2))
    print('=== BATTERY SUMMARY ===');print(summary.to_string(index=False));print('\n=== eBOSS ===');print(eboss.to_string(index=False));print('\n=== v1.5 stability ===');print(desi15.to_string(index=False));print('\n=== version ===');print(cv.to_string(index=False))
if __name__=='__main__': main()

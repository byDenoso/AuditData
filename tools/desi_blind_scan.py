#!/usr/bin/env python3
import argparse, json, math, gc
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.integrate import cumulative_trapezoid
from scipy.spatial import cKDTree

H0=67.66; OM=0.3111; OL=1-OM; C=299792.458; SEED=260826
RAD=np.array([50.,75.,100.,130.,160.,200.],float)
zg=np.linspace(0,1.5,30001); inv=1/np.sqrt(OM*(1+zg)**3+OL)
dcg=(C/H0)*cumulative_trapezoid(inv,zg,initial=0)
def dc(z): return np.interp(z,zg,dcg)
def iz(r): return np.interp(r,dcg,zg)

def load(paths,zmin,zmax):
    R=[]
    for reg,p in enumerate(paths):
        print('LOAD',p,flush=True)
        with fits.open(p,memmap=True) as h:
            d=h[1].data; names={n.upper():n for n in d.names}
            ra=np.asarray(d[names['RA']],float); de=np.asarray(d[names['DEC']],float); z=np.asarray(d[names['Z']],float)
            m=np.isfinite(ra)&np.isfinite(de)&np.isfinite(z)&(z>=zmin)&(z<=zmax)
            R.append((ra[m].astype('f4'),de[m].astype('f4'),z[m].astype('f4'),np.full(m.sum(),reg,'i1')))
            print(' kept',int(m.sum()),flush=True)
    return tuple(np.concatenate([x[k] for x in R]) for k in range(4))

def xyz(ra,de,z):
    r=dc(np.asarray(z,float)); a=np.deg2rad(ra); d=np.deg2rad(de); cd=np.cos(d)
    return np.column_stack((r*cd*np.cos(a),r*cd*np.sin(a),r*np.sin(d))).astype('f4')

def rdz(q):
    r=float(np.linalg.norm(q)); x,y,z=q
    return math.degrees(math.atan2(y,x))%360,math.degrees(math.asin(z/r)),float(iz(r)),r

def shuffle_xyz(ra,de,z,reg,rng):
    zz=z.copy()
    for g in np.unique(reg):
        ii=np.where(reg==g)[0]; zz[ii]=z[rng.permutation(ii)]
    return xyz(ra,de,zz)

def queries(ra,de,z,reg,n,rng):
    out=[]
    for g in np.unique(reg):
        ii=np.where(reg==g)[0]; ng=max(1000,round(n*len(ii)/len(reg)))
        ai=rng.choice(ii,ng,replace=True); zi=rng.choice(ii,ng,replace=True)
        out.append(xyz(ra[ai],de[ai],z[zi]))
    q=np.vstack(out)
    if len(q)>n: q=q[rng.choice(len(q),n,replace=False)]
    return q

def counts(tree,q):
    a=np.empty((len(q),len(RAD)),'f4')
    for j,r in enumerate(RAD):
        print(' COUNT',r,flush=True); a[:,j]=tree.query_ball_point(q,r,return_length=True,workers=-1)
    return a

def components(points,link):
    if len(points)==0:return []
    par=np.arange(len(points))
    def root(i):
        while par[i]!=i: par[i]=par[par[i]]; i=par[i]
        return i
    for a,b in cKDTree(points).query_pairs(link):
        a=root(a); b=root(b)
        if a!=b: par[b]=a
    g={}
    for i in range(len(points)):g.setdefault(root(i),[]).append(i)
    return list(g.values())

def pca(pts):
    if len(pts)<20:return [np.nan]*5
    X=pts-np.median(pts,axis=0); w,v=np.linalg.eigh(np.cov(X,rowvar=False)); v=v[:,np.argsort(w)[::-1]]
    P=X@v; e=np.percentile(P,95,axis=0)-np.percentile(P,5,axis=0); a,b,c=e
    return float(a),float(b),float(c),float(b/max(c,1e-6)),float(a/max(b,1e-6))

def detect(q,real,mu,var,tree):
    sig=(real-mu)/np.sqrt(np.maximum(mu+var,1)); rat=real/np.maximum(mu,1)
    R={r:i for i,r in enumerate(RAD)}
    core=[R[75.],R[100.],R[130.],R[160.]]
    pv=(rat[:,core]<.82).sum(1)
    vs=(-np.minimum(sig[:,R[75.]],0)-np.minimum(sig[:,R[100.]],0)-.6*np.minimum(sig[:,R[130.]],0)+5*np.maximum(0,1-rat[:,core].min(1)))
    vm=(pv>=3)&(rat[:,core].min(1)<.68)&(sig[:,R[75.]]<-2.2)&(sig[:,R[100.]]<-2.2)
    vi=np.where(vm)[0]
    if len(vi)>1200:vi=vi[np.argsort(vs[vi])[-1200:]]
    V=[]; outer={75.:130.,100.:160.,130.:200.,160.:200.}
    for gg in components(q[vi],95):
        ids=vi[np.asarray(gg)]; k=int(ids[np.argmax(vs[ids])]); cr=rat[k,core]; rr=float([75,100,130,160][int(np.argmin(cr))]); j=R[rr]; jo=R[outer[rr]]
        shell=(real[k,jo]-real[k,j])/max(mu[k,jo]-mu[k,j],1)
        a,de,z,dist=rdz(q[k]); V.append(dict(ra_deg=a,dec_deg=de,z=z,R_eff_mpc=rr,inner_ratio=float(rat[k,j]),inner_sigma=float(sig[k,j]),shell_ratio=float(shell),persistence=int(pv[k]),seed_cluster_size=len(ids),score=float(vs[k]),comoving_mpc=dist))
    V=pd.DataFrame(V)
    if len(V): V=V.sort_values('score',ascending=False).reset_index(drop=True)

    pw=(rat[:,[R[50.],R[75.],R[100.],R[130.]]]>1.22).sum(1)
    ws=np.maximum(sig[:,R[50.]],0)+np.maximum(sig[:,R[75.]],0)+.5*np.maximum(sig[:,R[100.]],0)+4*np.maximum(0,rat[:,R[75.]]-1)
    wm=(pw>=3)&(rat[:,R[75.]]>1.32)&(sig[:,R[75.]]>2.4)&(sig[:,R[100.]]>1.8)
    wi=np.where(wm)[0]
    if len(wi)>1200:wi=wi[np.argsort(ws[wi])[-1200:]]
    W=[]
    for gg in components(q[wi],105):
        ids=wi[np.asarray(gg)]; k=int(ids[np.argmax(ws[ids])]); cen=q[k].astype(float); ids2=tree.query_ball_point(cen,140); pts=tree.data[np.asarray(ids2)] if len(ids2) else np.empty((0,3)); a1,a2,a3,plan,elong=pca(pts)
        if not np.isfinite(plan) or plan<1.65 or a2<80 or a3>115: continue
        a,de,z,dist=rdz(cen); W.append(dict(ra_deg=a,dec_deg=de,z=z,R_probe_mpc=140.,ratio_75=float(rat[k,R[75.]]),sigma_75=float(sig[k,R[75.]]),sigma_100=float(sig[k,R[100.]]),persistence=int(pw[k]),axis1_mpc=a1,axis2_mpc=a2,axis3_mpc=a3,planarity=plan,elongation=elong,seed_cluster_size=len(ids),score=float(ws[k]*plan),comoving_mpc=dist))
    W=pd.DataFrame(W)
    if len(W): W=W.sort_values('score',ascending=False).reset_index(drop=True)
    return V,W

def validate(cands,paths,kind,rng):
    if len(cands)==0:return cands
    ra,de,z,reg=load(paths,.8,1.1); X=xyz(ra,de,z); t=cKDTree(X); mocks=[]
    for n in range(4):
        print(' ELG NULL',n+1,flush=True); mocks.append(cKDTree(shuffle_xyz(ra,de,z,reg,rng)))
    vals=[]
    for _,r in cands.iterrows():
        if not .8<=r.z<=1.1: vals.append((np.nan,np.nan,False,'outside')); continue
        rad=float(r.get('R_eff_mpc',r.get('R_probe_mpc',100))); d=float(r.comoving_mpc); A=math.radians(r.ra_deg); D=math.radians(r.dec_deg); c=np.array([d*math.cos(D)*math.cos(A),d*math.cos(D)*math.sin(A),d*math.sin(D)])
        nr=t.query_ball_point(c,rad,return_length=True); nn=np.array([x.query_ball_point(c,rad,return_length=True) for x in mocks],float); mu=nn.mean(); va=nn.var(ddof=1)
        if mu<12:vals.append((np.nan,np.nan,False,'low_support'));continue
        ratio=nr/mu; s=(nr-mu)/math.sqrt(max(mu+va,1)); ok=(ratio<.84 and s<-1.4) if kind=='V' else (ratio>1.16 and s>1.4); vals.append((ratio,s,ok,'ok'))
    cands=cands.copy(); cands['elg_ratio']=[x[0] for x in vals]; cands['elg_sigma']=[x[1] for x in vals]; cands['elg_confirm']=[x[2] for x in vals]; cands['elg_status']=[x[3] for x in vals]
    return cands

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--lrg-ngc',required=True);ap.add_argument('--lrg-sgc',required=True);ap.add_argument('--elg-ngc');ap.add_argument('--elg-sgc');ap.add_argument('--queries',type=int,default=50000);ap.add_argument('--mocks',type=int,default=5);ap.add_argument('--outdir',default='outputs');args=ap.parse_args()
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True);rng=np.random.default_rng(SEED)
    ra,de,z,reg=load([args.lrg_ngc,args.lrg_sgc],.4,1.1); X=xyz(ra,de,z); print('LRG',len(X),flush=True); q=queries(ra,de,z,reg,args.queries,rng);tree=cKDTree(X); real=counts(tree,q)
    sm=np.zeros_like(real,float);sq=np.zeros_like(real,float)
    for n in range(args.mocks):
        print('NULL',n+1,'/',args.mocks,flush=True); mt=cKDTree(shuffle_xyz(ra,de,z,reg,rng)); cc=counts(mt,q).astype(float);sm+=cc;sq+=cc*cc;del mt,cc;gc.collect()
    mu=sm/args.mocks;var=np.maximum(sq/args.mocks-mu*mu,0);V,W=detect(q,real,mu,var,tree)
    if args.elg_ngc and args.elg_sgc: V=validate(V,[args.elg_ngc,args.elg_sgc],'V',rng);W=validate(W,[args.elg_ngc,args.elg_sgc],'W',rng)
    if len(V):V.insert(0,'candidate_id',[f'NEXO-V{i+1:02d}' for i in range(len(V))])
    if len(W):W.insert(0,'candidate_id',[f'NEXO-W{i+1:02d}' for i in range(len(W))])
    V.to_csv(out/'void_candidates.csv',index=False);W.to_csv(out/'wall_candidates.csv',index=False)
    summary={'seed':SEED,'H0':H0,'Omega_m':OM,'LRG_objects':int(len(X)),'query_points':int(len(q)),'null_shuffles':args.mocks,'radii_mpc':RAD.tolist(),'void_candidates':int(len(V)),'wall_candidates':int(len(W)),'selection':'blind RA/DEC/z; redshift-shuffle null preserves angular footprint and n(z); no void/wall catalogue used for selection'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True);print('TOP VOIDS',flush=True);print(V.head(20).to_string(index=False) if len(V) else 'none',flush=True);print('TOP WALLS',flush=True);print(W.head(20).to_string(index=False) if len(W) else 'none',flush=True)
if __name__=='__main__':main()

#!/usr/bin/env python3
import argparse, itertools, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import linalg
from scipy.integrate import cumulative_trapezoid

from tools.shoes_host_jackknife import (
    build_precision_cache, fit_delete_cached, get_h0,
    physical_calibrator_key, anchor_constraint_rows
)

C_LIGHT=299792.458
PEER_LOCAL=71.92804067892604
SIG_PEER_LOCAL_BOOK=0.82
TARGET_H0=73.0
RNG_SEED=260827


def delta_mu(h2,h1):
    return -5.0*math.log10(float(h2)/float(h1))


def h_from_delta_mu(h1,dm):
    return float(h1)*10.0**(-float(dm)/5.0)


def subset_rows(keys, subset):
    keys=np.asarray(keys,str)
    wanted=set(subset)
    return np.where(np.array([k in wanted for k in keys],dtype=bool))[0]


def residual_fraction(effect,residual):
    return abs(float(effect))/abs(float(residual))


def Ez(z,om=.3114):
    z=np.asarray(z,float)
    return np.sqrt(om*(1+z)**3+(1-om))


def mu70(z,om=.3114):
    z=np.asarray(z,float)
    zg=np.linspace(0,max(.3,float(np.max(z))*1.02),12000)
    chi=(C_LIGHT/70.0)*cumulative_trapezoid(1/Ez(zg,om),zg,initial=0)
    d=(1+z)*np.interp(z,zg,chi)
    return 5*np.log10(np.maximum(d,1e-12))+25


def load_pantheon_cov(path,n):
    vals=np.loadtxt(path).reshape(-1)
    if len(vals)==n*n+1:
        vals=vals[1:]
    if len(vals)!=n*n:
        raise RuntimeError(f'Pantheon covariance length {len(vals)} != {n*n}')
    return vals.reshape(n,n)


def fit_h0_fixed_om(z,mu,cov,om=.3114):
    cf=linalg.cho_factor(cov+np.eye(len(cov))*1e-12,lower=True,check_finite=False)
    one=np.ones(len(z))
    cinv1=linalg.cho_solve(cf,one,check_finite=False)
    r=mu-mu70(z,om)
    cinvr=linalg.cho_solve(cf,r,check_finite=False)
    a=-float(one@cinvr)/float(one@cinv1)
    sig_a=math.sqrt(1/float(one@cinv1))
    h=70*10**(a/5)
    sig_h=math.log(10)/5*h*sig_a
    return h,sig_h,a


def load_shoes(y_path,L_path,C_path,q_path):
    yd=np.loadtxt(y_path,unpack=True,skiprows=1,dtype={'names':('Source','Data'),'formats':('U64',float)})
    sources=np.asarray(yd[0],str); y=np.asarray(yd[1],float)
    L=np.loadtxt(L_path,delimiter='\t'); C=np.loadtxt(C_path,delimiter='\t'); params=np.loadtxt(q_path,dtype=str).tolist()
    return sources,y,L,C,params


def identifiable_h0(h,s):
    return bool(np.isfinite(h) and np.isfinite(s) and 40.0<h<100.0 and 0.0<s<10.0)


def calibration_posterior_bookkeeping(h_obs,s_obs,n=250000,seed=RNG_SEED):
    rng=np.random.default_rng(seed)
    hloc=rng.normal(PEER_LOCAL,SIG_PEER_LOCAL_BOOK,size=n)
    hdraw=rng.normal(float(h_obs),float(s_obs),size=n)
    good=(hloc>0)&(hdraw>0)
    dm=-5*np.log10(hdraw[good]/hloc[good])
    q=np.quantile(dm,[.025,.16,.5,.84,.975])
    return {
      'delta_mu_p2p5':float(q[0]),'delta_mu_p16':float(q[1]),'delta_mu_median':float(q[2]),
      'delta_mu_p84':float(q[3]),'delta_mu_p97p5':float(q[4]),
      'P_delta_mu_lt_0':float(np.mean(dm<0)),'P_abs_delta_mu_lt_0p032':float(np.mean(np.abs(dm)<0.03212314923877396))
    }


def run_shoes_gate(a,out):
    sources,y,L,C,params=load_shoes(a.y,a.L,a.C,a.q)
    cache=build_precision_cache(y,L,C)
    cov0=linalg.inv(cache['A'],check_finite=False); q0=cov0@cache['g']; h0,s0,_=get_h0(q0,cov0,params)
    if abs(h0-73.04)>.08: raise RuntimeError(f'SH0ES baseline mismatch: {h0}')

    host_params=[p for p in params if p.startswith('mu_')]
    hosts=[p[3:] for p in host_params if p[3:] not in {'N4258','LMC','M31'}]
    keys=np.array([physical_calibrator_key(s,hosts) or '' for s in sources])
    sn_keys=sorted(set(keys)-{''})

    singles=[]
    for key in sn_keys:
        D=subset_rows(keys,[key])
        q,cov,p=fit_delete_cached(y,L,params,cache,D,())
        h,s,_=get_h0(q,cov,p)
        singles.append(dict(calibrator_SN=key,host=key.split('_',1)[0],H0=h,sigma_H0=s,delta_H0=h-h0,equiv_delta_mu=delta_mu(h,h0)))
    S=pd.DataFrame(singles)
    S['abs_delta']=S.delta_H0.abs()
    S.sort_values('delta_H0').to_csv(out/'single_calibrator_influence.csv',index=False)

    residual=TARGET_H0-PEER_LOCAL
    down=S.sort_values('delta_H0').calibrator_SN.tolist()
    absrank=S.sort_values('abs_delta',ascending=False).calibrator_SN.tolist()
    rng=np.random.default_rng(RNG_SEED)
    rows=[]; random_rows=[]
    for k in range(1,6):
        for label,subset in [('most_down',down[:k]),('largest_abs',absrank[:k])]:
            D=subset_rows(keys,subset)
            q,cov,p=fit_delete_cached(y,L,params,cache,D,())
            h,s,_=get_h0(q,cov,p)
            rows.append(dict(strategy=label,k=k,subset=';'.join(subset),H0=h,sigma_H0=s,delta_H0=h-h0,residual_fraction=residual_fraction(h-h0,residual),distance_to_PEERlocal=h-PEER_LOCAL))

        total=math.comb(len(sn_keys),k)
        if k<=3:
            subsets=list(itertools.combinations(sn_keys,k))
            exact=True
        else:
            nrand=min(5000,total)
            seen=set(); subsets=[]
            while len(subsets)<nrand:
                subset=tuple(sorted(rng.choice(sn_keys,size=k,replace=False).tolist()))
                if subset in seen: continue
                seen.add(subset); subsets.append(subset)
            exact=(nrand==total)
        vals=[]
        for subset in subsets:
            D=subset_rows(keys,subset)
            q,cov,p=fit_delete_cached(y,L,params,cache,D,())
            h,s,_=get_h0(q,cov,p)
            vals.append(h-h0)
        vals=np.array(vals)
        nrand=len(vals)
        obs_down=[r for r in rows if r['strategy']=='most_down' and r['k']==k][0]['delta_H0']
        obs_abs=[r for r in rows if r['strategy']=='largest_abs' and r['k']==k][0]['delta_H0']
        random_rows.append(dict(k=k,n_subsets=nrand,exact_enumeration=exact,total_combinations=total,median_delta=float(np.median(vals)),p05=float(np.quantile(vals,.05)),p95=float(np.quantile(vals,.95)),p_as_low_or_lower=float((np.sum(vals<=obs_down)+(0 if exact else 1))/(nrand+(0 if exact else 1))),p_abs_ge_topabs=float((np.sum(np.abs(vals)>=abs(obs_abs))+(0 if exact else 1))/(nrand+(0 if exact else 1))),p_abs_ge_full_residual=float((np.sum(np.abs(vals)>=residual)+(0 if exact else 1))/(nrand+(0 if exact else 1))))
    Cum=pd.DataFrame(rows); Cum.to_csv(out/'cumulative_calibrator_removals.csv',index=False)
    Rnd=pd.DataFrame(random_rows); Rnd.to_csv(out/'random_subset_nulls.csv',index=False)

    anchor_map={
      'N4258':anchor_constraint_rows(sources,'N4258'),
      'LMC':anchor_constraint_rows(sources,'LMC'),
      'MW':np.where(np.isin(sources,['MHW1_Gaia','MHW1_HST']))[0],
    }
    ar=[]
    names=list(anchor_map)
    for k in range(1,4):
        for combo in itertools.combinations(names,k):
            D=np.unique(np.concatenate([anchor_map[x] for x in combo])).astype(int)
            try:
                q,cov,p=fit_delete_cached(y,L,params,cache,D,())
                h,s,_=get_h0(q,cov,p)
                ident=identifiable_h0(h,s)
            except Exception:
                h=s=float('nan'); ident=False
            ar.append(dict(anchors_removed='+'.join(combo),n_rows=len(D),H0=h,sigma_H0=s,identifiable=ident,delta_H0=(h-h0 if ident else np.nan),residual_fraction=(residual_fraction(h-h0,residual) if ident else np.nan)))
    A=pd.DataFrame(ar).sort_values(['identifiable','delta_H0'],ascending=[False,True])
    A.to_csv(out/'anchor_subset_gate.csv',index=False)

    key='N5584_2007af'
    D=subset_rows(keys,[key])
    q,cov,p=fit_delete_cached(y,L,params,cache,D,())
    h_wo,s_wo,_=get_h0(q,cov,p)
    z_book=(h_wo-PEER_LOCAL)/math.sqrt(s_wo**2+SIG_PEER_LOCAL_BOOK**2)

    hier_rows=[]
    scenarios=[('baseline',h0,s0),('without_SN2007af',h_wo,s_wo)]
    for k in [2,3,5]:
        rr=Cum[(Cum.strategy=='most_down')&(Cum.k==k)].iloc[0]
        scenarios.append((f'most_down_k{k}',rr.H0,rr.sigma_H0))
    for name,h,s in scenarios:
        hier_rows.append({'scenario':name,'H0_ladder':h,'sigma_H0':s,**calibration_posterior_bookkeeping(h,s,seed=RNG_SEED+len(hier_rows))})
    H=pd.DataFrame(hier_rows); H.to_csv(out/'hierarchical_calibration_bookkeeping.csv',index=False)

    Aid=A[A.identifiable].copy()
    return dict(
      baseline_H0=h0, baseline_sigma=s0, n_calibrator_SNe=len(sn_keys),
      peer_local=PEER_LOCAL, full_residual=residual, required_delta_mu=delta_mu(TARGET_H0,PEER_LOCAL),
      without_SN2007af=dict(H0=h_wo,sigma_H0=s_wo,residual_to_peer_local=h_wo-PEER_LOCAL,z_bookkeeping=z_book),
      max_single_shift=float(S.abs_delta.max()), max_single_SN=str(S.loc[S.abs_delta.idxmax(),'calibrator_SN']),
      max_identifiable_anchor_subset_shift=float(Aid.delta_H0.abs().max()), max_identifiable_anchor_subset=str(Aid.loc[Aid.delta_H0.abs().idxmax(),'anchors_removed']),
      nonidentifiable_anchor_subsets=A[~A.identifiable].anchors_removed.tolist(),
      random_subset_min_p_full_residual=float(Rnd.p_abs_ge_full_residual.min()),
      hierarchical_bookkeeping=H.to_dict(orient='records')
    )


def run_pv_gate(a,out):
    d=pd.read_csv(a.pantheon,sep=r'\s+')
    C=load_pantheon_cov(a.pantheon_cov,len(d))
    zhd=d.zHD.to_numpy(float); zcmb=d.zCMB.to_numpy(float); mu=d.MU_SH0ES.to_numpy(float)
    cal=d.IS_CALIBRATOR.to_numpy(int)
    rows=[]
    for zmin in [.005,.01,.015,.023,.03,.05,.075]:
        m=(cal==0)&np.isfinite(zhd)&np.isfinite(zcmb)&np.isfinite(mu)&(zhd>=zmin)&(zhd<.15)
        idx=np.where(m)[0]; cov=C[np.ix_(idx,idx)]
        hhd,shd,_=fit_h0_fixed_om(zhd[idx],mu[idx],cov)
        hcmb,scmb,_=fit_h0_fixed_om(zcmb[idx],mu[idx],cov)
        rows.append(dict(zmin=zmin,N=len(idx),H0_zHD=hhd,sigma_zHD=shd,H0_zCMB=hcmb,sigma_zCMB=scmb,delta_zHD_minus_zCMB=hhd-hcmb,delta_fraction_of_peer_local_residual=(hhd-hcmb)/(TARGET_H0-PEER_LOCAL)))
    P=pd.DataFrame(rows)
    P.to_csv(out/'pv_zmin_gate.csv',index=False)
    return dict(max_abs_delta=float(P.delta_zHD_minus_zCMB.abs().max()),zmin_at_max=float(P.loc[P.delta_zHD_minus_zCMB.abs().idxmax(),'zmin']),delta_at_0p023=float(P.loc[np.isclose(P.zmin,.023),'delta_zHD_minus_zCMB'].iloc[0]))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--y',required=True); ap.add_argument('--L',required=True); ap.add_argument('--C',required=True); ap.add_argument('--q',required=True)
    ap.add_argument('--pantheon',required=True); ap.add_argument('--pantheon-cov',required=True)
    ap.add_argument('--out',default='paper_gate')
    a=ap.parse_args(); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    shoes=run_shoes_gate(a,out); pv=run_pv_gate(a,out)
    cum=pd.read_csv(out/'cumulative_calibrator_removals.csv'); rnd=pd.read_csv(out/'random_subset_nulls.csv')
    bestdown=cum[cum.strategy=='most_down'].sort_values('H0').iloc[0]
    verdict={
      'SH0ES':shoes,'PV':pv,
      'paper_gates':{
        'single_calibrator_explains_full_residual': bool(shoes['max_single_shift']>=shoes['full_residual']),
        'any_identifiable_anchor_subset_explains_full_residual': bool(shoes['max_identifiable_anchor_subset_shift']>=shoes['full_residual']),
        'best_downward_k5_reaches_PEERlocal': bool(bestdown.H0<=PEER_LOCAL),
        'best_downward_k5_H0': float(bestdown.H0),
        'random_subset_full_residual_min_p': float(rnd.p_abs_ge_full_residual.min()),
        'survives_without_SN2007af_at_2sigma_bookkeeping': bool(abs(shoes['without_SN2007af']['z_bookkeeping'])<2),
        'pv_effect_material_gt_25pct_residual': bool(pv['max_abs_delta']>.25*shoes['full_residual'])
      },
      'caveat':'Subset tests are structured robustness nulls on the observed ladder, not proof that selected SNe are bad. PEER+Local sigma and hierarchical delta-mu posteriors are bookkeeping sensitivity models, not a joint PEER+CF4+SH0ES posterior.'
    }
    (out/'summary.json').write_text(json.dumps(verdict,indent=2))
    print(json.dumps(verdict,indent=2))
    print('\nCUMULATIVE REMOVALS\n',cum.to_string(index=False))
    print('\nRANDOM SUBSET NULLS\n',rnd.to_string(index=False))
    print('\nANCHORS\n',pd.read_csv(out/'anchor_subset_gate.csv').to_string(index=False))
    print('\nHIERARCHICAL BOOKKEEPING\n',pd.read_csv(out/'hierarchical_calibration_bookkeeping.csv').to_string(index=False))
    print('\nPV ZMIN\n',pd.read_csv(out/'pv_zmin_gate.csv').to_string(index=False))

if __name__=='__main__': main()

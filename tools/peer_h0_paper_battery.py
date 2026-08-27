#!/usr/bin/env python3
import argparse, json, math, itertools
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import linalg
from scipy.optimize import minimize_scalar
from scipy.stats import norm


def required_delta_mu(predicted_h0, target_h0):
    return float(-5.0*np.log10(float(target_h0)/float(predicted_h0)))


def profile_common_offset(observed_h0, sigma_h0, predicted_h0, predicted_sigma=0.0):
    h=np.asarray(observed_h0,float); s=np.asarray(sigma_h0,float)
    w=1.0/s**2
    hm=float(np.sum(w*h)/np.sum(w)); sh=float(np.sqrt(1.0/np.sum(w)))
    dH=hm-float(predicted_h0)
    dmu=required_delta_mu(predicted_h0,hm)
    sig_tot=float(np.hypot(sh,predicted_sigma))
    z=dH/sig_tot if sig_tot>0 else np.nan
    return dict(H0_weighted=hm,sigma_weighted=sh,delta_h0=dH,delta_mu=dmu,z=z,p_two_sided=float(2*norm.sf(abs(z))) if np.isfinite(z) else np.nan)


def leaveout_gate(required_residual, shifts, max_fraction=0.75):
    a=np.abs(np.asarray(shifts,float)); req=abs(float(required_residual))
    frac=float(np.max(a)/req) if req>0 and len(a) else np.nan
    return dict(max_fraction_of_residual=frac,passes=bool(frac<=max_fraction))


def combine_effects(effects, covariance=None):
    e=np.asarray(effects,float)
    if covariance is None:
        return dict(identifiable=False,combined_shift=np.nan,reason='covariance_required')
    C=np.asarray(covariance,float)
    if C.shape!=(len(e),len(e)):
        raise ValueError('covariance shape mismatch')
    return dict(identifiable=True,combined_shift=float(np.sum(e)),sigma=float(np.sqrt(np.sum(C))))


def empirical_pvalue(observed, null, tail='greater'):
    x=np.asarray(null,float)
    if tail=='greater': k=int(np.sum(x>=observed))
    elif tail=='less': k=int(np.sum(x<=observed))
    elif tail=='two-sided': k=int(np.sum(np.abs(x)>=abs(observed)))
    else: raise ValueError('tail')
    return float((k+1)/(len(x)+1))


def load_shoes(y_path,L_path,C_path,q_path):
    yd=np.loadtxt(y_path,unpack=True,skiprows=1,dtype={'names':('Source','Data'),'formats':('U64',float)})
    sources=np.asarray(yd[0],str); y=np.asarray(yd[1],float)
    L=np.loadtxt(L_path,delimiter='\t'); C=np.loadtxt(C_path,delimiter='\t'); params=np.loadtxt(q_path,dtype=str).tolist()
    return sources,y,L,C,params


def shoes_cached_fit(sources,y,L,C,params,drop_sn_keys=()):
    from tools.shoes_host_jackknife import build_precision_cache, fit_delete_cached, get_h0, physical_calibrator_key
    cache=build_precision_cache(y,L,C)
    hosts=[p[3:] for p in params if p.startswith('mu_') and p[3:] not in {'N4258','LMC','M31'}]
    keys=np.array([physical_calibrator_key(s,hosts) or '' for s in sources])
    D=np.where(np.isin(keys,list(drop_sn_keys)))[0] if len(drop_sn_keys) else np.array([],dtype=int)
    q,cov,p=fit_delete_cached(y,L,params,cache,D,())
    h,s,_=get_h0(q,cov,p)
    return h,s,int(len(D))


def load_pantheon(data_path,cov_path):
    from tools.local_void_dark_energy_tests import load_cov
    d=pd.read_csv(data_path,sep=r'\s+')
    C=load_cov(cov_path,len(d))
    return d,C


def fit_pantheon_h0(d,C,zcol,zmin,zmax,om=0.3114):
    from tools.local_void_dark_energy_tests import mu_model
    zall=d[zcol].to_numpy(float); mu=d.MU_SH0ES.to_numpy(float)
    mask=(d.IS_CALIBRATOR.to_numpy(int)==0)&np.isfinite(zall)&np.isfinite(mu)&(zall>=zmin)&(zall<zmax)
    idx=np.where(mask)[0]; z=zall[idx]; y=mu[idx]; cc=C[np.ix_(idx,idx)]
    cf=linalg.cho_factor(cc+np.eye(len(cc))*1e-10,lower=True,check_finite=False)
    def f(H):
        r=y-mu_model(z,float(H),om)
        return float(r@linalg.cho_solve(cf,r,check_finite=False))
    r=minimize_scalar(f,bounds=(60,80),method='bounded',options={'xatol':1e-9})
    H=float(r.x)
    step=0.02; d2=(f(H+step)-2*f(H)+f(H-step))/(step*step)
    sig=math.sqrt(2/d2) if d2>0 else np.nan
    return dict(H0=H,sigma_H0=sig,N=int(len(z)),chi2=float(r.fun),z_definition=zcol,zmin=zmin,zmax=zmax)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--shoes-y',required=True); ap.add_argument('--shoes-L',required=True); ap.add_argument('--shoes-C',required=True); ap.add_argument('--shoes-q',required=True)
    ap.add_argument('--pantheon-data',required=True); ap.add_argument('--pantheon-cov',required=True)
    ap.add_argument('--out',default='peer_h0_paper_battery'); a=ap.parse_args()
    out=Path(a.out); out.mkdir(parents=True,exist_ok=True)

    # Frozen external inputs from prior, independently constructed batteries.
    H_PEER=70.391; SIG_PEER=0.801
    H_LOCAL_MED=71.92804067892604; H_LOCAL_MIN=71.71203503729637; H_LOCAL_MAX=72.17037753571651
    sig_local_env=((H_LOCAL_MAX-H_LOCAL_MIN)/2)/math.sqrt(3)
    sig_peer_scaled=SIG_PEER*(H_LOCAL_MED/H_PEER)
    sig_pred=math.hypot(sig_peer_scaled,sig_local_env)

    sources,y,L,C,params=load_shoes(a.shoes_y,a.shoes_L,a.shoes_C,a.shoes_q)
    base_h,base_s,_=shoes_cached_fit(sources,y,L,C,params,())
    if abs(base_h-73.0434)>0.08: raise RuntimeError(f'SH0ES reproduction failed {base_h}')

    ranked=['N5584_2007af','N1365_2012fr','N3982_1998aq']
    seq=[]
    for k in range(0,4):
        drop=ranked[:k]
        h,s,nrows=shoes_cached_fit(sources,y,L,C,params,drop)
        comb=math.hypot(s,sig_pred); z=(h-H_LOCAL_MED)/comb
        seq.append(dict(k_removed=k,removed=';'.join(drop),H0=h,sigma_H0=s,n_rows_dropped=nrows,residual_vs_PEER_local=h-H_LOCAL_MED,residual_sigma=z,p_two_sided=float(2*norm.sf(abs(z))),equiv_delta_mu=required_delta_mu(H_LOCAL_MED,h)))
    seqdf=pd.DataFrame(seq); seqdf.to_csv(out/'sequential_calibrator_removal.csv',index=False)

    # Pairwise and triple removal sensitivity. These are sensitivity tests, not p-values, because sets were selected after LOO ranking.
    combo=[]
    for k in [1,2,3]:
        for ss in itertools.combinations(ranked,k):
            h,s,nrows=shoes_cached_fit(sources,y,L,C,params,ss)
            combo.append(dict(k=k,removed=';'.join(ss),H0=h,sigma_H0=s,delta_vs_baseline=h-base_h,residual_vs_PEER_local=h-H_LOCAL_MED))
    pd.DataFrame(combo).to_csv(out/'calibrator_combination_sensitivity.csv',index=False)

    d,PC=load_pantheon(a.pantheon_data,a.pantheon_cov)
    pv=[]
    for zmin,zmax in [(0.005,.15),(.01,.15),(.015,.15),(.023,.15),(.03,.15),(.05,.15)]:
        hd=fit_pantheon_h0(d,PC,'zHD',zmin,zmax); cmb=fit_pantheon_h0(d,PC,'zCMB',zmin,zmax)
        pv.append(dict(zmin=zmin,zmax=zmax,N_zHD=hd['N'],H0_zHD=hd['H0'],sigma_zHD=hd['sigma_H0'],H0_zCMB=cmb['H0'],sigma_zCMB=cmb['sigma_H0'],PV_shift_HD_minus_CMB=hd['H0']-cmb['H0']))
    pvdf=pd.DataFrame(pv); pvdf.to_csv(out/'pantheon_pv_zcut_stress.csv',index=False)

    # Structured null: PEER+Local predicts the local ladder with no additional calibration offset.
    null_z=(base_h-H_LOCAL_MED)/math.hypot(base_s,sig_pred)
    null_p=float(2*norm.sf(abs(null_z)))
    dmu=required_delta_mu(H_LOCAL_MED,base_h)
    # Monte Carlo sanity check under the same Gaussian null.
    rng=np.random.default_rng(20260827); sm=math.hypot(base_s,sig_pred)
    mocks=H_LOCAL_MED+rng.normal(0,sm,200000)
    p_mc=empirical_pvalue(abs(base_h-H_LOCAL_MED),np.abs(mocks-H_LOCAL_MED),'greater')

    # Gate from single-SN influence, using exact prior jackknife values from the same matrix reproduction.
    single_shifts=[]
    for sn in ranked:
        h,s,_=shoes_cached_fit(sources,y,L,C,params,[sn]); single_shifts.append(h-base_h)
    gate=leaveout_gate(73.0-H_LOCAL_MED,np.array(single_shifts),0.75)

    summary={
      'frozen_background':{'H0_PEER':H_PEER,'sigma_PEER':SIG_PEER,'H0_PEER_plus_local_median':H_LOCAL_MED,'local_envelope':[H_LOCAL_MIN,H_LOCAL_MAX],'sigma_prediction_bookkeeping':sig_pred},
      'SH0ES_baseline':{'H0':base_h,'sigma':base_s,'delta_mu_required_from_PEER_local':dmu,'residual_sigma':null_z,'p_two_sided':null_p,'p_mc_two_sided_equivalent':p_mc},
      'dominant_calibrator_gate':gate,
      'single_calibrator_shifts':dict(zip(ranked,map(float,single_shifts))),
      'sequential_removal':seq,
      'PV_stress':pv,
      'claims':{
        'nonzero_calibration_offset_required': bool(null_p<0.05),
        'single_calibrator_dominates_over_75pct_residual': bool(not gate['passes']),
        'paper_superposition_hypothesis_survives_without_SN2007af': bool(abs(seq[1]['residual_sigma'])<2.0),
        'PV_is_material_gt_0p2_km_s_Mpc_some_cut': bool(np.max(np.abs(pvdf.PV_shift_HD_minus_CMB))>0.2),
        'effects_additive_without_covariance': False
      },
      'caveats':['Local-Hole contribution is a frozen external result from the prior density/flow battery, not re-fit here.','Sequential top-SN removals are post-selection sensitivity tests and are not assigned discovery p-values.','PV and calibrator-SN shifts are not summed because their covariance is not identified.']
    }
    (out/'summary.json').write_text(json.dumps(summary,indent=2,default=float))
    print(json.dumps(summary,indent=2,default=float))
    print('\nSEQUENTIAL\n',seqdf.to_string(index=False))
    print('\nPV STRESS\n',pvdf.to_string(index=False))

if __name__=='__main__': main()

#!/usr/bin/env python3
import argparse, json, math, re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import linalg


def h0_from_fivelogh0(x):
    return float(10.0**(float(x)/5.0))


def gls_fit(y,L,C):
    y=np.asarray(y,float); L=np.asarray(L,float); C=np.asarray(C,float)
    cf=linalg.cho_factor(C,check_finite=False)
    X=linalg.cho_solve(cf,L,check_finite=False)
    v=linalg.cho_solve(cf,y,check_finite=False)
    A=L.T@X; g=L.T@v
    cov=linalg.inv(A,check_finite=False)
    return cov@g,cov


def drop_host(y,L,C,sources,params,host):
    sources=np.asarray(sources,str); params=list(params)
    d=(sources==host)|np.char.startswith(sources,host+'_')
    keep=~d
    pdrop='mu_'+host
    pkeep=np.array([p!=pdrop for p in params])
    return y[keep],L[keep][:,pkeep],C[np.ix_(keep,keep)],[p for p,k in zip(params,pkeep) if k]


def anchor_constraint_rows(sources,anchor):
    sources=np.asarray(sources,str)
    if anchor=='MW':
        return np.where(np.isin(sources,['MHW1_HST','MHW1_Gaia']))[0]
    return np.where(sources=='mu_'+anchor)[0]


def build_precision_cache(y,L,C):
    cf=linalg.cho_factor(C,check_finite=False)
    P=linalg.cho_solve(cf,np.eye(C.shape[0]),check_finite=False)
    B=P@L
    py=P@y
    A=L.T@B
    g=L.T@py
    return dict(P=P,B=B,py=py,A=A,g=g)


def fit_delete_cached(y,L,params,cache,drop_rows,drop_params=()):
    D=np.array(sorted(set(int(x) for x in drop_rows)),dtype=int)
    if len(D)==0:
        A=cache['A'].copy(); g=cache['g'].copy()
    else:
        P=cache['P']; B=cache['B']; py=cache['py']
        LD=L[D,:]; PD=P[np.ix_(D,D)]; BD=B[D,:]; yD=y[D]; pyD=py[D]
        T=BD-PD@LD
        A0=cache['A']-LD.T@BD-BD.T@LD+LD.T@PD@LD
        g0=cache['g']-T.T@yD-LD.T@pyD
        S=linalg.inv(PD,check_finite=False)
        v=pyD-PD@yD
        A=A0-T.T@S@T
        g=g0-T.T@S@v
    keep=np.array([p not in set(drop_params) for p in params],dtype=bool)
    Ak=A[np.ix_(keep,keep)]; gk=g[keep]
    cov=linalg.inv(Ak,check_finite=False)
    q=cov@gk
    p=[x for x,k in zip(params,keep) if k]
    return q,cov,p


def get_h0(q,cov,params):
    j=params.index('5logH0')
    h=h0_from_fivelogh0(q[j])
    sig=float(math.log(10)/5*h*math.sqrt(max(cov[j,j],0)))
    return h,sig,float(q[j])


def physical_calibrator_key(source,hosts):
    for h in hosts:
        pref=h+'_'
        if source.startswith(pref):
            toks=source.split('_')
            if len(toks)>=3:
                return '_'.join(toks[:2])
    return None


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--y',required=True); ap.add_argument('--L',required=True); ap.add_argument('--C',required=True); ap.add_argument('--q',required=True)
    ap.add_argument('--out',default='shoes_jackknife')
    a=ap.parse_args(); out=Path(a.out); out.mkdir(parents=True,exist_ok=True)

    yd=np.loadtxt(a.y,unpack=True,skiprows=1,dtype={'names':('Source','Data'),'formats':('U64',float)})
    sources=np.asarray(yd[0],str); y=np.asarray(yd[1],float)
    L=np.loadtxt(a.L,delimiter='\t'); C=np.loadtxt(a.C,delimiter='\t'); params=np.loadtxt(a.q,dtype=str).tolist()
    if L.shape!=(len(y),len(params)): raise RuntimeError(f'L shape {L.shape}, expected {(len(y),len(params))}')
    if C.shape!=(len(y),len(y)): raise RuntimeError(f'C shape {C.shape}, expected {(len(y),len(y))}')

    cache=build_precision_cache(y,L,C)
    cov0=linalg.inv(cache['A'],check_finite=False); q0=cov0@cache['g']; h0,s0,f0=get_h0(q0,cov0,params)
    if abs(h0-73.04)>0.08: raise RuntimeError(f'baseline reproduction failed: H0={h0}')

    qa,ca,pa=fit_delete_cached(y,L,params,cache,np.where((sources=='M101')|np.char.startswith(sources,'M101_'))[0],['mu_M101'])
    yy,LL,CC,pp=drop_host(y,L,C,sources,params,'M101')
    qd,cd=gls_fit(yy,LL,CC)
    if np.max(np.abs(qa-qd))>2e-6: raise RuntimeError('cached deletion does not reproduce direct GLS')

    host_params=[p for p in params if p.startswith('mu_')]
    special={'N4258','LMC','M31'}
    hosts=[p[3:] for p in host_params if p[3:] not in special]

    host_rows=[]
    for host in hosts:
        D=np.where((sources==host)|np.char.startswith(sources,host+'_'))[0]
        q,cov,p=fit_delete_cached(y,L,params,cache,D,['mu_'+host])
        h,s,_=get_h0(q,cov,p)
        host_rows.append(dict(host=host,n_rows_dropped=len(D),n_ceph=int(np.sum(sources==host)),n_cal_rows=int(np.sum(np.char.startswith(sources,host+'_'))),H0=h,sigma_H0=s,delta_H0=h-h0,delta_sigma=(h-h0)/s0,equiv_delta_mu=-5*np.log10(h/h0)))
    H=pd.DataFrame(host_rows).sort_values('delta_H0')
    H.to_csv(out/'leave_one_host_out.csv',index=False)

    keys=np.array([physical_calibrator_key(s,hosts) or '' for s in sources])
    sn_rows=[]
    for key in sorted(set(keys)-{''}):
        D=np.where(keys==key)[0]
        q,cov,p=fit_delete_cached(y,L,params,cache,D,())
        h,s,_=get_h0(q,cov,p)
        host=key.split('_',1)[0]
        sn_rows.append(dict(calibrator_SN=key,host=host,n_rows_dropped=len(D),H0=h,sigma_H0=s,delta_H0=h-h0,delta_sigma=(h-h0)/s0,equiv_delta_mu=-5*np.log10(h/h0)))
    S=pd.DataFrame(sn_rows).sort_values('delta_H0')
    S.to_csv(out/'leave_one_calibrator_sn_out.csv',index=False)

    anchor_rows=[]
    for anchor in ['N4258','LMC','MW']:
        D=anchor_constraint_rows(sources,anchor)
        q,cov,p=fit_delete_cached(y,L,params,cache,D,())
        h,s,_=get_h0(q,cov,p)
        anchor_rows.append(dict(anchor=anchor,n_constraint_rows_dropped=len(D),constraint_sources=';'.join(sorted(set(sources[D]))),H0=h,sigma_H0=s,delta_H0=h-h0,delta_sigma=(h-h0)/s0,equiv_delta_mu=-5*np.log10(h/h0)))
    A=pd.DataFrame(anchor_rows)
    A.to_csv(out/'leave_one_anchor_constraint_out.csv',index=False)

    D=np.where(sources=='M31')[0]
    qm,cm,pm=fit_delete_cached(y,L,params,cache,D,['mu_M31'])
    hm,sm,_=get_h0(qm,cm,pm)

    Habs=H.assign(abs_delta=lambda x:abs(x.delta_H0)).sort_values('abs_delta',ascending=False)
    Sabs=S.assign(abs_delta=lambda x:abs(x.delta_H0)).sort_values('abs_delta',ascending=False)
    required_from_peer_local=73.0-71.92804067892604
    ext_mask=np.ones(len(sources),dtype=bool)
    for h in hosts+['N4258','LMC','M31']:
        ext_mask &= (sources!=h)
        ext_mask &= ~np.char.startswith(sources,h+'_')
    ext_counts=pd.Series(sources[ext_mask]).value_counts().head(30)
    ext_counts.rename_axis('source').reset_index(name='count').to_csv(out/'external_source_labels_top30.csv',index=False)

    from tools.h0_superposition_kill_tests import run_superposition_battery
    superposition=run_superposition_battery(y,L,C,params,sources,cache,hosts,S,out,nmocks=50000)

    summary={
      'data_source':'marcushogas/Cepheid-Distance-Ladder-Data SH0ES2022; modified text representation of Riess et al. 2022 matrices',
      'baseline':{'H0':h0,'sigma_H0':s0,'five_log_H0':f0,'n_obs':len(y),'n_params':len(params)},
      'host_jackknife':{'n_hosts':len(H),'H0_min':float(H.H0.min()),'H0_max':float(H.H0.max()),'max_abs_shift':float(Habs.abs_delta.iloc[0]),'max_host':str(Habs.host.iloc[0]),'top3_abs_shift_sum':float(Habs.abs_delta.head(3).sum()),'n_abs_gt_0p2':int((Habs.abs_delta>0.2).sum()),'n_abs_gt_0p5':int((Habs.abs_delta>0.5).sum()),'top8':Habs.head(8).to_dict(orient='records')},
      'calibrator_sn_jackknife':{'n_physical_sne':len(S),'H0_min':float(S.H0.min()),'H0_max':float(S.H0.max()),'max_abs_shift':float(Sabs.abs_delta.iloc[0]),'max_SN':str(Sabs.calibrator_SN.iloc[0]),'top5_abs_shift_sum':float(Sabs.abs_delta.head(5).sum()),'n_abs_gt_0p2':int((Sabs.abs_delta>0.2).sum()),'n_abs_gt_0p5':int((Sabs.abs_delta>0.5).sum()),'top10':Sabs.head(10).to_dict(orient='records')},
      'anchor_constraint_jackknife':A.to_dict(orient='records'),
      'M31_nonSN_host_control':{'n_rows_dropped':int(len(D)),'H0':hm,'sigma_H0':sm,'delta_H0':hm-h0},
      'peer_local_residual_to_73':required_from_peer_local,
      'superposition_kill_tests':superposition,
      'interpretation_rules':['Large individual leave-one-out shifts identify influential data, not automatically bad data.','Host and SN jackknives are correlated; their shifts must not be summed.','Anchor-constraint deletion removes only the geometric prior rows and retains the Cepheid data.','PEER and local-environment inputs are frozen external values in the superposition battery; they are not refit to SH0ES.']
    }
    (out/'summary.json').write_text(json.dumps(summary,indent=2,default=float))
    print(json.dumps(summary,indent=2,default=float)); print('\nTOP HOST INFLUENCE\n',Habs.head(12).to_string(index=False)); print('\nTOP SN INFLUENCE\n',Sabs.head(15).to_string(index=False)); print('\nANCHORS\n',A.to_string(index=False))

if __name__=='__main__': main()

#!/usr/bin/env python3
import io, itertools, math, re, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import linalg

PANTHEON_URL='https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat'
C_LIGHT=299792.458
Q0=-0.55
J0=1.0
H_PEER=70.391
H_PEER_LOCAL=71.92804067892604


def build_precision_cache(y,L,C):
    from tools.shoes_host_jackknife import build_precision_cache as f
    return f(y,L,C)


def fit_delete_cached(y,L,params,cache,drop_rows,drop_params=()):
    from tools.shoes_host_jackknife import fit_delete_cached as f
    return f(y,L,params,cache,drop_rows,drop_params)


def canonical_cid(name):
    s=str(name).strip()
    # SH0ES matrix labels append the survey/photometry ID, e.g. 2000dn_57 or 2008ar_65.
    s=re.sub(r'_\d+$','',s,flags=re.I)
    return re.sub(r'[^a-z0-9]','',s.lower())


def h0free_dl_shape(z,q0=Q0,j0=J0):
    z=np.asarray(z,float)
    return C_LIGHT*z*(1.0 + 0.5*(1.0-q0)*z - (1.0-q0-3.0*q0*q0+j0)*z*z/6.0)


def redshift_y_shift(z_old,z_new):
    a=float(h0free_dl_shape(z_old)); b=float(h0free_dl_shape(z_new))
    if a<=0 or b<=0: return float('nan')
    return -5.0*math.log10(b/a)


def anchor_rows(sources,anchor):
    s=np.asarray(sources,str)
    if anchor=='MW': return np.where(np.isin(s,['MHW1_HST','MHW1_Gaia']))[0]
    return np.where(s=='mu_'+anchor)[0]


def _delete_A_K(L,params,cache,drop_rows,drop_params=()):
    D=np.array(sorted(set(int(x) for x in drop_rows)),dtype=int)
    if len(D)==0:
        A=cache['A'].copy(); K=cache['B'].T.copy()
    else:
        P=cache['P']; B=cache['B']
        LD=L[D,:]; PD=P[np.ix_(D,D)]; BD=B[D,:]
        T=BD-PD@LD
        S=linalg.inv(PD,check_finite=False)
        A0=cache['A']-LD.T@BD-BD.T@LD+LD.T@PD@LD
        A=A0-T.T@S@T
        K=B.T-(LD.T+T.T@S)@P[D,:]
    drop=set(drop_params)
    keep=np.array([p not in drop for p in params],dtype=bool)
    return A[np.ix_(keep,keep)],K[keep,:],[p for p,k in zip(params,keep) if k]


def deleted_estimator_weight(L,params,cache,drop_rows,drop_params=(),target='5logH0'):
    A,K,p=_delete_A_K(L,params,cache,drop_rows,drop_params)
    cov=linalg.inv(A,check_finite=False)
    j=p.index(target)
    return cov[j,:]@K


def _fit_new_y(ynew,params,cache):
    cov=linalg.inv(cache['A'],check_finite=False)
    q=cov@(cache['B'].T@np.asarray(ynew,float))
    j=params.index('5logH0')
    f=float(q[j]); h=10.0**(f/5.0)
    sig=math.log(10)/5*h*math.sqrt(max(float(cov[j,j]),0.0))
    return h,sig,f


def _load_pantheon():
    with urllib.request.urlopen(PANTHEON_URL,timeout=60) as r:
        raw=r.read()
    p=pd.read_csv(io.BytesIO(raw),sep=r'\s+')
    p['canon']=p['CID'].astype(str).map(canonical_cid)
    return p


def _hf_crossmatch(sources,L,params,pantheon):
    j=params.index('5logH0')
    hf=np.where(np.abs(L[:,j])>1e-12)[0]
    pp=pantheon[pantheon['USED_IN_SH0ES_HF']==1].copy()
    lookup={}
    for c,g in pp.groupby('canon'):
        lookup[c]={k:float(np.nanmedian(g[k])) for k in ['zHD','zCMB','zHEL']}
    rows=[]
    for i in hf:
        c=canonical_cid(sources[i])
        if c in lookup:
            rows.append({'row':int(i),'source':sources[i],'canon':c,**lookup[c]})
    return np.asarray(hf,int),pd.DataFrame(rows)


def _physical_keys(sources,hosts):
    keys=[]
    for source in sources:
        key=''
        for h in hosts:
            pref=h+'_'
            if str(source).startswith(pref):
                toks=str(source).split('_')
                if len(toks)>=3: key='_'.join(toks[:2])
                break
        keys.append(key)
    return np.asarray(keys,str)


def _h0_from_fit(q,cov,params):
    j=params.index('5logH0'); h=10.0**(float(q[j])/5.0)
    s=math.log(10)/5*h*math.sqrt(max(float(cov[j,j]),0.0))
    return h,s


def _mock_max_influence(W,C,observed_shifts,nmocks,seed=20260827):
    V=W@C@W.T
    V=(V+V.T)/2
    vals,vecs=np.linalg.eigh(V); vals=np.clip(vals,0,None)
    R=vecs@np.diag(np.sqrt(vals))
    rng=np.random.default_rng(seed)
    obs=np.asarray(observed_shifts,float)
    obsmax=float(np.max(np.abs(obs))); obsdown=float(np.min(obs))
    ge=0; le=0; chunk=5000; done=0
    while done<nmocks:
        m=min(chunk,nmocks-done)
        z=rng.normal(size=(m,W.shape[0]))
        t=z@R.T
        h=73.0*10.0**(t/5.0)
        d=h[:,1:]-h[:,[0]]
        ge+=int(np.sum(np.max(np.abs(d),axis=1)>=obsmax))
        le+=int(np.sum(np.min(d,axis=1)<=obsdown))
        done+=m
    return {'observed_max_abs_delta_H0':obsmax,'observed_most_negative_delta_H0':obsdown,
            'p_max_abs':(ge+1)/(nmocks+1),'p_most_negative':(le+1)/(nmocks+1),'n_mocks':nmocks}


def run_superposition_battery(y,L,C,params,sources,cache,hosts,sn_table,out,nmocks=50000):
    out=Path(out)
    p=_load_pantheon()
    hf_all,hf=_hf_crossmatch(sources,L,params,p)
    coverage=len(hf)/max(len(hf_all),1)
    if coverage<0.80:
        matched=set(hf['row'].tolist()) if len(hf) else set()
        unmatched=[str(sources[i]) for i in hf_all if int(i) not in matched][:30]
        raise RuntimeError(f'Pantheon crossmatch coverage too low: {len(hf)}/{len(hf_all)}={coverage:.3f}; examples={unmatched}')

    cov0=linalg.inv(cache['A'],check_finite=False); q0=cov0@cache['g']
    h0,s0=_h0_from_fit(q0,cov0,params)
    gap=h0-H_PEER_LOCAL

    zrows=[]
    for cut in [0.023,0.03,0.04,0.05,0.06,0.075,0.10]:
        D=hf.loc[hf.zHD<cut,'row'].to_numpy(int)
        q,cov,pp=fit_delete_cached(y,L,params,cache,D,())
        h,s=_h0_from_fit(q,cov,pp)
        zrows.append({'zmin':cut,'n_hf_rows_dropped':len(D),'H0':h,'sigma_H0':s,'delta_H0':h-h0})
    Z=pd.DataFrame(zrows); Z.to_csv(out/'zmin_ladder.csv',index=False)

    fr=[]
    for newcol in ['zCMB','zHEL']:
        yy=np.array(y,float,copy=True); shifts=[]
        for r in hf.itertuples(index=False):
            ds=redshift_y_shift(r.zHD,getattr(r,newcol))
            if np.isfinite(ds): yy[int(r.row)]+=ds; shifts.append(ds)
        h,s,_=_fit_new_y(yy,params,cache)
        fr.append({'frame':'zHD->'+newcol,'n_rows':len(shifts),'median_delta_y_mag':float(np.median(shifts)),
                   'H0':h,'sigma_H0':s,'delta_H0':h-h0})
    F=pd.DataFrame(fr); F.to_csv(out/'redshift_frame_swap.csv',index=False)

    keys=_physical_keys(sources,hosts)
    sorted_down=sn_table.sort_values('delta_H0').calibrator_SN.tolist()
    rng=np.random.default_rng(20260827)
    uniq=sorted(set(keys)-{''})
    target=[]
    for k in [1,2,3,5,10]:
        chosen=sorted_down[:k]
        D=np.where(np.isin(keys,chosen))[0]
        q,cov,pp=fit_delete_cached(y,L,params,cache,D,())
        h,s=_h0_from_fit(q,cov,pp); td=h-h0
        draws=[]
        for _ in range(1000):
            rand=rng.choice(uniq,size=k,replace=False)
            Dr=np.where(np.isin(keys,rand))[0]
            qr,cr,pr=fit_delete_cached(y,L,params,cache,Dr,())
            hr,_=_h0_from_fit(qr,cr,pr); draws.append(hr-h0)
        draws=np.asarray(draws)
        target.append({'k':k,'removed':';'.join(chosen),'H0':h,'sigma_H0':s,'delta_H0':td,
                       'fraction_gap_closed':min(max((-td)/gap,0),2),'random_subset_p_one_sided':float((np.sum(draws<=td)+1)/(len(draws)+1)),
                       'random_delta_p05':float(np.quantile(draws,.05)),'random_delta_median':float(np.median(draws))})
    T=pd.DataFrame(target); T.to_csv(out/'targeted_calibrator_removals.csv',index=False)

    ar=[]; anchors=['N4258','LMC','MW']
    for k in [1,2,3]:
        for comb in itertools.combinations(anchors,k):
            D=np.unique(np.concatenate([anchor_rows(sources,a) for a in comb])).astype(int)
            try:
                q,cov,pp=fit_delete_cached(y,L,params,cache,D,())
                h,s=_h0_from_fit(q,cov,pp); status='ok'
            except Exception as e:
                h=s=float('nan'); status=type(e).__name__
            ar.append({'anchors_removed':'+'.join(comb),'n_rows':len(D),'H0':h,'sigma_H0':s,'delta_H0':h-h0 if np.isfinite(h) else np.nan,'status':status})
    A=pd.DataFrame(ar); A.to_csv(out/'anchor_subset_tests.csv',index=False)

    basew=deleted_estimator_weight(L,params,cache,[],(),target='5logH0')
    sn_ws=[basew]; sn_obs=[]
    for key in uniq:
        D=np.where(keys==key)[0]
        sn_ws.append(deleted_estimator_weight(L,params,cache,D,(),target='5logH0'))
        rr=sn_table[sn_table.calibrator_SN==key]
        sn_obs.append(float(rr.delta_H0.iloc[0]))
    snnull=_mock_max_influence(np.vstack(sn_ws),C,sn_obs,nmocks,seed=20260827)

    host_ws=[basew]; host_obs=[]
    for host in hosts:
        D=np.where((sources==host)|np.char.startswith(sources,host+'_'))[0]
        host_ws.append(deleted_estimator_weight(L,params,cache,D,['mu_'+host],target='5logH0'))
        q,cov,pp=fit_delete_cached(y,L,params,cache,D,['mu_'+host])
        hh,_=_h0_from_fit(q,cov,pp); host_obs.append(hh-h0)
    hostnull=_mock_max_influence(np.vstack(host_ws),C,host_obs,nmocks,seed=20260828)
    N=pd.DataFrame([{'family':'calibrator_SN',**snnull},{'family':'host',**hostnull}])
    N.to_csv(out/'mock_influence_nulls.csv',index=False)

    summary={
      'H0_SH0ES_matrix':h0,'sigma_SH0ES_matrix':s0,'H0_PEER_frozen':H_PEER,'H0_PEER_plus_local_frozen':H_PEER_LOCAL,
      'SH0ES_minus_PEER_local':gap,'residual_73_minus_PEER_local':73.0-H_PEER_LOCAL,
      'hubble_flow_matrix_rows':int(len(hf_all)),'pantheon_crossmatched_rows':int(len(hf)),'crossmatch_fraction':coverage,
      'zmin_H0_range':[float(Z.H0.min()),float(Z.H0.max())],
      'frame_swap':F.to_dict(orient='records'),
      'targeted_calibrator_removals':T.to_dict(orient='records'),
      'anchor_subsets':A.to_dict(orient='records'),
      'mock_nulls':N.to_dict(orient='records'),
      'caveats':['PV frame swaps hold the published SH0ES covariance fixed and are sensitivity tests, not replacement likelihoods.',
                 'Targeted removal uses data-ranked calibrators; random-subset p values quantify concentration, not evidence that removed SNe are bad data.',
                 'Parametric mock null assumes the published linear model and covariance are the data-generating process.',
                 'PEER and Local-Hole values are frozen external inputs; this battery does not refit PEER or the density field.']}
    pd.Series(summary).to_json(out/'superposition_summary.json',indent=2)
    return summary

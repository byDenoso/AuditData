#!/usr/bin/env python3
import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

C = 299792.458
OM = 0.33
CMB_RA = 167.94
CMB_DEC = -6.94
CMB_V = 369.82
DEFAULT_H0 = 70.391
SCAN_D = [25.0, 50.0, 75.0, 100.0, 125.0, 150.0, 175.0, 200.0]
SCAN_R = [30.0, 50.0, 75.0, 100.0, 125.0, 150.0, 180.0, 220.0]
PROFILE_R = [25.0, 40.0, 60.0, 80.0, 100.0, 125.0, 150.0, 180.0, 220.0, 260.0]
ICRS_TO_GAL = np.array([
    [-0.0548755604, -0.8734370902, -0.4838350155],
    [ 0.4941094279, -0.4448296300,  0.7469822445],
    [-0.8676661490, -0.1980763734,  0.4559837762],
])


def E(z):
    z = np.asarray(z, float)
    return np.sqrt(OM * (1.0 + z) ** 3 + (1.0 - OM))


def chi_of_z(z, H0=DEFAULT_H0):
    z = np.asarray(z, float)
    zmax = max(0.2, float(np.nanmax(z)) * 1.03)
    zg = np.linspace(0.0, zmax, 30000)
    chi = (C / H0) * cumulative_trapezoid(1.0 / E(zg), zg, initial=0.0)
    return np.interp(z, zg, chi)


def eq_unitvec(ra, dec):
    a = np.deg2rad(np.asarray(ra, float))
    d = np.deg2rad(np.asarray(dec, float))
    q = np.cos(d)
    return np.c_[q * np.cos(a), q * np.sin(a), np.sin(d)]


def gal_unitvec(lon, lat):
    l = np.deg2rad(np.asarray(lon, float))
    b = np.deg2rad(np.asarray(lat, float))
    q = np.cos(b)
    return np.c_[q * np.cos(l), q * np.sin(l), np.sin(b)]


def xyz_to_lb(v):
    v = np.asarray(v, float)
    r = np.linalg.norm(v, axis=-1)
    u = v / np.maximum(r[..., None], 1e-15)
    l = np.rad2deg(np.arctan2(u[..., 1], u[..., 0])) % 360.0
    b = np.rad2deg(np.arcsin(np.clip(u[..., 2], -1.0, 1.0)))
    return l, b


def fibonacci_sphere(n):
    n = int(n)
    i = np.arange(n, dtype=float)
    z = 1.0 - 2.0 * (i + 0.5) / n
    phi = (2.0 * math.pi * i / ((1.0 + math.sqrt(5.0)) / 2.0)) % (2.0 * math.pi)
    r = np.sqrt(np.maximum(0.0, 1.0 - z * z))
    return np.c_[r * np.cos(phi), r * np.sin(phi), z]


def shuffle_positions(pos, rng):
    pos = np.asarray(pos, float)
    rad = np.linalg.norm(pos, axis=1)
    u = pos / np.maximum(rad[:, None], 1e-15)
    return u[rng.permutation(len(u))] * rad[:, None]


def scan_counts_against_shuffles(pos, centers, radius, nnull=40, seed=1):
    pos = np.asarray(pos, float)
    centers = np.asarray(centers, float)
    obs_tree = cKDTree(pos)
    obs = obs_tree.query_ball_point(centers, float(radius), return_length=True).astype(float)
    null = np.empty((int(nnull), len(centers)), dtype=float)
    rng = np.random.default_rng(seed)
    for j in range(int(nnull)):
        sp = shuffle_positions(pos, rng)
        tr = cKDTree(sp)
        null[j] = tr.query_ball_point(centers, float(radius), return_length=True)
    mean = null.mean(axis=0)
    med = np.median(null, axis=0)
    sd = null.std(axis=0, ddof=1)
    ratio = obs / np.maximum(med, 1.0)
    z = (obs - mean) / np.maximum(sd, 1e-9)
    p = (1.0 + np.sum(null <= obs[None, :], axis=0)) / (len(null) + 1.0)
    return dict(obs=obs, null_mean=mean, null_median=med, null_std=sd, ratio=ratio, z=z, p_lower=p, null_counts=null)


def outflow_regression(pos, vpec, center):
    pos = np.asarray(pos, float)
    vpec = np.asarray(vpec, float)
    center = np.asarray(center, float)
    nhat = pos / np.maximum(np.linalg.norm(pos, axis=1)[:, None], 1e-12)
    x = np.sum((pos - center[None, :]) * nhat, axis=1)
    m = np.isfinite(x) & np.isfinite(vpec)
    x = x[m]
    y = vpec[m]
    if len(x) < 4:
        return dict(n=int(len(x)), slope=np.nan, intercept=np.nan, slope_se=np.nan, spearman_r=np.nan, spearman_p=np.nan)
    X = np.c_[np.ones(len(x)), x]
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    resid = y - X @ beta
    dof = max(1, len(x) - 2)
    s2 = float(np.sum(resid * resid) / dof)
    cov = s2 * np.linalg.inv(X.T @ X)
    se = float(np.sqrt(max(cov[1, 1], 0.0)))
    rho, rp = spearmanr(x, y)
    return dict(n=int(len(x)), slope=float(beta[1]), intercept=float(beta[0]), slope_se=se, spearman_r=float(rho), spearman_p=float(rp))


def read_2mrs(path):
    colspecs = [(0,16),(17,26),(27,36),(37,46),(47,56),(57,63),(78,84),(173,178)]
    names = ['ID','RAdeg','DEdeg','GLON','GLAT','Kcmag','Ktmag','cz']
    d = pd.read_fwf(path, colspecs=colspecs, names=names, na_values=['','     '])
    for c in names[1:]:
        d[c] = pd.to_numeric(d[c], errors='coerce')
    return d


def cmb_corrected_cz(df):
    eq = eq_unitvec(df.RAdeg.to_numpy(float), df.DEdeg.to_numpy(float))
    apex = eq_unitvec([CMB_RA], [CMB_DEC])[0]
    return df.cz.to_numpy(float) + CMB_V * (eq @ apex)


def make_2mrs_positions(df, frame='cmb', maglim=11.75, H0=DEFAULT_H0):
    m = np.isfinite(df.cz) & np.isfinite(df.GLON) & np.isfinite(df.GLAT) & np.isfinite(df.Kcmag) & (df.Kcmag <= maglim)
    q = df.loc[m].copy()
    if frame == 'cmb':
        cz = cmb_corrected_cz(q)
    else:
        cz = q.cz.to_numpy(float)
    z = np.clip(cz / C, 1e-5, None)
    r = chi_of_z(z, H0)
    u = gal_unitvec(q.GLON, q.GLAT)
    pos = u * r[:, None]
    keep = (r >= 5.0) & (r <= 450.0)
    return q.loc[keep].reset_index(drop=True), pos[keep], r[keep]


def build_centers(ndir=384):
    u = fibonacci_sphere(ndir)
    rows = []
    centers = []
    for D in SCAN_D:
        for j, v in enumerate(u):
            c = D * v
            l, b = xyz_to_lb(c[None, :])
            rows.append(dict(center_id=f'D{int(D):03d}_U{j:04d}', D_Mpc=D, l_deg=float(l[0]), b_deg=float(b[0]), x=float(c[0]), y=float(c[1]), z=float(c[2])))
            centers.append(c)
    return pd.DataFrame(rows), np.asarray(centers)


def full_scan(pos, centers_df, centers, radii=SCAN_R, nnull=80, seed=26082631):
    tree = cKDTree(pos)
    obs = np.vstack([tree.query_ball_point(centers, R, return_length=True) for R in radii]).T.astype(float)
    null = np.empty((nnull, len(centers), len(radii)), dtype=np.float32)
    rng = np.random.default_rng(seed)
    for j in range(nnull):
        sp = shuffle_positions(pos, rng)
        tr = cKDTree(sp)
        for k, R in enumerate(radii):
            null[j, :, k] = tr.query_ball_point(centers, R, return_length=True)
    mean = null.mean(axis=0)
    med = np.median(null, axis=0)
    sd = null.std(axis=0, ddof=1)
    z = (obs - mean) / np.maximum(sd, 1e-6)
    ratio = obs / np.maximum(med, 1.0)
    p = (1.0 + np.sum(null <= obs[None, :, :], axis=0)) / (nnull + 1.0)
    rows = []
    for i, c in centers_df.iterrows():
        persistence = int(np.sum((ratio[i] < 0.85) & (z[i] < -2.0) & (med[i] >= 20)))
        for k, R in enumerate(radii):
            rows.append({**c.to_dict(), 'R_Mpc':float(R), 'N_obs':int(obs[i,k]), 'N_null_mean':float(mean[i,k]), 'N_null_median':float(med[i,k]), 'N_null_sd':float(sd[i,k]), 'ratio':float(ratio[i,k]), 'z_local':float(z[i,k]), 'p_local':float(p[i,k]), 'persistence_scales':persistence})
    grid = pd.DataFrame(rows)
    # Global look-elsewhere based on the same standardized candidate grid.
    nullz = (null - mean[None, :, :]) / np.maximum(sd[None, :, :], 1e-6)
    mask = med >= 20
    data_min = float(np.nanmin(np.where(mask, z, np.nan)))
    null_min = np.nanmin(np.where(mask[None, :, :], nullz, np.nan), axis=(1,2))
    p_global = float((1.0 + np.sum(null_min <= data_min)) / (nnull + 1.0))
    return grid, null, dict(data_min_z=data_min, null_min_z_median=float(np.median(null_min)), null_min_z_p05=float(np.quantile(null_min,.05)), p_global_min_z=p_global, nnull=int(nnull), ncandidate_cells=int(mask.sum()))


def select_components(grid, maxn=12):
    q = grid[(grid.N_null_median >= 20) & (grid.ratio < 0.82) & (grid.z_local < -2.5) & (grid.persistence_scales >= 2)].copy()
    if q.empty:
        q = grid[(grid.N_null_median >= 20)].nsmallest(60, 'z_local').copy()
    q['rank_score'] = q.z_local - 0.6 * q.persistence_scales
    q = q.sort_values(['rank_score','ratio'])
    keep = []
    for _, r in q.iterrows():
        c = np.array([r.x, r.y, r.z], float)
        duplicate = False
        for k in keep:
            ck = np.array([k['x'], k['y'], k['z']], float)
            sep = np.linalg.norm(c - ck)
            if sep < max(45.0, 0.55 * min(float(r.R_Mpc), float(k['R_Mpc']))):
                duplicate = True
                break
        if not duplicate:
            keep.append(r.to_dict())
        if len(keep) >= maxn:
            break
    out = pd.DataFrame(keep)
    if not out.empty:
        out.insert(0, 'candidate', [f'LNV{i+1:02d}' for i in range(len(out))])
        out['observer_inside'] = out.D_Mpc < out.R_Mpc
    return out


def profile_candidate(pos, center, nnull=80, seed=1):
    center = np.asarray(center, float)[None, :]
    tree = cKDTree(pos)
    obs = np.array([tree.query_ball_point(center, R, return_length=True)[0] for R in PROFILE_R], float)
    null = np.empty((nnull, len(PROFILE_R)), float)
    rng = np.random.default_rng(seed)
    for j in range(nnull):
        tr = cKDTree(shuffle_positions(pos, rng))
        null[j] = [tr.query_ball_point(center, R, return_length=True)[0] for R in PROFILE_R]
    med = np.median(null, axis=0)
    mean = null.mean(axis=0)
    sd = null.std(axis=0, ddof=1)
    ratio = obs / np.maximum(med, 1.0)
    z = (obs - mean) / np.maximum(sd, 1e-6)
    shell_obs = np.diff(np.r_[0, obs])
    shell_med = np.diff(np.r_[0, med])
    shell_ratio = shell_obs / np.maximum(shell_med, 1.0)
    return pd.DataFrame(dict(R_Mpc=PROFILE_R, N_obs=obs.astype(int), N_null_median=med, ratio_cumulative=ratio, z_cumulative=z, shell_ratio=shell_ratio))


def robustness_for_candidate(df2mrs, row, nnull=36):
    center = np.array([row.x, row.y, row.z], float)[None, :]
    R = float(row.R_Mpc)
    records = []
    seed0 = 26082670 + int(str(row.candidate).replace('LNV',''))
    for frame in ['cmb','helio']:
        for mag in [11.25, 11.50, 11.75]:
            _, pos, _ = make_2mrs_positions(df2mrs, frame=frame, maglim=mag)
            s = scan_counts_against_shuffles(pos, center, R, nnull=nnull, seed=seed0 + int(100*mag) + (0 if frame=='cmb' else 1000))
            records.append(dict(candidate=row.candidate, test='frame_mag', frame=frame, maglim=mag, R_Mpc=R, ratio=float(s['ratio'][0]), z=float(s['z'][0]), p=float(s['p_lower'][0]), null_median=float(s['null_median'][0])))
    # 3x3x3 center perturbations around the selected center, full sample, CMB frame.
    _, pos, _ = make_2mrs_positions(df2mrs, frame='cmb', maglim=11.75)
    offs = np.array([[a,b,c] for a in [-25.,0.,25.] for b in [-25.,0.,25.] for c in [-25.,0.,25.]])
    centers = center + offs
    s = scan_counts_against_shuffles(pos, centers, R, nnull=nnull, seed=seed0+9999)
    records.append(dict(candidate=row.candidate, test='center_grid_summary', frame='cmb', maglim=11.75, R_Mpc=R, ratio=float(s['ratio'][13]), z=float(s['z'][13]), p=float(s['p_lower'][13]), null_median=float(s['null_median'][13]), center_grid_fraction_ratio_lt_0p85=float(np.mean(s['ratio']<.85)), center_grid_ratio_min=float(np.min(s['ratio'])), center_grid_ratio_max=float(np.max(s['ratio']))))
    return pd.DataFrame(records)


def read_cf4_groups(path):
    d = pd.read_csv(path, sep='\t', comment='#', engine='python')
    # VizieR may append units/meta lines; force numeric and discard non-data rows.
    needed = ['RAJ2000','DEJ2000','GLON','GLAT','Dist','Vpec']
    for c in needed + ['e_DMzp','Ngal']:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors='coerce')
    d = d.dropna(subset=[c for c in needed if c in d.columns]).copy()
    if 'e_DMzp' in d.columns:
        d = d[(d.e_DMzp <= 1.0) | (~np.isfinite(d.e_DMzp))]
    d = d[(d.Dist >= 5) & (d.Dist <= 350) & np.isfinite(d.Vpec) & (np.abs(d.Vpec) <= 6000)]
    u = gal_unitvec(d.GLON, d.GLAT)
    d['x'] = u[:,0] * d.Dist.to_numpy(float)
    d['y'] = u[:,1] * d.Dist.to_numpy(float)
    d['z'] = u[:,2] * d.Dist.to_numpy(float)
    return d.reset_index(drop=True)


def cf4_test_candidate(cf4, row, nperm=300):
    center = np.array([row.x, row.y, row.z], float)
    pos = cf4[['x','y','z']].to_numpy(float)
    rr = np.linalg.norm(pos-center[None,:], axis=1)
    testR = max(100.0, 1.5*float(row.R_Mpc))
    m = rr <= testR
    q = cf4.loc[m].copy()
    p = pos[m]
    fit = outflow_regression(p, q.Vpec.to_numpy(float), center)
    if fit['n'] < 30:
        return {**fit, 'candidate':row.candidate, 'test_R_Mpc':testR, 'perm_p_upper':np.nan, 'perm_z':np.nan, 'null_slope_mean':np.nan, 'null_slope_sd':np.nan}
    rng = np.random.default_rng(26082700 + int(str(row.candidate).replace('LNV','')))
    dist = q.Dist.to_numpy(float)
    v = q.Vpec.to_numpy(float)
    bins = np.floor(dist / 25.0).astype(int)
    null = []
    for _ in range(nperm):
        vv = v.copy()
        for b in np.unique(bins):
            ix = np.flatnonzero(bins==b)
            vv[ix] = vv[rng.permutation(ix)]
        null.append(outflow_regression(p, vv, center)['slope'])
    null = np.asarray(null, float)
    p_up = float((1 + np.sum(null >= fit['slope'])) / (len(null)+1))
    z = float((fit['slope'] - np.nanmean(null)) / max(np.nanstd(null,ddof=1),1e-9))
    return {**fit, 'candidate':row.candidate, 'test_R_Mpc':testR, 'perm_p_upper':p_up, 'perm_z':z, 'null_slope_mean':float(np.nanmean(null)), 'null_slope_sd':float(np.nanstd(null,ddof=1))}


def mu_model(z, H0):
    d = (1+np.asarray(z,float))*chi_of_z(z,H0)
    return 5*np.log10(np.maximum(d,1e-12))+25


def enclosed_top_hat(r, R, delta, comp=2.0):
    r = np.asarray(r,float)
    out = np.zeros_like(r)
    q1 = r <= R
    out[q1] = delta
    R2 = comp*R
    d2 = -delta*R**3/(R2**3-R**3)
    q2 = (r>R)&(r<=R2)
    out[q2] = (delta*R**3 + d2*(r[q2]**3-R**3))/np.maximum(r[q2]**3,1e-12)
    return out


def void_velocity(relpos, H0, R, delta):
    r = np.linalg.norm(relpos,axis=-1)
    db = enclosed_top_hat(r,R,delta,2.0)
    f = OM**0.55
    return -(H0*f/3.0)*db[...,None]*relpos


def pantheon_incremental_h0(pan, row, bgal=1.1):
    m = (pan.IS_CALIBRATOR==0)&(pan.zHD>=.023)&(pan.zHD<.15)&np.isfinite(pan.RA)&np.isfinite(pan.DEC)
    q = pan.loc[m]
    eq = eq_unitvec(q.RA,q.DEC)
    ug = (ICRS_TO_GAL @ eq.T).T
    delta_g = float(row.ratio)-1.0
    delta_m = delta_g/bgal
    R = float(row.R_Mpc)
    center = np.array([row.x,row.y,row.z],float)
    out=[]
    for H in [70.391,71.62]:
        z=q.zHD.to_numpy(float); rr=chi_of_z(z,H); src=ug*rr[:,None]
        vo=void_velocity((-center)[None,:],H,R,delta_m)[0]
        vs=void_velocity(src-center[None,:],H,R,delta_m)
        vlos=np.sum((vs-vo[None,:])*ug,axis=1)
        zapp=z+(1+z)*vlos/C
        mut=mu_model(z,H); mub=mu_model(zapp,H)
        c=float(np.mean(mut-mub)); hf=H*10**(-c/5.0)
        out.append(dict(candidate=row.candidate,H0_background=H,H0_incremental=float(hf-H),H0_app=float(hf),boost_percent=float(100*(hf/H-1)),delta_g_localized=delta_g,delta_m_assumed=delta_m,R_Mpc=R,observer_offset_Mpc=float(row.D_Mpc)))
    return out


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--2mrs',dest='mrs',required=True)
    ap.add_argument('--cf4',required=True)
    ap.add_argument('--pantheon',required=True)
    ap.add_argument('--out',default='nearby_scan')
    ap.add_argument('--ndir',type=int,default=384)
    ap.add_argument('--nnull',type=int,default=80)
    args=ap.parse_args()
    out=Path(args.out); out.mkdir(parents=True,exist_ok=True)

    d2=read_2mrs(args.mrs)
    _,pos,_=make_2mrs_positions(d2,frame='cmb',maglim=11.75)
    centers_df, centers=build_centers(args.ndir)
    grid,null,global_stats=full_scan(pos,centers_df,centers,nnull=args.nnull)
    grid.to_csv(out/'scan_grid.csv',index=False)
    top=select_components(grid,12)
    top.to_csv(out/'top_candidates.csv',index=False)

    prof=[]; robust=[]
    for _,r in top.iterrows():
        p=profile_candidate(pos,[r.x,r.y,r.z],nnull=args.nnull,seed=26082640+int(r.candidate.replace('LNV','')))
        p.insert(0,'candidate',r.candidate); prof.append(p)
        robust.append(robustness_for_candidate(d2,r,nnull=max(24,args.nnull//2)))
    profiles=pd.concat(prof,ignore_index=True) if prof else pd.DataFrame()
    robustness=pd.concat(robust,ignore_index=True) if robust else pd.DataFrame()
    profiles.to_csv(out/'candidate_profiles.csv',index=False)
    robustness.to_csv(out/'robustness_tests.csv',index=False)

    cf4=read_cf4_groups(args.cf4)
    cfrows=[cf4_test_candidate(cf4,r) for _,r in top.iterrows()]
    cft=pd.DataFrame(cfrows)
    cft.to_csv(out/'cf4_flow_tests.csv',index=False)

    pan=pd.read_csv(args.pantheon,sep=r'\s+')
    hrows=[]
    for _,r in top.iterrows(): hrows.extend(pantheon_incremental_h0(pan,r))
    ht=pd.DataFrame(hrows); ht.to_csv(out/'incremental_h0_tests.csv',index=False)

    # Join compact verdict per candidate.
    verdict=top.copy()
    if not cft.empty:
        verdict=verdict.merge(cft[['candidate','n','slope','slope_se','spearman_r','perm_p_upper','perm_z']],on='candidate',how='left')
    rb=[]
    if not robustness.empty:
        for c,g in robustness.groupby('candidate'):
            fm=g[g.test=='frame_mag']
            cg=g[g.test=='center_grid_summary']
            rb.append(dict(candidate=c,robust_ratio_min=float(fm.ratio.min()),robust_ratio_max=float(fm.ratio.max()),robust_all_ratio_lt_0p9=bool((fm.ratio<.9).all()),center_grid_fraction_ratio_lt_0p85=float(cg.center_grid_fraction_ratio_lt_0p85.iloc[0]) if len(cg) else np.nan))
    if rb: verdict=verdict.merge(pd.DataFrame(rb),on='candidate',how='left')
    verdict['density_gate']=(verdict.ratio<.82)&(verdict.persistence_scales>=2)
    if 'perm_p_upper' in verdict:
        verdict['cf4_outflow_gate']=(verdict.perm_p_upper<.05)&(verdict.slope>0)
    else:
        verdict['cf4_outflow_gate']=False
    verdict['multimessenger_like_gate']=verdict.density_gate & verdict.cf4_outflow_gate
    verdict.to_csv(out/'candidate_verdicts.csv',index=False)

    # Best profile morphology summaries.
    morph=[]
    for c,g in profiles.groupby('candidate') if not profiles.empty else []:
        core=g[g.ratio_cumulative<.85]
        coreR=float(core.R_Mpc.max()) if len(core) else np.nan
        after=g[g.R_Mpc>coreR] if np.isfinite(coreR) else g.iloc[0:0]
        cross=after[after.ratio_cumulative>=1.0]
        compR=float(cross.R_Mpc.iloc[0]) if len(cross) else np.nan
        morph.append(dict(candidate=c,core_R_lt0p85=coreR,compensation_cross_R=compR,min_cumulative_ratio=float(g.ratio_cumulative.min()),max_shell_ratio=float(g.shell_ratio.max())))
    morph=pd.DataFrame(morph); morph.to_csv(out/'morphology_summary.csv',index=False)

    summary={
        'scope':'blind local-universe scan around Milky Way; candidate centers 25-200 Mpc, radii 30-220 Mpc',
        'detector':'2MRS spatial counts vs angular-permutation nulls that preserve exact observer-centric radial distribution and angular mask directions',
        'important_blind_spot':'observer-centered isotropic radial underdensity is intentionally preserved by the null and therefore is not rediscovered; use published Local Hole radial profile for that component',
        '2MRS_rows':int(len(d2)),
        '2MRS_used':int(len(pos)),
        'CF4_groups_used':int(len(cf4)),
        'global_look_elsewhere':global_stats,
        'n_components':int(len(top)),
        'n_density_gate':int(verdict.density_gate.sum()) if len(verdict) else 0,
        'n_cf4_outflow_gate':int(verdict.cf4_outflow_gate.sum()) if len(verdict) else 0,
        'n_joint_gate':int(verdict.multimessenger_like_gate.sum()) if len(verdict) else 0,
        'top_candidates':verdict.head(12).to_dict(orient='records'),
        'published_centered_Local_Hole_baseline':{'delta_g_inner_100_hinv_Mpc':-0.20,'delta_g_cumulative_150_hinv_Mpc':-0.13,'role':'separate radial baseline, not selected by this angular-localization scan'},
        'PEER_rule':'PEER/H0 is post-selection only. incremental_h0_tests.csv estimates only the anisotropic component above the observer-centric radial distribution.'
    }
    (out/'summary.json').write_text(json.dumps(summary,indent=2,default=float))
    print(json.dumps(summary,indent=2,default=float))
    print('\nTOP CANDIDATES\n',verdict.to_string(index=False))
    print('\nCF4 FLOW\n',cft.to_string(index=False))
    print('\nH0 INCREMENTAL\n',ht.to_string(index=False))


if __name__=='__main__':
    main()

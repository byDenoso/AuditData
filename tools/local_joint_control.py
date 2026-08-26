#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import rankdata

from tools.local_universe_nearby_scan import (
    SCAN_D, SCAN_R, build_centers, make_2mrs_positions, read_2mrs,
    read_cf4_groups,
)


def tail_ranks(x, lower=True):
    x = np.asarray(x, float)
    n = len(x)
    r = rankdata(x, method='average')
    if lower:
        return r / (n + 1.0)
    return (n + 1.0 - r) / (n + 1.0)


def joint_score(p_density, p_flow):
    pdn = np.asarray(p_density, float)
    pfl = np.asarray(p_flow, float)
    return -np.log10(np.maximum(pdn * pfl, 1e-300))


def random_rotation(rng):
    q = rng.normal(size=4)
    q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array([
        [1-2*(y*y+z*z), 2*(x*y-z*w), 2*(x*z+y*w)],
        [2*(x*y+z*w), 1-2*(x*x+z*z), 2*(y*z-x*w)],
        [2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x*x+y*y)],
    ])


def fast_flow_slope(pos, vpec, center, idx):
    if len(idx) < 30:
        return np.nan, len(idx)
    p = pos[idx]
    y = vpec[idx]
    nh = p / np.maximum(np.linalg.norm(p, axis=1)[:, None], 1e-12)
    x = np.sum((p - center[None, :]) * nh, axis=1)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]; y = y[m]
    if len(x) < 30:
        return np.nan, len(x)
    dx = x - x.mean(); dy = y - y.mean()
    den = np.sum(dx*dx)
    if den <= 0:
        return np.nan, len(x)
    return float(np.sum(dx*dy)/den), len(x)


def compute_maps(mrs_path, cf4_path, ndir=384):
    d2 = read_2mrs(mrs_path)
    _, pos2, _ = make_2mrs_positions(d2, frame='cmb', maglim=11.75)
    centers_df, centers = build_centers(ndir)
    dirs = centers[:ndir] / SCAN_D[0]

    t2 = cKDTree(pos2)
    density = {}
    for D in SCAN_D:
        ids = np.flatnonzero(np.isclose(centers_df.D_Mpc.to_numpy(float), D))
        cc = centers[ids]
        for R in SCAN_R:
            density[(D,R)] = t2.query_ball_point(cc, R, return_length=True).astype(float)

    cf4 = read_cf4_groups(cf4_path)
    posv = cf4[['x','y','z']].to_numpy(float)
    vv = cf4.Vpec.to_numpy(float)
    tv = cKDTree(posv)
    flow = {}
    nflow = {}
    testR_for_R = {R:max(100.0,1.5*R) for R in SCAN_R}
    for D in SCAN_D:
        ids = np.flatnonzero(np.isclose(centers_df.D_Mpc.to_numpy(float), D))
        cc = centers[ids]
        for R in SCAN_R:
            tr = testR_for_R[R]
            neighborhoods = tv.query_ball_point(cc, tr)
            slopes = np.empty(ndir, float)
            ns = np.empty(ndir, int)
            for j, idx in enumerate(neighborhoods):
                slopes[j], ns[j] = fast_flow_slope(posv, vv, cc[j], np.asarray(idx, int))
            flow[(D,R)] = slopes
            nflow[(D,R)] = ns
    return dirs, density, flow, nflow, len(pos2), len(cf4)


def build_cell_table(dirs, density, flow, nflow):
    rows=[]
    # directions are in Galactic Cartesian coords
    l = np.rad2deg(np.arctan2(dirs[:,1],dirs[:,0])) % 360
    b = np.rad2deg(np.arcsin(np.clip(dirs[:,2],-1,1)))
    for D in SCAN_D:
        for R in SCAN_R:
            cnt = density[(D,R)]
            sl = flow[(D,R)]
            ns = nflow[(D,R)]
            pdn = tail_ranks(cnt, lower=True)
            # NaNs get neutral/worst p=1
            good=np.isfinite(sl)
            pf=np.ones(len(sl),float)
            if good.sum()>2:
                pf[good]=tail_ranks(sl[good],lower=False)
            score=joint_score(pdn,pf)
            score[(~good)|(sl<=0)] = 0.0
            for j in range(len(dirs)):
                rows.append(dict(D_Mpc=D,R_Mpc=R,dir_id=j,l_deg=float(l[j]),b_deg=float(b[j]),N_2MRS=int(cnt[j]),p_density_rank=float(pdn[j]),CF4_slope=float(sl[j]) if np.isfinite(sl[j]) else np.nan,CF4_n=int(ns[j]),p_flow_rank=float(pf[j]),joint_score=float(score[j])))
    return pd.DataFrame(rows)


def rotation_global_null(dirs, table, nrot=300, seed=26082699):
    d_tree=cKDTree(dirs)
    # reshape keyed arrays: group order D then R, direction fastest
    groups=[(D,R) for D in SCAN_D for R in SCAN_R]
    pd_maps=[]; pf_maps=[]; sl_maps=[]
    for D,R in groups:
        g=table[(table.D_Mpc==D)&(table.R_Mpc==R)].sort_values('dir_id')
        pd_maps.append(g.p_density_rank.to_numpy(float))
        pf_maps.append(g.p_flow_rank.to_numpy(float))
        sl_maps.append(g.CF4_slope.to_numpy(float))
    pd_maps=np.asarray(pd_maps); pf_maps=np.asarray(pf_maps); sl_maps=np.asarray(sl_maps)
    data_score=joint_score(pd_maps,pf_maps)
    data_score[(~np.isfinite(sl_maps))|(sl_maps<=0)]=0
    data_max=float(np.nanmax(data_score))
    rng=np.random.default_rng(seed)
    mx=np.empty(nrot,float)
    for k in range(nrot):
        Rm=random_rotation(rng)
        rotated=dirs @ Rm.T
        _, mapping=d_tree.query(rotated,k=1)
        pf=pf_maps[:,mapping]
        sl=sl_maps[:,mapping]
        s=joint_score(pd_maps,pf)
        s[(~np.isfinite(sl))|(sl<=0)]=0
        mx[k]=np.nanmax(s)
    pglob=float((1+np.sum(mx>=data_max))/(nrot+1))
    return dict(data_max_joint_score=data_max,null_max_median=float(np.median(mx)),null_max_p95=float(np.quantile(mx,.95)),null_max_p99=float(np.quantile(mx,.99)),p_global_rotation=pglob,nrot=nrot),mx


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--2mrs',dest='mrs',required=True)
    ap.add_argument('--cf4',required=True)
    ap.add_argument('--out',default='joint_control')
    ap.add_argument('--ndir',type=int,default=384)
    ap.add_argument('--nrot',type=int,default=300)
    a=ap.parse_args()
    out=Path(a.out); out.mkdir(parents=True,exist_ok=True)
    dirs,density,flow,nflow,n2,ncf=compute_maps(a.mrs,a.cf4,a.ndir)
    tab=build_cell_table(dirs,density,flow,nflow)
    tab.to_csv(out/'joint_cell_table.csv',index=False)
    stats,mx=rotation_global_null(dirs,tab,a.nrot)
    pd.DataFrame({'null_max_joint_score':mx}).to_csv(out/'rotation_null_maxima.csv',index=False)
    top=tab.sort_values('joint_score',ascending=False).head(30).copy()
    top.to_csv(out/'top_joint_cells.csv',index=False)

    # Locate prior LNV06 coordinates approximately in this same grid.
    target=np.array([45.55122344050117,-18.052953038073817])
    l=np.deg2rad(tab.l_deg.to_numpy(float)); b=np.deg2rad(tab.b_deg.to_numpy(float))
    lt,bt=np.deg2rad(target)
    ang=np.arccos(np.clip(np.sin(b)*np.sin(bt)+np.cos(b)*np.cos(bt)*np.cos(l-lt),-1,1))
    q=tab[(tab.D_Mpc==50)&(tab.R_Mpc==50)].copy()
    lq=np.deg2rad(q.l_deg.to_numpy(float)); bq=np.deg2rad(q.b_deg.to_numpy(float))
    aq=np.arccos(np.clip(np.sin(bq)*np.sin(bt)+np.cos(bq)*np.cos(bt)*np.cos(lq-lt),-1,1))
    prior=q.iloc[int(np.argmin(aq))].to_dict()

    summary={'method':'real-sky rank calibration within each D,R plus random SO(3) relative rotations of CF4 flow maps against 2MRS density maps','2MRS_used':n2,'CF4_groups_used':ncf,'global_rotation':stats,'prior_LNV06_same_cell':prior,'top10':top.head(10).to_dict(orient='records'),'interpretation_rule':'This tests spatial association of low density with positive peculiar-velocity outflow while preserving the observed structure of each map. It is not a LambdaCDM cosmic-rarity calibration.'}
    (out/'summary.json').write_text(json.dumps(summary,indent=2,default=float))
    print(json.dumps(summary,indent=2,default=float))
    print('\nTOP JOINT\n',top.head(20).to_string(index=False))

if __name__=='__main__':
    main()

#!/usr/bin/env python3
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import linalg

from tools.shoes_host_jackknife import (
    anchor_constraint_rows,
    build_precision_cache,
    fit_delete_cached,
    get_h0,
    physical_calibrator_key,
)

PEER_H0 = 70.391
PEER_SIGMA = 0.801
LOCAL_MEDIAN = 71.92804067892604
LOCAL_LO = 71.71203503729637
LOCAL_HI = 72.17037753571651
TARGET_H0 = 73.0
RNG_SEED = 271828


def equiv_dmu(h2, h1):
    return float(-5.0 * math.log10(float(h2) / float(h1)))


def combined_local_sigma(peer_h0=PEER_H0, peer_sigma=PEER_SIGMA,
                         local_median=LOCAL_MEDIAN, local_lo=LOCAL_LO,
                         local_hi=LOCAL_HI):
    # PEER uncertainty propagates multiplicatively through the frozen observer boost.
    sigma_peer_local = float(peer_sigma) * float(local_median) / float(peer_h0)
    # The Local-Hole range is not a posterior. Use a uniform-envelope bookkeeping sigma.
    sigma_hole = ((float(local_hi) - float(local_lo)) / 2.0) / math.sqrt(3.0)
    return float(math.sqrt(sigma_peer_local**2 + sigma_hole**2))


def rank_signed_influence(rows, direction='down'):
    if direction not in {'down', 'up'}:
        raise ValueError('direction must be down or up')
    reverse = direction == 'up'
    return sorted(rows, key=lambda r: float(r['delta_H0']), reverse=reverse)


def load_shoes(y_path, L_path, C_path, q_path):
    yd = np.loadtxt(y_path, unpack=True, skiprows=1,
                    dtype={'names': ('Source', 'Data'), 'formats': ('U64', float)})
    sources = np.asarray(yd[0], str)
    y = np.asarray(yd[1], float)
    L = np.loadtxt(L_path, delimiter='\t')
    C = np.loadtxt(C_path, delimiter='\t')
    params = np.loadtxt(q_path, dtype=str).tolist()
    if L.shape != (len(y), len(params)):
        raise RuntimeError(f'L shape {L.shape}, expected {(len(y), len(params))}')
    if C.shape != (len(y), len(y)):
        raise RuntimeError(f'C shape {C.shape}, expected {(len(y), len(y))}')
    return y, L, C, params, sources


def baseline_fit(y, L, C, params):
    cache = build_precision_cache(y, L, C)
    cov = linalg.inv(cache['A'], check_finite=False)
    q = cov @ cache['g']
    h0, sig, five = get_h0(q, cov, params)
    if abs(h0 - 73.04) > 0.08:
        raise RuntimeError(f'SH0ES baseline reproduction failed: {h0}')
    return cache, h0, sig, five


def physical_sn_groups(sources, hosts):
    keys = np.array([physical_calibrator_key(s, hosts) or '' for s in sources])
    groups = {}
    for key in sorted(set(keys) - {''}):
        groups[key] = np.where(keys == key)[0]
    return groups


def fit_drop_rows(y, L, params, cache, rows):
    q, cov, p = fit_delete_cached(y, L, params, cache, rows, ())
    return get_h0(q, cov, p)[:2]


def single_sn_influence(y, L, params, cache, groups, h0_base, sig_base):
    rows = []
    for key, D in groups.items():
        h, s = fit_drop_rows(y, L, params, cache, D)
        rows.append({
            'calibrator_SN': key,
            'host': key.split('_', 1)[0],
            'n_rows_dropped': len(D),
            'H0': h,
            'sigma_H0': s,
            'delta_H0': h - h0_base,
            'delta_sigma_baseline': (h - h0_base) / sig_base,
            'equiv_delta_mu': equiv_dmu(h, h0_base),
        })
    return pd.DataFrame(rows).sort_values('delta_H0').reset_index(drop=True)


def host_influence(y, L, params, sources, cache, h0_base, sig_base):
    special = {'N4258', 'LMC', 'M31'}
    hosts = [p[3:] for p in params if p.startswith('mu_') and p[3:] not in special]
    out = []
    for host in hosts:
        D = np.where((sources == host) | np.char.startswith(sources, host + '_'))[0]
        q, cov, p = fit_delete_cached(y, L, params, cache, D, ['mu_' + host])
        h, s, _ = get_h0(q, cov, p)
        out.append({
            'host': host,
            'n_rows_dropped': len(D),
            'n_ceph': int(np.sum(sources == host)),
            'n_cal_rows': int(np.sum(np.char.startswith(sources, host + '_'))),
            'H0': h,
            'sigma_H0': s,
            'delta_H0': h - h0_base,
            'delta_sigma_baseline': (h - h0_base) / sig_base,
            'equiv_delta_mu': equiv_dmu(h, h0_base),
        })
    return pd.DataFrame(out).sort_values('delta_H0').reset_index(drop=True), hosts


def anchor_influence(y, L, params, sources, cache, h0_base, sig_base):
    specs = [
        ('N4258', anchor_constraint_rows(sources, 'N4258')),
        ('LMC', anchor_constraint_rows(sources, 'LMC')),
        ('MW', np.where(np.isin(sources, ['MHW1_Gaia', 'MHW1_HST']))[0]),
    ]
    out = []
    for name, D in specs:
        h, s = fit_drop_rows(y, L, params, cache, D)
        out.append({
            'anchor': name,
            'n_constraint_rows_dropped': len(D),
            'constraint_sources': ';'.join(sorted(set(sources[D]))),
            'H0': h,
            'sigma_H0': s,
            'delta_H0': h - h0_base,
            'delta_sigma_baseline': (h - h0_base) / sig_base,
            'equiv_delta_mu': equiv_dmu(h, h0_base),
        })
    return pd.DataFrame(out)


def greedy_sequence(y, L, params, cache, groups, max_k=5, direction='down'):
    chosen = []
    remaining = set(groups)
    rows = []
    dropped = np.array([], dtype=int)
    for k in range(1, max_k + 1):
        best = None
        for key in sorted(remaining):
            D = np.unique(np.concatenate([dropped, groups[key]])).astype(int)
            h, s = fit_drop_rows(y, L, params, cache, D)
            objective = h if direction == 'down' else -h
            candidate = (objective, key, h, s, D)
            if best is None or candidate[0] < best[0]:
                best = candidate
        _, key, h, s, D = best
        chosen.append(key)
        remaining.remove(key)
        dropped = D
        rows.append({
            'direction': direction,
            'k': k,
            'added_SN': key,
            'subset': ';'.join(chosen),
            'n_rows_dropped': len(dropped),
            'H0': h,
            'sigma_H0': s,
        })
    return pd.DataFrame(rows)


def random_subset_nulls(y, L, params, cache, groups, greedy_df, n_draws=1000, seed=RNG_SEED):
    rng = np.random.default_rng(seed)
    keys = np.array(sorted(groups), dtype=object)
    out = []
    for direction in ['down', 'up']:
        gd = greedy_df[greedy_df.direction == direction]
        for _, grow in gd.iterrows():
            k = int(grow.k)
            observed = float(grow.H0)
            vals = np.empty(n_draws, float)
            for i in range(n_draws):
                subset = rng.choice(keys, size=k, replace=False)
                D = np.unique(np.concatenate([groups[str(s)] for s in subset])).astype(int)
                vals[i] = fit_drop_rows(y, L, params, cache, D)[0]
            if direction == 'down':
                empirical = (1 + int(np.sum(vals <= observed))) / (n_draws + 1)
            else:
                empirical = (1 + int(np.sum(vals >= observed))) / (n_draws + 1)
            out.append({
                'direction': direction,
                'k': k,
                'greedy_H0': observed,
                'random_mean_H0': float(np.mean(vals)),
                'random_std_H0': float(np.std(vals, ddof=1)),
                'random_p05_H0': float(np.quantile(vals, 0.05)),
                'random_p50_H0': float(np.quantile(vals, 0.50)),
                'random_p95_H0': float(np.quantile(vals, 0.95)),
                'empirical_tail_fraction': empirical,
                'n_draws': n_draws,
                'note': 'Sensitivity percentile only; greedy subset was optimized after seeing the data and this is not a discovery p-value.',
            })
    return pd.DataFrame(out)


def build_paper_gate(h0_base, sig_base, local_sigma, sn_df, host_df, anchor_df,
                     greedy_df, random_df):
    resid = h0_base - LOCAL_MEDIAN
    resid_sig = math.sqrt(sig_base**2 + local_sigma**2)
    z = resid / resid_sig
    max_sn = float(sn_df.delta_H0.abs().max())
    max_host = float(host_df.delta_H0.abs().max())
    max_anchor = float(anchor_df.delta_H0.abs().max())
    down = greedy_df[greedy_df.direction == 'down'].copy()
    enters_local = down[(down.H0 >= LOCAL_LO) & (down.H0 <= LOCAL_HI)]
    first_k = int(enters_local.k.min()) if len(enters_local) else None
    random_at_first = None
    if first_k is not None:
        rr = random_df[(random_df.direction == 'down') & (random_df.k == first_k)]
        if len(rr):
            random_at_first = float(rr.empirical_tail_fraction.iloc[0])
    rows = [
        {
            'gate': 'G0_baseline_reproduction',
            'result': f'H0={h0_base:.4f}±{sig_base:.4f}',
            'status': 'PASS' if abs(h0_base-73.04) < 0.08 else 'FAIL',
            'meaning': 'Public SH0ES matrix solution reproduced.'
        },
        {
            'gate': 'G1_residual_required',
            'result': f'PEER+Local vs SH0ES residual={resid:.3f}±{resid_sig:.3f}, z={z:.2f}',
            'status': 'NO' if abs(z) < 2 else 'YES',
            'meaning': 'If NO, a nonzero calibration residual is not statistically required.'
        },
        {
            'gate': 'G2_single_SN_explains_residual',
            'result': f'max single-SN |ΔH0|={max_sn:.3f} vs residual={abs(resid):.3f}',
            'status': 'NO' if max_sn < abs(resid) else 'POSSIBLE',
            'meaning': 'Influence is not equivalent to bias.'
        },
        {
            'gate': 'G3_single_host_explains_residual',
            'result': f'max single-host |ΔH0|={max_host:.3f} vs residual={abs(resid):.3f}',
            'status': 'NO' if max_host < abs(resid) else 'POSSIBLE',
            'meaning': 'Tests concentration at host level.'
        },
        {
            'gate': 'G4_single_anchor_explains_residual',
            'result': f'max anchor-prior |ΔH0|={max_anchor:.3f} vs residual={abs(resid):.3f}',
            'status': 'NO' if max_anchor < abs(resid) else 'POSSIBLE',
            'meaning': 'Only geometric priors are removed; anchor Cepheids remain.'
        },
        {
            'gate': 'G5_small_SN_subset_reaches_local_band',
            'result': f'first greedy k={first_k}; random-subset tail fraction={random_at_first}' if first_k else 'no k<=5 reaches frozen Local-Hole band',
            'status': 'SENSITIVITY_SIGNAL' if first_k is not None else 'NO',
            'meaning': 'Greedy removal is adversarial and can only establish sensitivity concentration.'
        },
    ]
    if abs(z) < 2 and first_k is not None:
        classification = 'paper-grade sensitivity result'
    elif abs(z) >= 2 and first_k is not None and random_at_first is not None and random_at_first < 0.05:
        classification = 'mechanism candidate; independent replication required'
    else:
        classification = 'not promoted'
    rows.append({
        'gate': 'G6_promotion',
        'result': classification,
        'status': classification,
        'meaning': 'Mechanism evidence requires an independently predicted correction, not optimized data deletion.'
    })
    return pd.DataFrame(rows), z, classification


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--y', required=True)
    ap.add_argument('--L', required=True)
    ap.add_argument('--C', required=True)
    ap.add_argument('--q', required=True)
    ap.add_argument('--out', default='peer_superposition_paper')
    ap.add_argument('--random-draws', type=int, default=1000)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    y, L, C, params, sources = load_shoes(args.y, args.L, args.C, args.q)
    cache, h0_base, sig_base, five = baseline_fit(y, L, C, params)

    host_df, hosts = host_influence(y, L, params, sources, cache, h0_base, sig_base)
    groups = physical_sn_groups(sources, hosts)
    sn_df = single_sn_influence(y, L, params, cache, groups, h0_base, sig_base)
    anchor_df = anchor_influence(y, L, params, sources, cache, h0_base, sig_base)

    down = greedy_sequence(y, L, params, cache, groups, max_k=5, direction='down')
    up = greedy_sequence(y, L, params, cache, groups, max_k=5, direction='up')
    greedy_df = pd.concat([down, up], ignore_index=True)
    random_df = random_subset_nulls(y, L, params, cache, groups, greedy_df,
                                    n_draws=args.random_draws, seed=RNG_SEED)

    local_sigma = combined_local_sigma()
    gate_df, z_resid, classification = build_paper_gate(
        h0_base, sig_base, local_sigma, sn_df, host_df, anchor_df, greedy_df, random_df
    )

    resid_h0 = h0_base - LOCAL_MEDIAN
    resid_mu = equiv_dmu(h0_base, LOCAL_MEDIAN)
    # Linearized uncertainty in delta mu from the two H0 measurements.
    resid_mu_sigma = (5.0 / math.log(10.0)) * math.sqrt(
        (sig_base / h0_base)**2 + (local_sigma / LOCAL_MEDIAN)**2
    )

    sn_df.to_csv(out / 'calibrator_sn_influence.csv', index=False)
    host_df.to_csv(out / 'host_influence.csv', index=False)
    anchor_df.to_csv(out / 'anchor_influence.csv', index=False)
    greedy_df.to_csv(out / 'greedy_sn_removal.csv', index=False)
    random_df.to_csv(out / 'random_subset_nulls.csv', index=False)
    gate_df.to_csv(out / 'paper_gate.csv', index=False)

    summary = {
        'frozen_inputs': {
            'PEER_H0': PEER_H0,
            'PEER_sigma': PEER_SIGMA,
            'Local_H0_median': LOCAL_MEDIAN,
            'Local_H0_lo': LOCAL_LO,
            'Local_H0_hi': LOCAL_HI,
            'Target_reference_H0': TARGET_H0,
        },
        'shoes_baseline': {
            'H0': h0_base,
            'sigma_H0': sig_base,
            'five_log_H0': five,
            'n_obs': len(y),
            'n_params': len(params),
            'n_physical_calibrator_SNe': len(groups),
            'n_SN_hosts': len(hosts),
        },
        'superposition_residual': {
            'delta_H0_SH0ES_minus_PEERLocal': resid_h0,
            'sigma_combined_H0': math.sqrt(sig_base**2 + local_sigma**2),
            'z': z_resid,
            'equiv_delta_mu': resid_mu,
            'sigma_delta_mu_linearized': resid_mu_sigma,
            'zero_residual_rejected_2sigma': bool(abs(z_resid) >= 2),
        },
        'single_object_bounds': {
            'max_abs_SN_shift': float(sn_df.delta_H0.abs().max()),
            'max_abs_host_shift': float(host_df.delta_H0.abs().max()),
            'max_abs_anchor_shift': float(anchor_df.delta_H0.abs().max()),
        },
        'greedy_down': down.to_dict(orient='records'),
        'greedy_up': up.to_dict(orient='records'),
        'promotion': classification,
        'interpretation_rules': [
            'PEER and Local-Hole inputs are frozen before this battery.',
            'Greedy deletions are adversarial sensitivity bounds, not preferred cuts.',
            'Random-subset tail fractions are sensitivity percentiles, not discovery p-values.',
            'No correlated leave-one-out shifts are summed to claim a correction.',
            'A high influence score does not imply a bad datum.',
        ],
    }
    (out / 'summary.json').write_text(json.dumps(summary, indent=2, default=float))
    print(json.dumps(summary, indent=2, default=float))
    print('\nPAPER GATE\n' + gate_df.to_string(index=False))
    print('\nGREEDY DOWN\n' + down.to_string(index=False))
    print('\nRANDOM NULLS\n' + random_df[random_df.direction == 'down'].to_string(index=False))


if __name__ == '__main__':
    main()

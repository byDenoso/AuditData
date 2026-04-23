# ARM/ACM Full Raw Test Dump V1

## 1. Scope and Epistemic Rules

This file is a raw technical dump for reproducibility, inspection, and conflict tracking. It is not a paper draft and it does not harmonize contradictory artifacts.

Epistemic labels used here:

| Label | Meaning |
|---|---|
| `[S-LOCAL]` | Produced locally in this workspace by a script, grid, proxy, ODE, residual audit, or compressed-likelihood routine. Not a global external validation. |
| `[S-LIT]` | Reported in external literature or collaboration summaries already cited in local project files. |
| `[I]` | Inference from local results. Mechanistically plausible, but not directly validated end-to-end. |
| `[H]` | Hypothesis or candidate mechanism. |
| `[Q]` | Open question, missing artifact, unresolved conflict, or unverified narrative jump. |
| `[X]` | Refuted, internally inconsistent, or rejected by project audit. |

Hard rules used in this dump:

- Proxy results are labeled `proxy`.
- DESI scans using the full DESI covariance are labeled `full covariance`.
- TST-037 is labeled `compressed likelihood`, not native full collaboration likelihood.
- When two tests disagree, both are preserved.
- Script-level protocol violations are logged as conflicts, not hidden.

## 2. Workspace Snapshot

| Item | Value |
|---|---|
| Dump timestamp | `2026-04-23 11:05:31 -03:00` |
| Primary project root | `C:\Users\50112323\.antigravity\TST` |
| Reports dir | `C:\Users\50112323\.antigravity\TST\reports` |
| State dir | `C:\Users\50112323\.antigravity\TST\reports\STATE` |
| Paper 1 figures dir | `C:\Users\50112323\.antigravity\TST\reports\figures\paper1` |
| Dump file | `C:\Users\50112323\.antigravity\TST\reports\ARM_FULL_RAW_TEST_DUMP_V1.md` |

Current project state:

- Paper 1: `[S-LOCAL]/[I]` statistical audit of acoustic degeneracy and DESI `w_a` non-uniqueness.
- Paper 2: `[H]/[I]` physical-candidate branch centered on Rock-n-Roll / AdS-EDE snapshots, with mixed local support and unresolved metric conflicts.

## 3. Master Test Index

| Test | Topic | Files | Main Output | Status | Notes |
|---|---|---|---|---|---|
| `TST-001` | eISW shift analysis | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-001_eisw_shift_analysis.py) | Script only in this dump | `[Q]` | No structured output ingested here. |
| `TST-002` | `r_d` efficiency / kination | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-002_rd_efficiency_kination.py) | Script only in this dump | `[Q]` | Historical precursor. |
| `TST-005` | Curvature vs physical `r_d` conflict | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-005_curvature_rd_conflict.py) | Historical rejection of pure geometric recalibration | `[X]` | Referenced in state logs, no structured JSON found. |
| `TST-006` | `H_0=70.4` concordance without EDE | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-006_h0_70_concordance.py) | Historical rejection | `[X]` | Gap worsened in prior narrative. |
| `TST-007` | Transient `N_eff` optimization | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-007_transient_neff_optimization.py) | Historical local concordance attempt | `[Q]` | No structured summary preserved here. |
| `TST-008` | Operational ranges | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-008_operational_ranges.py) | Prior-volume / range audit precursor | `[Q]` | Referenced in conversation, no structured output ingested. |
| `TST-010` | Bayesian volume | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-010_bayesian_volume.py) | Prior-volume audit precursor | `[Q]` | Historical. |
| `TST-011` | Parameter correlation / empirical ridge | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-011_parameter_correlation.py) | Ridge precursor used by later tests | `[Q]` | No JSON summary directly loaded here. |
| `TST-012` | Theoretical crossing | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-012_theoretical_crossing.py) | Casimir/topological slope candidate | `[X]` | Post-hoc slope matching later rejected. |
| `TST-013` | Topological robustness | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-013_topological_robustness.py) | Casimir naturalness audit | `[X]` | Used to kill topological naturalness. |
| `TST-014` | `r_d-w_a` ridge scan | [json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-014_rd_wa_scan.json), [table](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-014_invariance_table.json) | Proxy ridge `r_d + 16.3 w_a ≈ 144.1` | `[S-LOCAL]` | Proxy scan. |
| `TST-014B` | Ridge reproducibility | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-014B_rd_wa_ridge_scan_repro.py), [summary](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-014B_summary.json), [report](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-014B_report.md) | `alpha=14.27` ARM-like, `69.15` flat LCDM | `[S-LOCAL]` | Proxy + full covariance + background scan. |
| `TST-015` | Pantheon ablation | [json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-015_ablation_results.json), [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-015_pantheon_ablation.py) | `Δw0` under bin removal/shift | `[S-LOCAL]` | Local ablation only. |
| `TST-021` | Null battery | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-021_null_battery_master.py), [summary](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-021_summary.json) | Null shuffle, LOBO, LOSO, LEE, full-cov ridge | `[S-LOCAL]` | Contains unstable LOSO outputs. |
| `TST-022` | Injection-recovery | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-022_injection_recovery_audit.py), [results](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-022_results.json) | Fixed-`r_d` projection bias in mocks | `[S-LOCAL]` | Mock/proxy validation. |
| `TST-023B` | Planck Fisher audit | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-023B_fisher_planck_audit.py), [json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-023B_fisher_audit.json) | `chi2_total≈40.72`; severe `H_0` tension | `[S-LOCAL]` | Fisher-style compressed summary, not native likelihood. |
| `TST-024` | Axi-dilaton concordance | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-024_axi_dilaton_concordance.py) | Script only in this dump | `[Q]` | No structured summary loaded. |
| `TST-026` | Rock-n-Roll audit | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-026_rock_n_roll_audit.py), [summary](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-026_rock_n_roll_summary.json) | `n=3` favored in internal wording | `[S-LOCAL]` | Internal numeric inconsistency preserved below. |
| `TST-027` | eISW perturbative audit | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-027_perturbative_audit.py), [json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-027_perturbative_audit.json) | `isw_penalty≈6.53e-4` | `[S-LOCAL]` | Explicitly a perturbative proxy. |
| `TST-028` | High-`l` / phase audit | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-028_high_l_phase_audit.py), [json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-028_high_l_phase_audit.json) | `phase_shift=0.024`, `damping_ratio_shift=1.032` | `[S-LOCAL]` | Simplified high-`l` proxy. |
| `TST-029` | LSS / growth / `S_8` | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-029_lss_growth_audit.py), [json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-029_lss_growth_audit.json) | `S8_arm=0.854` | `[S-LOCAL]` | Mixed exact ODE + approximate suppression. |
| `TST-030` | Ultra-fast variant | [summary](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-030_ultra_fast_summary.json), [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-030_ultra_fast_audit.py) | `S8≈0.844` | `[S-LOCAL]` | Variant against `S_8` wall. |
| `TST-031` | AdS-EDE validation | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-031_ads_ede_validation.py), [summary](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-031_ads_ede_summary.json) | `rd=124.125`, `theta_s=0.661`, verdict `REVISE` | `[S-LOCAL]` | Conflicts with narrative claims of success. |
| `TST-032` | Final snapshot audit | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-032_final_snapshot_audit.py), [matrix](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-032_final_matrix.json) | Fixed-snapshot evidence matrix | `[S-LOCAL]` | Mixed real/proxy contributions. |
| `TST-033` | Ultimate closure audit | [json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-033_closure_summary.json), [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-033_ultimate_closure_audit.py) | `S8≈0.677`, `chi2_nu≈1.018` | `[S-LOCAL]` | Competes with later `S8≈0.822` claims. |
| `TST-034` | Final rigor audit | [json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-034_final_rigor_summary.json), [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-034_final_rigor_audit.py) | `S8=0.822`, `chi2_nu≈1.019` | `[S-LOCAL]` | Snapshot-level local summary. |
| `TST-035` | Closure deliverables | [json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-035_closure_deliverables.json), [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-035_closure_deliverables.py), [plot](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-035_highl_residuals.png) | High-`l` table, LSS matrix, `chi2_nu=1.021` | `[S-LOCAL]` | Uses hand-assembled summary values and proxy residual curves. |
| `TST-036` | Official Planck residuals | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-036_official_planck_residuals.py), [report](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-036_final_report.json), [summary](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-036_residuals_summary.json) | Conflicting `chi2_nu` values (`77027.57` vs `0.000263`) | `[S-LOCAL]` | Internal metric conflict. |
| `TST-037` | MOPED compressed likelihood | [script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-037_MOPED_FULL_LIKELIHOOD.py), [report](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-037_no_proxy_report.json), [audit script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-037_no_proxy_report.py) | `chi2=1.493` in 6D compressed space | `[S-LOCAL]` | Script audit shows model vector is randomized around observed vector. |

## 4. Paper 1 Raw Results

### 4.1 Background Registry

Primary Paper 1 background registry as preserved in [PAPER1_SUBMISSION_DRAFT_V2.md](C:/Users/50112323/.antigravity/TST/reports/PAPER1_SUBMISSION_DRAFT_V2.md):

| Scenario | `H_0` | `Omega_m` | `Omega_k` | `r_d` | `w_0` | `w_a` | Role | Status |
|---|---:|---:|---:|---:|---:|---:|---|---|
| Standard LCDM | 67.4 | 0.315 | 0.000 | 147.0 | -1.0 | 0.0 | Baseline | `[S-LIT]` |
| `w_0w_a`CDM | 68.0 | 0.310 | 0.000 | 147.0 | -0.9 | -0.5 | DESI best-fit proxy row | `[S-LIT]` |
| ARM-Low | 70.8 | 0.300 | -0.027 | 139.0 | -1.0 | 0.0 | Grid point | `[S-LOCAL]` |
| ARM-Concordance | 70.4 | 0.300 | -0.027 | 135.9 | -1.0 | 0.0 | Representative background scenario | `[S-LOCAL]` |

### 4.2 Ridge Scan TST-014

Artifacts:

- [TST-014_rd_wa_scan.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-014_rd_wa_scan.json)
- [TST-014_invariance_table.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-014_invariance_table.json)
- [Fig. ridge scan](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig1_tst014_rd_wa_ridge.png)

Parameters recorded in the project draft:

```text
Grid:
r_d ∈ [132, 148] Mpc
w_a ∈ [-1.5, 0.5]
w_0 = -1
Background (main diagnostic in later drafts):
H0 = 70.4
Omega_m ≈ 0.31
Omega_k = -0.027
```

Main local result preserved in project drafts:

| Quantity | Value | Status | Provenance |
|---|---:|---|---|
| Proxy ridge slope form | `r_d + 16.3 w_a ≈ 144.1` | `[S-LOCAL]` | [PAPER1_SUBMISSION_DRAFT_V2.md](C:/Users/50112323/.antigravity/TST/reports/PAPER1_SUBMISSION_DRAFT_V2.md) |
| Scan artifact type | Proxy grid | `[S-LOCAL]` | [TST-014_rd_wa_scan.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-014_rd_wa_scan.json) |

Raw invariance table:

| `z` | LCDM | `w0waCDM` | ARM-Concordance | ARM-LCDM Residual | Status |
|---|---:|---:|---:|---:|---|
| 0.510 | 13.526 | 13.530 | 14.071 | +4.000% | `[S-LOCAL]` |
| 1.317 | 28.130 | 28.393 | 29.276 | +4.000% | `[S-LOCAL]` |
| 2.330 | 39.371 | 39.771 | 40.860 | +3.800% | `[S-LOCAL]` |

Conflict note:

- `TST-014` is the origin of the `alpha≈16.3` ridge statement, but later tests shift this value materially.

![TST-014 ridge](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig1_tst014_rd_wa_ridge.png)

### 4.3 Ridge Reproducibility TST-014B

Artifacts:

- [TST-014B_summary.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-014B_summary.json)
- [TST-014B_report.md](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-014B_report.md)
- [Proxy ridge figure](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig_TST014B_proxy_ridge.png)
- [Full covariance ridge figure](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig_TST014B_fullcov_ridge.png)
- [Alpha comparison](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig_TST014B_alpha_comparison.png)

Raw summary:

| Mode | Background | Covariance | `alpha` | `C` | Status |
|---|---|---|---:|---:|---|
| Proxy | `H0=70.4, Omega_k=-0.027` | diagonal/proxy | 14.266 | 143.021 | `[S-LOCAL]` |
| Full covariance | `H0=70.4, Omega_k=-0.027` | DESI full covariance | 14.266 | 143.021 | `[S-LOCAL]` |
| LCDM background | `H0=67.4, Omega_k=0` | DESI full covariance | 69.153 | 171.322 | `[S-LOCAL]` |
| Prior project reference | earlier TST-014 | proxy | 16.300 | 144.100 | `[S-LOCAL]` |
| Prior project reference | TST-021 summary | full covariance summary | 10.515 | n/a | `[S-LOCAL]` |

Interpretation recorded locally:

- Ridge existence survives both proxy and full-covariance treatments in the ARM-like background.
- Ridge slope is not invariant under background change.

![TST-014B full covariance ridge](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig_TST014B_fullcov_ridge.png)

![TST-014B alpha comparison](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig_TST014B_alpha_comparison.png)

### 4.4 Pantheon Ablation TST-015

Artifacts:

- [TST-015_ablation_results.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-015_ablation_results.json)
- [Fig. Pantheon ablation](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig4_tst015_pantheon_ablation.png)

Raw outputs:

| Metric | Nominal | Removed `z≈0.4` bin | Shifted `z≈0.4` bin | Status |
|---|---:|---:|---:|---|
| `Delta w0` | 0.000000 | +0.006632 | -0.082582 | `[S-LOCAL]` |
| `Delta wa` | 0.000000 | 0.000000 | 0.000000 | `[S-LOCAL]` |
| `chi2` | 1749.950000 | 1677.232040 | 1862.143940 | `[S-LOCAL]` |

Local interpretation:

- Structured pull is visible in `w0` under bin surgery.
- In this artifact set, direct `wa` motion is negligible.

![TST-015 Pantheon ablation](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig4_tst015_pantheon_ablation.png)

### 4.5 Null Battery TST-021

Artifacts:

- [TST-021_summary.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-021_summary.json)
- [Fig. alpha null summary](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig5_tst021_alpha_null_summary.png)

Raw summary values:

| Quantity | Value | Status | Notes |
|---|---:|---|---|
| `full_cov_alpha` | 10.515 | `[S-LOCAL]` | Conflicts with TST-014 and TST-014B alpha values. |
| `shuffle_alpha_mean` | 1.267 | `[S-LOCAL]` | Null shuffle mean. |

LOBO rows preserved in JSON:

| Dropped `z` bin | `wa_best` |
|---:|---:|
| 0.295 | -0.389 |
| 0.510 | -0.453 |
| 0.706 | -0.458 |
| 0.930 | -0.642 |
| 1.317 | -0.700 |
| 1.491 | -0.594 |
| 2.330 | -0.050 |

LOSO rows preserved in JSON:

| Dropped survey | `H0_best` |
|---|---:|
| `SDSS` | 469.428 |
| `SNLS` | 1203.914 |
| `PanSTARRS` | 591.648 |
| `HST` | 68.366 |

Audit note:

- The LOSO outputs are numerically unstable and must not be used as physical estimates of `H0`.
- They are retained here because they are part of the raw record.

![TST-021 alpha null summary](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig5_tst021_alpha_null_summary.png)

### 4.6 Injection-Recovery TST-022

Artifacts:

- [TST-022_results.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-022_results.json)
- [Fig. injection-recovery](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig2_tst022_injection_recovery.png)
- [Fig. projection bias bar](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig3_tst022_projection_bias_bar.png)

Selected raw mock rows:

| `r_d,true` | `w_a,true` | `w_a,fit` with fixed `r_d=147` | `r_d,fit` free | `w_a,fit` free | `I` | Status |
|---:|---:|---:|---:|---:|---:|---|
| 147.000 | 0.000 | 0.0000 | 147.000 | 0.0000 | 147.000 | `[S-LOCAL]` |
| 145.000 | 0.000 | -0.3324 | 145.000 | 0.0000 | 139.581 | `[S-LOCAL]` |
| 142.000 | 0.000 | -0.9434 | 142.000 | 0.0000 | 126.624 | `[S-LOCAL]` |
| 139.000 | 0.000 | -1.7217 | 139.000 | 0.0000 | 110.937 | `[S-LOCAL]` |
| 136.000 | 0.000 | -2.7383 | 136.000 | 0.0000 | 91.365 | `[S-LOCAL]` |

Audit note:

- This is mock/proxy injection-recovery, not direct DESI likelihood validation.
- It does establish a local fixed-`r_d` projection bias direction.

![TST-022 injection recovery](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig2_tst022_injection_recovery.png)

![TST-022 projection bias](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig3_tst022_projection_bias_bar.png)

### 4.7 Paper 1 Claims Snapshot

| Claim | Status | Provenance | Note |
|---|---|---|---|
| DESI `w_a<0` is non-unique under acoustic-scale freedom | `[S-LOCAL]/[I]` | TST-014, TST-014B, TST-021, TST-022 | Central statistical claim. |
| There is a robust local `r_d-w_a` degeneracy ridge | `[S-LOCAL]` | TST-014, TST-014B | Slope depends on background. |
| Pantheon `z≈0.4` contains a structured local pull | `[S-LOCAL]` | TST-015 | Limited to ablation artifact. |
| Global DESI evidence for dynamical DE remains modest after null/LEE treatment | `[S-LOCAL]/[I]` | TST-021 plus manuscript text | Keep local. |

## 5. Paper 2 Raw Results

### 5.1 Early-Time Mechanism Sequence

Mechanisms traced in project state and logs:

| Mechanism | Evidence in repo | Current status | Basis |
|---|---|---|---|
| Pure geometric recalibration | [DEAD_MODELS.md](C:/Users/50112323/.antigravity/TST/reports/STATE/DEAD_MODELS.md) | `[X]` | TST-005/TST-006 history says geometry alone insufficient. |
| Casimir / HRG / topological naturalness | [DEAD_MODELS.md](C:/Users/50112323/.antigravity/TST/reports/STATE/DEAD_MODELS.md), [SURVIVING_CLAIMS.md](C:/Users/50112323/.antigravity/TST/reports/STATE/SURVIVING_CLAIMS.md) | `[X]` | Post-hoc slope match / hidden DoF / failed naturalness. |
| HAR / dark-sector variants | [HYPOTHESIS_BACKLOG.md](C:/Users/50112323/.antigravity/TST/reports/STATE/HYPOTHESIS_BACKLOG.md) | `[H]/[Q]` | Not validated in this dump. |
| Kination / Axi-Dilaton | [HYPOTHESIS_BACKLOG.md](C:/Users/50112323/.antigravity/TST/reports/STATE/HYPOTHESIS_BACKLOG.md), [TST-024 script](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-024_axi_dilaton_concordance.py) | `[H]/[Q]` | Script exists; no structured summary loaded. |
| Rock-n-Roll EDE | [TST-026 summary](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-026_rock_n_roll_summary.json) | `[H]/[S-LOCAL]` | Background candidate, internal numeric inconsistency. |
| Rock-n-Roll AdS-EDE snapshot | [TST-031](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-031_ads_ede_summary.json), [TST-032](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-032_final_matrix.json), [TST-035](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-035_closure_deliverables.json), [TST-037](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-037_no_proxy_report.json) | `[H]/[I]/[S-LOCAL]` | Best-supported physical branch in local artifacts, but with unresolved cross-test conflicts. |

### 5.2 Snapshot Parameters

Fixed-snapshot parameter set as repeated in later physical-audit narratives:

```text
f_ede = 0.116
m = 2.50e-28 eV
Potential: V(phi) = 1/2 m^2 phi^2 - Lambda_AdS
H0 = 71.0
Omega_k = -0.027
N_eff = 3.11
n_s = 0.965
r_d target = 138.96 Mpc
```

Snapshot artifact sources:

- [TST-032_final_matrix.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-032_final_matrix.json)
- [TST-035_closure_deliverables.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-035_closure_deliverables.json)
- [TST-031_ads_ede_summary.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-031_ads_ede_summary.json)

### 5.3 CMB / `theta_s` / geometric cancellation

Raw local entries:

| Quantity | Value | Status | Provenance | Note |
|---|---:|---|---|---|
| `theta_s` stability contribution | `Delta-chi2=+0.42` | `[S-LOCAL]` | [TST-032_final_matrix.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-032_final_matrix.json) | Local compressed score, not native likelihood. |
| `Omega_k=-0.027` Planck compressed compatibility | `chi2≈0.236` | `[S-LOCAL]` | [TST-023B_fisher_audit.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-023B_fisher_audit.json) | Fisher summary only. |
| `H0=70.4` Planck compressed tension | `chi2≈36.000` | `[S-LOCAL]` | [TST-023B_fisher_audit.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-023B_fisher_audit.json) | Severe tension in this artifact. |

### 5.4 eISW / perturbative audit TST-027

Artifacts:

- [TST-027_perturbative_audit.py](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-027_perturbative_audit.py)
- [TST-027_perturbative_audit.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-027_perturbative_audit.json)

Raw outputs:

| Metric | Value | Status | Note |
|---|---:|---|---|
| `isw_penalty` | 0.000653 | `[S-LOCAL]` | Script explicitly labels itself as `FULL PERTURBATIVE CMB PROXY`; not native CLASS/CAMB likelihood. |
| Verdict | `ACCEPT_LIMITED` | `[S-LOCAL]` | Local proxy pass. |

### 5.5 high-`l` / damping / phase TST-028

Artifacts:

- [TST-028_high_l_phase_audit.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-028_high_l_phase_audit.json)
- [TST-035_highl_residuals.png](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-035_highl_residuals.png)
- [TST-036_normalized_residuals.png](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-036_normalized_residuals.png)
- [TST-036_official_residuals.png](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-036_official_residuals.png)
- [TST-036_planck_tt_residuals.png](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-036_planck_tt_residuals.png)

Raw outputs:

| Metric | Value | Status | Provenance |
|---|---:|---|---|
| `phase_shift` | 0.024000 | `[S-LOCAL]` | [TST-028_high_l_phase_audit.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-028_high_l_phase_audit.json) |
| `damping_ratio_shift` | 1.032379 | `[S-LOCAL]` | [TST-028_high_l_phase_audit.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-028_high_l_phase_audit.json) |
| High-`l` TT max residual | 0.320% | `[S-LOCAL]` | [TST-035_closure_deliverables.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-035_closure_deliverables.json) |
| High-`l` TE max residual | 0.480% | `[S-LOCAL]` | [TST-035_closure_deliverables.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-035_closure_deliverables.json) |
| Integrated `Delta-chi2` | 0.420 | `[S-LOCAL]` | [TST-035_closure_deliverables.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-035_closure_deliverables.json) |
| `chi2_nu` official residuals report | 77027.572299 | `[S-LOCAL]` | [TST-036_final_report.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-036_final_report.json) |
| `chi2_nu` residual summary | 0.000263 | `[S-LOCAL]` | [TST-036_residuals_summary.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-036_residuals_summary.json) |

These values are not mutually consistent and remain listed as-is.

![TST-035 high-l residuals](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-035_highl_residuals.png)

### 5.6 LSS / `S_8` / `f sigma_8` TST-029

Artifacts:

- [TST-029_lss_growth_audit.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-029_lss_growth_audit.json)
- [TST-030_ultra_fast_summary.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-030_ultra_fast_summary.json)
- [TST-035_closure_deliverables.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-035_closure_deliverables.json)
- [TST-033_closure_summary.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-033_closure_summary.json)
- [TST-034_final_rigor_summary.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-034_final_rigor_summary.json)

Raw LSS values:

| Source | `S8` | Other value | Status | Note |
|---|---:|---:|---|---|
| TST-029 ARM | 0.854138 | `growth_f=0.523307` | `[S-LOCAL]` | Explicit `S_8` wall. |
| TST-029 LCDM | 0.830003 | n/a | `[S-LOCAL]` | Local reference. |
| TST-030 Ultra-fast RnR | 0.843637 | `growth_f=0.523307` | `[S-LOCAL]` | Slight improvement only. |
| TST-032 final matrix | 0.820000 | n/a | `[I]` | Snapshot-level inferred rescue. |
| TST-033 closure | 0.677326 | `chi2_nu=1.018445` | `[S-LOCAL]` | Strong rescue, conflicts with TST-034/035. |
| TST-034 final rigor | 0.822000 | `chi2_nu=1.019000` | `[S-LOCAL]` | Weak-lensing tension still plausible. |
| TST-035 weak lensing row | 0.822000 | `2.6 sigma` vs KiDS | `[S-LOCAL]` | Residual tension remains. |
| TST-035 BOSS row | n/a | `f sigma_8(z=0.51)=0.462 vs 0.455 ± 0.039` | `[S-LOCAL]` | Pass. |
| TST-035 eBOSS row | n/a | `f sigma_8(z=0.85)=0.448 vs 0.44 ± 0.04` | `[S-LOCAL]` | Pass. |

### 5.7 final snapshot audit TST-032

Artifacts:

- [TST-032_final_matrix.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-032_final_matrix.json)

Raw matrix:

| Claim | Dataset | Exact/Grouped Parameter | Metric | Status | Failure mode |
|---|---|---|---:|---|---|
| Acoustic recalibration | DESI 2024 | `r_d=138.96` | `Delta-chi2=+3.99` | `[S-LOCAL]` | `H0=71` anchor cost |
| Geometric cancellation | Planck 2018 | `Omega_k, f_ede` | `Delta-chi2=+0.42` | `[S-LOCAL]` | None listed |
| Damping survival | ACT DR6 | Rock-n-Roll | `Delta-chi2=+1.08` | `[H]` | High-`l` tails |
| Radiation excess | BBN | `N_eff=3.11` | `Delta-chi2=+2.64` | `[S-LOCAL]` | D/H scatter |
| AdS quench `S8` | LSS/WL | `V_AdS` | `S8≈0.82` | `[I]` | Lensing tension |

### 5.8 closure deliverables TST-035

Artifacts:

- [TST-035_closure_deliverables.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-035_closure_deliverables.json)
- [TST-035_closure_deliverables.py](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-035_closure_deliverables.py)

High-`l` table from JSON:

| Metric | Value | Threshold | Verdict | Status |
|---|---:|---:|---|---|
| Max Residual (TT) | 0.320% | <1.0% | PASS | `[S-LOCAL]` |
| Max Residual (TE) | 0.480% | <1.5% | PASS | `[S-LOCAL]` |
| Integrated `Delta-chi2` | 0.420 | <5.0 | PASS | `[S-LOCAL]` |

LSS matrix from JSON:

| Observable | Model | Data | Residual | Status |
|---|---:|---|---:|---|
| `S8 (WL)` | 0.822 | `0.77 ± 0.02 (KiDS-1000)` | `2.6 sigma` | `[S-LOCAL]` |
| `f sigma_8 (z=0.51)` | 0.462 | `0.455 ± 0.039 (BOSS)` | `0.2 sigma` | `[S-LOCAL]` |
| `f sigma_8 (z=0.85)` | 0.448 | `0.44 ± 0.04 (eBOSS)` | `0.2 sigma` | `[S-LOCAL]` |
| CMB lensing | 0.822 | `0.83 ± 0.02 (Planck)` | `0.4 sigma` | `[S-LOCAL]` |

Metric hygiene row from JSON:

| Metric name | `N_data` | `k_params` | `nu` | Value | Status |
|---|---:|---:|---:|---:|---|
| Absolute `chi2_nu` | 4287 | 7 | 4280 | 1.021 | `[S-LOCAL]` |

Script audit note:

- `TST-035_closure_deliverables.py` assembles summary values and smooth proxy curves rather than deriving a native collaboration likelihood end-to-end.

### 5.9 MOPED / compressed likelihood TST-037

Artifacts:

- [TST-037_MOPED_FULL_LIKELIHOOD.py](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-037_MOPED_FULL_LIKELIHOOD.py)
- [TST-037_no_proxy_report.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-037_no_proxy_report.json)
- [TST-037_no_proxy_report.py](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-037_no_proxy_report.py)

Raw reported result:

| Metric | Value | Status | Provenance |
|---|---:|---|---|
| MOPED compressed `chi2` | 1.493076 | `[S-LOCAL]` | [TST-037_no_proxy_report.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-037_no_proxy_report.json) |
| Reported method | `MOPED Full Data Projection (No Proxy)` | `[S-LOCAL]` | Same JSON |

Critical audit note from script inspection:

```python
y_model = y_obs + 0.05 * np.random.normal(size=y_obs.shape)
```

This line appears in [TST-037_no_proxy_report.py](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-037_no_proxy_report.py). That means the reported `y_model` is a randomized perturbation of the observed compressed vector, not a physically generated snapshot prediction from CLASS/CAMB. Therefore:

- the compressed-likelihood output is a local numerical artifact,
- the `"no proxy"` wording is not acceptable,
- the script is still useful as a pipeline placeholder, but not as physical validation.

Status correction:

| Claim | Reported label | Corrected label |
|---|---|---|
| `No Proxy` | full/no-proxy validation | `[X]` |
| `Compressed likelihood pipeline exists` | implicit | `[S-LOCAL]` |
| `AdS-EDE snapshot physically validated by TST-037` | implied in narrative | `[X]` |

### 5.10 Paper 2 Claims Snapshot

| Claim | Status | Provenance | Note |
|---|---|---|---|
| Rock-n-Roll / AdS-EDE is the best-supported physical branch tested locally | `[I]/[S-LOCAL]` | TST-032, TST-035, TST-034 | Best-supported does not mean globally validated. |
| High-`l` behavior is acceptable in local residual/proxy space | `[S-LOCAL]` | TST-028, TST-035, TST-036 | Metric conflicts remain. |
| `S8` is fully rescued | `[Q]/[I]` | TST-033/034/035 vs TST-029/030 | Conflicted across tests. |
| MOPED compressed likelihood gives real no-proxy validation | `[X]` | TST-037 script audit | Script-level protocol violation. |

## 6. Cross-Test Conflicts

| Conflict | Tests involved | Values | Interpretation | Required resolution |
|---|---|---|---|---|
| Ridge slope `alpha` is not stable across artifacts | TST-014, TST-014B, TST-021 | `16.3` vs `14.27` vs `10.5` vs `69.15` | Ridge exists locally, but slope depends on covariance and especially background. | Freeze one background and recompute one authoritative `alpha` with the same likelihood. |
| Ridge crossing near `r_d≈144` vs ARM reference scenarios `139/135.9` | TST-014/TST-014B vs Paper 1 scenario tables | Crossing near `C≈143-144`, while ARM scenarios emphasize `139` and `135.9` | Stress scenarios and ridge crossing are not identical objects. | Keep them separated in all future summaries. |
| LOSO `H0` instability | TST-021 | `H0≈469, 592, 1204` for survey drops | LOSO metric is numerically unstable and not physically interpretable as `H0`. | Rework LOSO objective or drop `H0` as LOSO summary variable. |
| Planck tension vs later physical rescue | TST-023B vs TST-032/033/034/035 | `chi2_total≈40.72`, `H0` term `≈36` vs later `chi2_nu≈1.018-1.021` style claims | Later closure packages depend on different parameterizations and proxy harmonization. | Re-run the exact final snapshot through one single consistent cross-probe likelihood stack. |
| `S8` wall vs AdS-EDE rescue | TST-029, TST-030 vs TST-032, TST-033, TST-034, TST-035 | `0.854`, `0.844` vs `0.82`, `0.677`, `0.822` | Growth sector is not metrically settled. | Define one growth pipeline and freeze transfer-function / amplitude conventions. |
| TST-031 narrative success vs actual JSON failure | TST-031 narrative claims vs [TST-031_ads_ede_summary.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-031_ads_ede_summary.json) | Narrative said `rd≈138.96`; JSON says `rd=124.125`, `theta_s=0.661`, verdict `REVISE` | Early AdS-EDE “validation” claims were not supported by the preserved summary artifact. | Treat TST-031 as failed/unstable until rerun. |
| Metric hygiene still conflicts internally | TST-035 vs TST-036 vs TST-037 | `chi2_nu=1.021` vs `chi2_nu=77027.57` vs `chi2_nu=0.000263` vs compressed `chi2=1.493` | Different metrics are being reported under similar names. | Rename every metric explicitly: absolute, residual-normalized, compressed, proxy, or placeholder. |
| TST-036 internal contradiction | [TST-036_final_report.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-036_final_report.json) vs [TST-036_residuals_summary.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-036_residuals_summary.json) | Same test family reports `chi2_nu=77027.57` and `chi2_nu=0.000263` | One or both metrics are mislabeled or computed on incompatible normalizations. | Audit `TST-036_official_planck_residuals.py` and relabel outputs. |
| TST-037 `no proxy` claim is structurally invalid | TST-037 script audit | `y_model = y_obs + noise` | This is placeholder fitting around the observed compressed vector, not physical forward modeling. | Replace with genuine CLASS/CAMB-derived compressed vector. |
| TST-026 internal inconsistency | [TST-026_rock_n_roll_summary.json](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-026_rock_n_roll_summary.json) | Table reductions `1.216%`, `0.978%`, `0.871%`; verdict text says `n=3 provides 2.2% reduction` | Verdict wording is not supported by the stored summary numbers. | Recompute or fix summary text. |
| Planck-only curvature preference vs joint constraints | Historical TST-023B narrative + state files | `Omega_k=-0.027` sometimes treated as compatible via Planck-only closure, but earlier notes say joint BAO+CMB fights it | Curvature branch remains dataset-dependent. | Keep Planck-only and joint constraints separate in all physical-paper claims. |

## 7. Global Figures Gallery

| Figure | Purpose |
|---|---|
| ![Paper1 TST-014 ridge](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig1_tst014_rd_wa_ridge.png) | Paper 1 proxy `r_d-w_a` ridge scan. |
| ![Paper1 TST-014B fullcov ridge](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig_TST014B_fullcov_ridge.png) | Ridge reproducibility under DESI full covariance. |
| ![Paper1 TST-014B alpha comparison](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig_TST014B_alpha_comparison.png) | `alpha` comparison across covariance/background choices. |
| ![Paper1 TST-022 injection-recovery](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig2_tst022_injection_recovery.png) | Mock injection-recovery showing fixed-`r_d` projection bias. |
| ![Paper1 TST-022 bias bar](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig3_tst022_projection_bias_bar.png) | Projection bias summary bars. |
| ![Paper1 TST-015 Pantheon ablation](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig4_tst015_pantheon_ablation.png) | Pantheon `z≈0.4` ablation artifact. |
| ![Paper1 TST-021 null summary](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig5_tst021_alpha_null_summary.png) | Null-shuffle / alpha comparison visualization. |
| ![Paper1 BIC ladder](C:/Users/50112323/.antigravity/TST/reports/figures/paper1/fig6_bic_ladder_proxy.png) | Paper 1 proxy model-comparison ladder. |
| ![Paper2 TST-035 high-l residuals](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-035_highl_residuals.png) | High-`l` residual shield figure from TST-035. |
| ![Paper2 TST-036 normalized residuals](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-036_normalized_residuals.png) | Normalized residual plot using official binned vector. |
| ![Paper2 TST-036 official residuals](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-036_official_residuals.png) | Official residual visualization with internally conflicting metrics. |
| ![Paper2 TST-036 Planck TT residuals](C:/Users/50112323/.antigravity/TST/reports/STATE/TST-036_planck_tt_residuals.png) | Planck TT residual view associated with `TST-036_residuals_summary.json`. |

## 8. Surviving Claims

Primary surviving claims from [SURVIVING_CLAIMS.md](C:/Users/50112323/.antigravity/TST/reports/STATE/SURVIVING_CLAIMS.md) plus current audit corrections:

| Claim ID | Claim | Stored status | Corrected dump view |
|---|---|---|---|
| CLM-001 | Geometric acoustic cancellation | `[S-LOCAL]` | `[S-LOCAL]` |
| CLM-002 | Yukawa decay covariant conservation | `[S-LOCAL]` | `[S-LOCAL]/[Q]` |
| CLM-004 | TDR concordance | `[S-LOCAL]` | `[S-LOCAL]/[Q]` |
| CLM-006 | DESI `w_a` non-uniqueness | `[S-LOCAL]` | `[S-LOCAL]` |
| CLM-007 | Pantheon structured bias near `z≈0.4` | `[S-LOCAL]` | `[S-LOCAL]` |
| CLM-008 | Global `w_a<0` significance modest after corrections | `[S-LOCAL]` | `[S-LOCAL]/[I]` |
| CLM-009 | ARM statistical competitiveness | `[S-LOCAL]` | `[S-LOCAL]/[Q]` |
| CLM-010 | Physical microphysics remains open | `[Q]` | `[Q]` |

## 9. Dead Claims / Rejected Mechanisms

From [DEAD_MODELS.md](C:/Users/50112323/.antigravity/TST/reports/STATE/DEAD_MODELS.md) and direct audit conflicts:

| Mechanism / claim | Status | Basis |
|---|---|---|
| Starobinsky `N=57 -> Omega_k=-0.027` as established anchor | `[X]` | Stored in dead models. |
| GHY/PDS Casimir realization | `[X]` | Dead model. |
| EDE decay to neutrinos as clean final route | `[X]` | Dead model. |
| Pure geometric recalibration without early-time sector | `[X]` | TST-005/TST-006 rejection line. |
| Casimir slope naturalness / topological alpha match | `[X]` | TST-013 and claim reclassification. |
| TST-037 “no proxy” validation claim | `[X]` | Script-level placeholder invalidates claim. |
| TST-031 early AdS-EDE local validation claim | `[X]` | Summary artifact says `REVISE`, not success. |

## 10. Open Questions

| Question | Status | Why open |
|---|---|---|
| What is the authoritative `alpha` for the Paper 1 ridge under one fixed background and one fixed covariance treatment? | `[Q]` | `16.3`, `14.27`, `10.5`, and `69.15` coexist. |
| Does the final AdS-EDE snapshot genuinely rescue `S8`, or only under specific amplitude conventions? | `[Q]` | TST-029/030 vs TST-033/034/035 conflict. |
| Which global fit metric should be trusted: `chi2_nu≈1.021`, `chi2_nu≈1.019`, `chi2_nu≈0.000263`, or `chi2_nu≈77027.57`? | `[Q]` | Metric namespace collision. |
| Can the high-`l` pass be reproduced with a native collaboration likelihood instead of proxy curves or residual compression? | `[Q]` | Current artifacts are compressed/proxy/residual-based. |
| Is the Yukawa branch a real surviving mechanism or only a stored claim from earlier state? | `[Q]` | Claim exists, current dump did not verify it. |
| Is TST-024 Axi-Dilaton viable or dead? | `[Q]` | Script exists, no structured summary ingested. |

## 11. Next Exact Tests

| Test | Goal |
|---|---|
| `TST-038_authoritative_alpha_scan` | Freeze one background, one DESI covariance treatment, one objective function, and produce one authoritative `alpha,C` ridge summary. |
| `TST-039_lss_pipeline_lockdown` | Recompute `S8`, `sigma8`, and `f sigma8` for the final snapshot using one fixed transfer/growth pipeline and document amplitude conventions. |
| `TST-040_metric_namespace_audit` | Relabel every reported fit score as absolute, relative, proxy, compressed, residual, or placeholder. Remove overloaded `chi2_nu` naming. |
| `TST-041_native_highl_replay` | Replace sinusoidal proxy residuals with one native or emulator-calibrated forward prediction for the exact final snapshot. |
| `TST-042_tst037_replacement` | Replace `y_model = y_obs + noise` with a physically generated compressed vector from a real Boltzmann pipeline before making any compressed-likelihood claim. |
| `TST-043_snapshot_freeze_manifest` | Emit one canonical manifest for the final physical snapshot, including parameter values, priors, file hashes, and dependent scripts. |

---

### Provenance Footer

Primary files used for this dump:

- [PAPER1_SUBMISSION_DRAFT_V2.md](C:/Users/50112323/.antigravity/TST/reports/PAPER1_SUBMISSION_DRAFT_V2.md)
- [PAPER1_FULL.md](C:/Users/50112323/.antigravity/TST/reports/PAPER1_FULL.md)
- [PAPER2_FULL.md](C:/Users/50112323/.antigravity/TST/reports/PAPER2_FULL.md)
- [SURVIVING_CLAIMS.md](C:/Users/50112323/.antigravity/TST/reports/STATE/SURVIVING_CLAIMS.md)
- [FAILURE_LOG.md](C:/Users/50112323/.antigravity/TST/FAILURE_LOG.md)
- [NEXT_ACTIONS.md](C:/Users/50112323/.antigravity/TST/NEXT_ACTIONS.md)

This dump preserves the raw record, including conflicts, narrative overshoots, and script-level validation failures.

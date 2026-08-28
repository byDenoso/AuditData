# PEER H0 Superposition Paper Battery Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Test whether the residual between a frozen PEER global H0 and the local distance-ladder H0 is better described as a superposition of observer-environment and calibration/sample effects than as one monolithic H0 offset.

**Architecture:** Keep the physical layers separate. Freeze PEER H0=70.391±0.801 km/s/Mpc, use the previously derived Local-Hole observer envelope 71.712035–72.170378 with median 71.928041, reproduce the public SH0ES 2022 GLS solution from its y/L/C/q matrices, then run exact correlated deletions of calibrator SNe, hosts, and geometric-anchor priors. Add adversarial subset-removal nulls and random-subset controls so influence concentration is measured rather than inferred from summed jackknife shifts.

**Tech Stack:** Python 3.12, NumPy, SciPy, pandas, pytest, GitHub Actions.

**Spec:** This conversation's PEER + Local Hole + SH0ES calibration-superposition hypothesis, frozen before the final battery.

## Global Constraints

- PEER background stays fixed at H0=70.391 and sigma=0.801.
- Local-Hole observer result stays fixed at median 71.92804067892604, envelope [71.71203503729637, 72.17037753571651].
- SH0ES data source is the public SH0ES2022 y_R22/L_R22/C_R22/q_R22 matrix representation.
- Never sum correlated leave-one-out shifts as if independent.
- A large leave-one-out shift means influence, not bad data.
- H0 must not be used to select Local-Hole geometry.
- Adversarial greedy deletions are sensitivity bounds, not preferred data cuts.

---

### Task 1: Add paper-battery contracts

**Files:**
- Create: `tests/test_peer_superposition_battery.py`
- Modify: none

**Interfaces:**
- Consumes: `tools.shoes_host_jackknife.fit_delete_cached`, `get_h0`, `physical_calibrator_key`.
- Produces: tests for `equiv_dmu`, `combined_local_sigma`, and deterministic subset ranking.

- [ ] Write failing tests for magnitude conversion, uncertainty propagation, and stable ranking by signed H0 influence.
- [ ] Run `PYTHONPATH=. pytest -q tests/test_peer_superposition_battery.py` and verify RED.
- [ ] Commit the tests.

### Task 2: Implement exact superposition battery

**Files:**
- Create: `tools/peer_superposition_battery.py`
- Test: `tests/test_peer_superposition_battery.py`

**Interfaces:**
- Consumes: SH0ES y/L/C/q matrices and frozen PEER/local constants.
- Produces: `summary.json`, `calibrator_sn_influence.csv`, `host_influence.csv`, `anchor_influence.csv`, `greedy_upward_sn_removal.csv`, `random_subset_nulls.csv`, `paper_gate.csv`.

- [ ] Implement helpers required by Task 1.
- [ ] Reproduce SH0ES baseline H0≈73.0434±1.0072.
- [ ] Compute PEER+Local vs SH0ES residual in H0 and equivalent magnitude, with uncertainty propagation.
- [ ] Compute exact single host/SN/anchor influence using correlated GLS deletion.
- [ ] Build an adversarial greedy sequence that repeatedly removes the calibrator SN whose removal most lowers H0, recomputing the GLS after every deletion.
- [ ] For each greedy subset size k=1..5, draw 1000 random physical-calibrator-SN subsets of the same size and compute an empirical lower-tail p value for the greedy H0 shift.
- [ ] Run mirror greedy removals that increase H0 as a symmetry/control check.
- [ ] Write a paper gate that distinguishes: residual required/not required, single-object explanation, concentrated-tail explanation, anchor explanation, and sensitivity-only result.
- [ ] Run the contract tests and verify GREEN.
- [ ] Commit implementation.

### Task 3: Add reproducible workflow

**Files:**
- Create: `.github/workflows/nexo-peer-superposition-paper.yml`

**Interfaces:**
- Consumes: public SH0ES2022 matrix files.
- Produces: uploaded artifact `nexo-peer-superposition-paper`.

- [ ] Download y_R22.txt, L_R22.txt, C_R22.txt, q_R22.txt from the public source.
- [ ] Run both test files relevant to the battery.
- [ ] Run `tools/peer_superposition_battery.py`.
- [ ] Print `summary.json` and `paper_gate.csv` for auditability.
- [ ] Upload all outputs as an artifact.
- [ ] Commit workflow.

### Task 4: Scientific verification and promotion gate

**Files:**
- No new source files.

**Interfaces:**
- Consumes: workflow artifact and logs.
- Produces: final scientific verdict.

- [ ] Verify workflow success and exact baseline reproduction.
- [ ] Verify greedy subsets are compared against random subsets of identical size.
- [ ] Verify no correlated single-object shifts are summed to claim a solution.
- [ ] Check whether PEER+Local residual remains <2σ before any adversarial cut.
- [ ] Check whether any one SN, host, or anchor can explain the full residual.
- [ ] Check whether a small SN subset can move SH0ES into PEER+Local range and whether that subset is extreme under the random-subset null.
- [ ] Classify the result as `paper-grade mechanism evidence`, `paper-grade sensitivity result`, or `not promoted`.

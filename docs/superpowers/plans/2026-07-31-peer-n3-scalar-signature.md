# PEER n=3 Scalar Perturbation Signature Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Identify a reproducible CMB signature caused by the linear perturbations of the canonical PEER scalar with n=3.

**Architecture:** Generate CAMB 1.6.6 spectra for the canonical scalar, an exact-background no-perturbation build, a background-matched AxionEffectiveFluid control, and n=2.8/3.2 controls. Project the exact scalar perturbation residual against standard cosmological derivative directions, then rank the surviving TT/TE/EE/lensing fingerprint by multipole band.

**Tech Stack:** Python 3.12, CAMB 1.6.6, NumPy, SciPy, Matplotlib, pytest, GitHub Actions Ubuntu 24.04.

## Global Constraints

- Canonical slice: n=3, fde_zc=0.0880, log10(zc)=3.81, theta_i=2.89155, H0=70.795 km/s/Mpc.
- Exact no-perturbation control must reuse the scalar background and differ only by num_perturb_equations=0.
- Background equality gate: H(z), Omega_PEER(z), H0, rdrag and thetastar must agree to numerical tolerance.
- Nuisance projection basis: lnAs, ns, ombh2, omch2, tau, H0 and Alens.
- Results are a theoretical CMB fingerprint, not an observed detection or likelihood preference.

---

### Task 1: Analysis primitives and tests

**Files:**
- Create: `peer_n3_signature/core.py`
- Create: `tests/test_peer_n3_signature.py`

**Interfaces:**
- Produces `gaussian_covariance_ttee`, `weighted_project`, `dominant_period`, and `assert_background_equal`.

- [ ] **Step 1: Write failing tests** for exact span projection, covariance positivity, oscillation-period recovery and background mismatch rejection.
- [ ] **Step 2: Run** `pytest -q tests/test_peer_n3_signature.py` and verify RED.
- [ ] **Step 3: Implement the minimal primitives.**
- [ ] **Step 4: Re-run the focused suite and verify GREEN.**

### Task 2: CAMB spectrum generator

**Files:**
- Create: `peer_n3_signature/generate.py`

**Interfaces:**
- Produces one NPZ spectrum file and one JSON metadata file per model variant.

- [ ] **Step 1: Generate** stock scalar n=2.8, 3.0 and 3.2 spectra, LCDM, finite-difference nuisance spectra and the optimized AxionEffectiveFluid control.
- [ ] **Step 2: Export** lensed/unlensed TT, TE, EE, BB, lensing potential, backgrounds and derived parameters.
- [ ] **Step 3: Verify** finite outputs through a smoke run at lmax=300.

### Task 3: Exact-background perturbation ablation

**Files:**
- Create: `scripts/build_camb_nopert.sh`

**Interfaces:**
- Produces an isolated CAMB 1.6.6 environment where only TQuintessence perturbation equations are disabled.

- [ ] **Step 1: Download** the official CAMB 1.6.6 source and record its git commit.
- [ ] **Step 2: Patch only** `TQuintessence_Init` so `num_perturb_equations=0` while preserving background evolution.
- [ ] **Step 3: Build** the patched source in an isolated virtual environment and generate n=3 spectra.
- [ ] **Step 4: Require** the exact-background gate before analysis.

### Task 4: Signature extraction and report

**Files:**
- Create: `peer_n3_signature/analyze.py`
- Create: `peer_n3_signature/__init__.py`

**Interfaces:**
- Consumes generated NPZ/JSON files.
- Produces `signature_metrics.json`, `signature_spectra.csv`, figures, `REPORT.md`, manifests and SHA-256 hashes.

- [ ] **Step 1: Form** scalar-full minus scalar-no-perturbations and scalar-full minus matched-fluid residuals.
- [ ] **Step 2: Project** the exact residual against the nuisance derivative basis using the Gaussian TT/TE/EE covariance.
- [ ] **Step 3: Measure** channel/band signal norms, dominant oscillation period, zero crossings, peak locations and correlation with dC/dn.
- [ ] **Step 4: Write** machine-readable results, plots and scientific boundary text.

### Task 5: Reproducible execution lane

**Files:**
- Create: `.github/workflows/peer-n3-scalar-signature.yml`

**Interfaces:**
- Produces a single immutable evidence artifact.

- [ ] **Step 1: Pin** runtime dependencies.
- [ ] **Step 2: Run** tests, stock generation, no-perturbation build and final analysis.
- [ ] **Step 3: Validate** the evidence contract and hashes.
- [ ] **Step 4: Upload** the complete artifact with 90-day retention.

# Native-Windows CAMB reproducibility artifact

Public hardening artifact for **ASCOM-D-26-00323 — Reproducible CAMB-Based Cosmological Inference on Native Windows**.

This repository exists to let an external reviewer audit the environment identity, frozen content hashes, numerical tolerances, validation gates, and paper-facing replay without relying on workstation-specific paths.

## Scope

This artifact covers only the **native Windows reference environment**:

- CPython 3.12.13;
- CAMB source commit `3ef0272d6f7ba1231128872e56e6d4c12af8267b`;
- Cobaya 3.6.2;
- released CMB and late-time likelihood interfaces used by the registered gates;
- numerical parity, exact reuse, workload-aware benchmarks, reconstruction, relocation, and paper replay.

It **does not** contain NOVA's cross-environment certification, backend authority, promotion/rollback governance, JAX/Rust development, or nomadic reconstruction claims.

## Evidence policy

Historical results from the submitted manuscript are stored as `historical_pass_needs_r1_replay`. They are references, not fresh R1 evidence. A paper-facing result is emitted as `verified` only after a new receipt is produced by the hardening run.

Current headline reference points include:

- Windows M0: `log L = -395.4831594065799`, `chi2 = 790.966`;
- ACT kernel historical speedup: `5.5659x`, objective drift `1.42e-13`;
- `r_drag` continuation: `14 -> 1` CAMB background calls;
- likelihood identity grid: 64 points, registered tolerance `1e-9`;
- bounded minimizer battery: 8 contexts, registered tolerance `1e-7`;
- profile identity: 7 `Omega_m` points, registered tolerance `1e-7`;
- balanced Cobaya historical benchmark: `1.8111x`, shape calls `14 -> 2`.

## Fast start

```powershell
./bootstrap.ps1
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\python.exe scripts\verify_artifact.py
.\.venv\Scripts\python.exe scripts\doctor.py --strict
.\.venv\Scripts\python.exe scripts\run_validation.py --strict
.\.venv\Scripts\python.exe scripts\reproduce_tables.py
.\.venv\Scripts\python.exe scripts\reproduce_figures.py
```

`run_validation.py --strict` fails closed while a required R1 gate remains unverified. That is intentional.

See [REPRODUCING.md](REPRODUCING.md), [docs/SCOPE.md](docs/SCOPE.md), and [validation/validation_matrix.csv](validation/validation_matrix.csv).

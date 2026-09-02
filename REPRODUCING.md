# Reproducing ASCOM-D-26-00323

## 1. Clean Windows host

Use a fresh Windows user or VM. Do not reuse the original `C:\Users\Dener` or `D:\CODEX` paths.

Required starting tools:

- Git
- CPython 3.12.13
- a Fortran/C toolchain only if the pinned CAMB source must be rebuilt locally
- Microsoft MPI runtime when exercising MPI gates

## 2. Clone and bootstrap

```powershell
git clone <public-repository-url>
cd <repository>
./bootstrap.ps1
```

The bootstrap creates `.venv`, installs the public Python dependency set, then installs CAMB from the frozen upstream commit. It does not infer private workstation paths.

## 3. Materialize external scientific data

Large or collaboration-managed payloads are not silently vendored. `manifests/external-components.json` records the component, identity, redistribution status, and expected hash when known.

Before likelihood execution, place each allowed payload under `external-data/` and run:

```powershell
.\.venv\Scripts\python.exe scripts\verify_artifact.py --external-data external-data
```

Hash mismatch or a missing required payload is a hard failure.

## 4. Environment doctor

```powershell
.\.venv\Scripts\python.exe scripts\doctor.py --strict --output results\doctor.json
```

The doctor records OS, architecture, Python identity, package versions, CAMB source identity when discoverable, and external payload status. A non-Windows host cannot produce a native-Windows acceptance receipt.

## 5. R1 validation battery

```powershell
.\.venv\Scripts\python.exe scripts\run_validation.py --strict
```

The matrix covers clean bootstrap, environment identity, theory/likelihood parity, optimization identity, workload benchmarks, relocation, cache behavior, repeated runs, failure modes, and paper replay.

The submitted-paper numbers are comparison targets only. Fresh receipts go in `results/r1/` and must be marked `verified` only after the registered tolerance is met.

## 6. Paper replay

```powershell
.\.venv\Scripts\python.exe scripts\reproduce_tables.py
.\.venv\Scripts\python.exe scripts\reproduce_figures.py
```

These scripts refuse to promote historical/pending values into the paper-facing surface.

## 7. Release

After every required gate is `verified`, freeze the exact commit and create tag:

```text
v1.0-ascom-major-revision
```

The journal manuscript should cite the immutable tag/commit, not `main`.

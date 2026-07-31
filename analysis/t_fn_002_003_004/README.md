# PEER T-FN-002 / T-FN-003 / T-FN-004

Canonical status: `LOCAL_N3_MODE_ROBUST_NO_PLANCK_ONLY_DETECTION`.

## Largest supported claim

The Planck 2018 TTTEEE plik-lite likelihood contains a local physical basin in `(f_PEER,n)` whose posterior mode remains near `n=3` across eleven independently anchored exact CAMB 1.6.6 tangent bases.

## Results

- Parent theory surface: 81/81 exact CAMB spectra.
- New tangent sweep: 8 bases and 112 exact derivative spectra.
- Total tangent bases tested: 11.
- Initial boundary calibration: 6,000,000 mocks.
- Final tangent sweep calibration: 22,000,000 mocks.
- Total mock draws executed: 28,000,000.
- Fresh tests after package extraction: 7/7 pass.

Posterior geometry across all eleven bases:

- MAP `f_PEER = 0.065–0.086`;
- MAP `n = 3.00–3.06`;
- common local 68% interval: `f_PEER = 0.0533–0.0851`;
- common local 68% interval: `n = 2.943–3.262`.

Boundary-null calibration:

- empirical p envelope: `0.000040–0.095928`;
- conservative supremum p: `0.0959279`;
- one-sided Gaussian equivalent: `Z = 1.305`;
- basis-robust rejection at 5%: `false`.

The dependence on tangent anchor is smooth rather than an isolated bad basis. The predeclared cubic geometry gate gives `R² = 0.998948` and maximum standardized residual `1.730`.

## Claim boundary

`f_PEER` and `n` are evaluated through continuous interpolation of 81 exact CAMB spectra. Standard cosmological directions are analytically marginalized using exact local CAMB tangent bases. This is not a full chain that calls CAMB at every posterior state. It supports local identifiability of the `n≈3` mode, but does not support a Planck-only detection of `f_PEER>0`.

## Decision

Preserve the `n≈3` local-mode result in the microphysics paper. Do not quote Planck-only sigma or detection language. The next material gate for the main model paper is matched multi-probe/high-l inference. A fully nonlinear Planck-only CAMB chain is only required if a formal Planck-only posterior width is needed.

## Provenance

- Parent package SHA-256: `0c421b04320e31067219a60f22914c4896de187745160e89cc2a52db1768d6c7`.
- Tangent sweep PR: `#17`.
- Head commit: `3b1431b755c696eca2e58ccbd2bb62c18c03db69`.
- Workflow run: `30669335390`.
- Workflow artifact: `8808219920`.
- Artifact digest: `sha256:de240801d82b390a42d8fe5e6d276ca413c7aeea47bd12677d1f8d27dcf998d0`.
- Reproducibility package SHA-256: `22a53924bc622d159f34f9a8fe7973328b06f8d6636f240a988ba8e57f341e42`.

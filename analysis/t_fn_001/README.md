# T-FN-001 exact Planck f_PEER × n profile

Canonical status: `LOCAL_N3_BASIN_ROBUST_NULL_DELTA_BASIS_SENSITIVE`.

Largest supported claim: the exact physical Planck likelihood surface contains a robust local basin centered on n=3. Profiling f_PEER does not shift the preferred native index away from n=3 under canonical, null-anchored, or best-point cosmological tangent bases.

Execution evidence:
- exact theory grid: 81/81 CAMB 1.6.6 spectra;
- Planck likelihood: native 2018 TTTEEE plik-lite v22, 613 bandpowers;
- robustness: three complete 81-cell surfaces, 243/243 valid fits;
- optimizer: ten deterministic starts per cell, 2,430 attempts;
- best f_PEER across bases: 0.060–0.088;
- best n in every basis: 3.0;
- common local 68% intersection: f_PEER=0.080–0.100 and n=2.90–3.10;
- canonical point (0.088,3.0): Delta chi-square=0.000–0.246;
- observed improvement over f_PEER=0: Delta chi-square=3.400–8.474.

The null magnitude is tangent-basis-sensitive. No p-value, sigma, or discovery claim is assigned because n is unidentified at f_PEER=0, the null lies on a boundary, and the standard-cosmology profile is locally linearized.

Provenance:
- grid PR15, workflow 30664736113, artifact 8806572789;
- tangent PR16, workflow 30665307777, artifact 8806742786;
- local reproducibility package SHA-256: `0504024491f9ce39c469cf5cecd869002079455b4e3548149e0b3077535ccd6e`;
- generator lineage commit: `16d1df948ba17762ff430b9de293cba799ed1e7e`.

Next gate: full nonlinear n-free sampling with exact CAMB updates across the f-n ridge, followed by boundary-null calibration with mocks.

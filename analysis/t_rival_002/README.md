# T-RIVAL-002: model-specific nuisance tangent comparison

Canonical decision: `NO_COHERENT_MODEL_SPECIFIC_PREFERENCE`.

This test removes the main profiling asymmetry in T-RIVAL-001. The canonical PEER scalar n=3 and the acoustically near-matched w=1/2 effective-fluid rival are each profiled using their own exact CAMB 1.6.6 finite-difference nuisance tangent basis.

Observable state:
- 11 scalar tangent bases;
- one fluid-native tangent basis with 14 derivative spectra;
- Planck 2018 TTTEEE plik-lite, ACT DR6 CMB-only TT/TE/EE, and SPT-3G D1 lite TT/TE/EE;
- 66 valid joint model fits, all full rank and without bound pinning;
- 99 conditional channel-attribution blocks;
- 9/9 tests after clean package extraction; 453/453 package hashes valid.

Result:
- Planck Delta chi-square scalar minus fluid: -4.294 to -3.066, scalar favored in 11/11 bases;
- ACT DR6: +1.750 to +5.032, fluid favored in 11/11 bases;
- SPT-3G: -0.535 to +2.779, fluid favored in 8/11 and scalar in 3/11;
- independent-profile diagnostic sum: -2.996 to +4.745, crossing zero.

The cross-experiment conflict persists after each model receives its own nuisance tangent basis. Therefore the T-RIVAL-001 negative conclusion was not caused by profiling the fluid with scalar derivatives.

Claim boundary: the effective-fluid rival is acoustically near-matched, not identical in the full background history. Cross-experiment sums remain diagnostic because shared nuisances and cross-covariances are not included.

Next material gate: port one exact native AdS-EDE or NEDE rival with its own background, perturbations, and tangent basis, then repeat the matched Planck+ACT+SPT comparison.

Reproducibility package SHA-256: `6476f9ea74865543a1daf043c823faca46f5e0825b79001fbd18137f3343fe57`.

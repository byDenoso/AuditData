# T-GROWTH-002 + T-LENS-001

Canonical decision: `PLANCK_SPT_PEER_GROWTH_AND_LENSING_CONSISTENT`.

The Planck+SPT preferred PEER cell, f_PEER=0.060 with n=3, log10(zc)=3.81 and theta_i=2.89155, was evaluated with exact CAMB 1.6.6 linear and HMCode/mead2020 nonlinear calculations under three tau priors.

Observable state:
- 30 exact CAMB runs;
- own-best-fit LCDM, PEER and native axion-EDE;
- matched f_PEER=0 and matched axion-shape controls at fixed PEER late cosmology;
- ten redshifts over z=0-2;
- physical lensing predictions always use A_lens=1;
- 4/4 tests in the workspace and 4/4 after clean package extraction;
- 92/92 package files hash-verified and ZIP integrity passed.

Growth result:
- maximum matched |Delta D/D| = 1.269e-6;
- maximum matched |Delta f/f| = 1.688e-6;
- PEER S8=0.8309-0.8344 across tau priors;
- diagonal BOSS DR12 three-point chi-square=1.562-1.717.

Lensing result:
- PEER pull against the published Planck lensing-only amplitude combination=1.086-1.218 sigma;
- native axion pull=1.110-1.241 sigma;
- matched native axion versus PEER nonlinear C_L^phiphi shape residual over L=8-400 <3.42e-4;
- matched f_PEER=0 versus PEER retains an approximately 4.25% shape residual after amplitude removal, showing a real early-transfer/lensing morphology rather than a modified late growth law.

The initial bounded execution stopped after 15/30 spectra because of the wall-clock limit. All checkpoints were finite and validated. The run resumed from cache and completed the remaining 15 spectra; no partial result was promoted.

Claim boundary: the BOSS check is diagonal and the Planck lensing check uses the published one-dimensional lensing-only combination plus theoretical C_L^phiphi morphology. These are validated closure diagnostics, not substitutes for full correlated full-shape or Planck lensing likelihood sampling.

Reproducibility package SHA-256: `770c66ad46be7d23f08ed003412515acdefc8baf548f5b32f420f75925c72c48`.

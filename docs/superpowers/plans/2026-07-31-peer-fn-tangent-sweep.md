# T-FN-004 exact tangent-basis sweep

## Objective

Determine whether the boundary-null p-value variation found by T-FN-003 is an isolated artifact of the canonical, null, and best-point tangent anchors or a smooth property of the PEER f-n likelihood geometry.

## Frozen inputs

- Parent exact grid: T-FN-001, 81 CAMB 1.6.6 spectra.
- Planck likelihood: TTTEEE plik-lite v22.
- Fixed PEER coordinates: log10(zc)=3.81, theta_i=2.89155.
- Standard directions: lnAs, ns, ombh2, omch2, tau, H0, Alens.
- Finite-difference steps identical to T-FN-001.

## New exact anchors

- (f,n)=(0.02,3.0), (0.04,3.0), (0.08,2.9), (0.08,3.0), (0.08,3.4), (0.10,3.0), (0.12,3.0), (0.14,3.0).

Existing exact bases at f=0, f=0.06,n=3, and f=0.088,n=3 remain independent controls.

## Completion contract

- 8/8 anchor jobs complete.
- 14 derivative spectra and 14 metadata files per anchor.
- 112/112 finite CAMB spectra.
- CAMB version exactly 1.6.6.
- Complete hashes and immutable workflow artifact.
- Scientific interpretation only after local replay of T-FN-002 and T-FN-003 for every basis.

## Decision rule

- If p varies smoothly with anchor position, classify tangent sensitivity as a genuine local-geometry effect and use the supremum p as the conservative result.
- If one or two anchors are isolated outliers, inspect derivative conditioning and regenerate those bases before changing the scientific claim.
- No workflow success alone is promoted into a scientific result.

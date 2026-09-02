# Scope boundary: Native Windows vs NOVA

## Owned here

- native Windows host integration;
- CPython/CAMB/Cobaya environment reconstruction;
- likelihood/data identities used by the paper;
- numerical equivalence at registered points, minima and profiles;
- exact reuse and workload-aware performance;
- clean-room bootstrap, relocation, cache and failure-mode tests;
- paper table/figure replay.

## Explicitly outside this repository

- NOVA TheoryRequest/TheoryProduct architecture as a general cross-environment abstraction;
- nomadic hydration as an independent framework contribution;
- backend promotion/rollback governance;
- backend authority semantics;
- JAX or Rust solver-development claims;
- cross-environment certification beyond the optional matched Windows-Linux extension.

Shared CAMB/Cobaya vocabulary and historical lineage do not transfer NOVA claims into ASCOM-D-26-00323.

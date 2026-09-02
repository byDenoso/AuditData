# Provenance

The public artifact was reconstructed from the paper-facing records rather than from a raw workstation directory.

Canonical historical inputs used to define the R1 contract:

- `Reproducible_CAMB_Inference_on_Native_Windows.pdf` (submitted-paper snapshot);
- `Supplementary_Information_Native_Windows_CAMB.docx`;
- `PEER_CLOSURE_EXECUTION_20260729_v2.md` for the native Windows M0 gate;
- `pip_freeze_cobaya_camb_py312.txt` for the historical package inventory;
- `LH-Data` PR75/PR83 for the runtime/benchmark lineage.

Important hardening decision: PR75/PR83 contain Linux-oriented execution history and explicitly left native CPython 3.12/Windows acceptance open for some benchmarks. Those historical performance values are therefore comparison targets, not fresh R1 Windows evidence.

The frozen CAMB source commit is public upstream commit `3ef0272d6f7ba1231128872e56e6d4c12af8267b`.

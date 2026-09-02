# ASCOM-00323 R1 hardening implementation plan

1. Freeze paper/NOVA scope boundary.
2. Normalize the historical Windows environment and remove absolute paths from public manifests.
3. Freeze upstream CAMB source commit and content hashes.
4. Build fail-closed validation/comparison helpers with tests first.
5. Encode every paper claim as a validation gate and R1 receipt path.
6. Add Windows bootstrap, doctor, artifact verification and paper replay scripts.
7. Publish on an isolated public branch without modifying the existing public repo main branch.
8. Run static/unit verification locally.
9. Run native-Windows R1 gates; keep publication readiness false until every required gate is verified.
10. Freeze tag `v1.0-ascom-major-revision` only after the full R1 matrix closes.

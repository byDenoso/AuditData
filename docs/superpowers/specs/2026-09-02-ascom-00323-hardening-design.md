# ASCOM-00323 R1 hardening design

Goal: publish a reviewer-auditable Native Windows reproduction artifact without importing NOVA's distinct scientific contribution.

Architecture: a small portable harness plus frozen manifests and an explicit external-data boundary. Heavy/restricted scientific payloads are verified by hash, not hidden behind local paths. Historical results enter as non-promoted references. Fresh Windows receipts are the only source allowed to populate the revised paper.

Happy path: clone -> bootstrap -> external-data verification -> doctor -> validation matrix -> paper replay -> immutable release.

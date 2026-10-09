# doc230-cbio043-embryo-arrest
partial: transcriptome only; batch-dominated. Audited (RESULTS.md): D9 resveratrol contrast permutation p=0.005, D3 source separability confirmed; no verdict changed.
partial: transcriptome only; arrest confounded with data source (D3 control shows batch dominates). See PREREG.md, RESULTS.md, src/run.py. Not an extension of mega27-14-digital-embryo (different study).

## Reproduce / readability audit (2026-10-09)
Data: sources and SHA256 are in data_manifest/. Scripts read fixed paths under /tmp (download the files there first, check the SHA256, then run python3 src/run.py and python3 src/audit.py). Not independently re-run from a clean machine as part of this readability check; the numbers in RESULTS.md come from the original runs on 2 CPU / 2 GB. Figshare Data S2.tsv (zip) is unzipped into /tmp/emb/.

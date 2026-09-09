---
name: dna-only mode
description: Added --type dna-only flag to support DNA-only oncoanalyser cohorts; found and fixed a silent data-loss bug in the RNA integration join
type: project
---

Added `params.type` (`both` default, or `dna-only`) so the pipeline can run on nf-core/oncoanalyser output that has no RNA (Isofox, SAGE RNA-append) — motivated by the CARDINAL/DRAC discussion where oncoanalyser can be run DNA-only or DNA+RNA, WGS or WES.

**Why:** `GENOMIC_MUTATIONS` (`subworkflows/local/genomic_mutations/main.nf`) joined the SAGE RNA-append VCF against the DNA MAF with a plain `.join()` (inner join). For any subject lacking RNA, that join produced nothing — the subject's mutation output silently vanished from the pipeline (not degraded, just absent). This was found by tracing the actual channel logic, not by trusting the user's initial guess that "Neo, CUPPA, CHORD, Sigs" were the RNA-dependent steps — those names don't exist anywhere in this repo (they're nf-core/oncoanalyser's own upstream tools); Sigs here is actually DNA-only (PAVE somatic VCF), confirmed via `docs/mutational_signatures.md`.

**How to apply:** Any future RNA-related change should assume `params.type == 'dna-only'` is a real, tested code path — see [[oncoanalyser migration]]. `dna-only` is a global flag (not per-subject); it skips RNA channel construction in `workflows/genomic.nf`, skips `INTEGRATE_RNA_VARIANTS` entirely in `genomic_mutations`, and omits the `expression` meta-file entry in `genomic_aggregate_output` so no orphan `meta_expression.txt` is packaged. Tests added in `tests/subworkflows/genomic_mutations.nf.test` and `tests/subworkflows/genomic_aggregate_output.nf.test` — both deliberately supply RNA/TPM input under `type=dna-only` to prove it's an override, not just a fallback for missing files.

Not yet done: WES/targeted vs WGS support (separate ask, out of scope for this change) — see nf-core/oncoanalyser `--mode` (`wgs`/`wts`/`wgts`/`targeted`) if that's revisited.

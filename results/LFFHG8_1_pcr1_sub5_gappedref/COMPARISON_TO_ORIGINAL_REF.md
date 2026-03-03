# Comparison: original vs gapped reference

## Inputs
- Reads: `data/reads/LFFHG8_fastq/LFFHG8_1_pcr1_sub5.fastq`
- Original reference: `data/references/gblock_f1r1.gb`
- Gapped reference: `data/references/gblock_f1r1_gapped.gb`

## High-level result
- Mapped read count is unchanged: `4 / 4586` (0.09%) for both references.
- The mapped read IDs are unchanged: `795_LFFHG8_1`, `3446_LFFHG8_1`, `3685_LFFHG8_1`, `3741_LFFHG8_1`.

## Minority-read quality comparison
Per-read comparison is in `mapped4_ref_compare.tsv`.

Observed changes:
- `3446_LFFHG8_1`: essentially unchanged.
- `795_LFFHG8_1`: essentially unchanged.
- `3685_LFFHG8_1`: identity improves substantially (0.8563 -> 0.9338).
- `3741_LFFHG8_1`: alignment extends further (query coverage 0.2454 -> 0.3957), with lower identity due to longer aligned span.

Aggregate over the 4 mapped reads:
- Mean identity: `0.9532` (original) -> `0.9569` (gapped)
- Mean query coverage: `0.7204` (original) -> `0.7580` (gapped)
- Length-weighted identity: `0.9589` (original) -> `0.9653` (gapped)

## Interpretation
The gapped/ambiguous-barcode reference does not increase the number of target-like reads, but it provides a modest improvement in alignment quality/coverage for the existing minority of mapped reads.

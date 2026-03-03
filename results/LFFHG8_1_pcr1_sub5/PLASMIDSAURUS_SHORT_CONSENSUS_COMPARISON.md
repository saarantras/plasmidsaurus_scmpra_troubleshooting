# Plasmidsaurus short consensus vs in-repo analyses

Date: 2026-03-03

## Input short consensus
`GAGTTCAGACGTGTGCTCTTCCGATCTGGCAACAGTATATTGAGGAATGAGGGGTCAGAAGGACACTCCAACCACCTCCTTCATCATTTTTTTTCTCCCTAGTCATTTTTAGGGATTGACTCACTGTTAGAGAACGTAGGGAGGTGGTTGGACTGACCCCTCATAGATCGGAAGAGCACACGTCTGAACTCA`

This sequence exactly matches:
- `data/references/LFFHG8_plasmidsaurus_draft.fasta`

## Comparison to the 1611 bp contaminant consensus
Compared against:
- `results/LFFHG8_1_pcr1_sub5/LFFHG8_1_pcr1_sub5.consensus.fasta`

Findings:
- No full-length 192 bp exact match.
- Best shared exact block is at the adapter-like end (31 bp at the left flank).
- Human 18S-like hallmark motifs used to identify the 1611 bp contaminant are absent from this 192 bp sequence.

## Comparison to gBlock reference
Compared against:
- `data/references/gblock_f1r1_gapped.fasta`

Findings:
- No full-length 192 bp match.
- Only short adapter/primer segments match (e.g., 27 bp `AGATCGGAAGAGCACACGTCTGAACTC`).

## Read-level support in this FASTQ
Dataset:
- `data/reads/LFFHG8_fastq/LFFHG8_1_pcr1_sub5.fastq`

Summary:
- 0 reads contain the full 192 bp sequence exactly.
- Many reads contain the two flank regions with variable internal sequence/length.
- Using BLAST HSP aggregation at >=85% identity, strongest reads typically cover ~121-128 bp of the 192 bp query (two blocks: left and right ends).

Interpretation:
- This is likely a related but distinct adapter-flanked contaminant family in the same library, rather than an alternative consensus call of the same 1611 bp human-rDNA-dominant sequence.

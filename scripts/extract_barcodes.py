#!/usr/bin/env python3
"""
Extract barcode sequences from aligned reads at N-masked reference positions.

For each read in the BAM, walks the CIGAR to find query bases at the target
reference interval, then counts and reports barcode frequencies.

Usage:
    python scripts/extract_barcodes.py \
        --bam results/MCSS56_1_sub8/MCSS56_1_sub8.bam \
        --ref-start 32 --ref-end 48 \
        --out-prefix results/MCSS56_1_sub8/barcodes \
        [--min-mapq 10] [--require-full]
"""
import argparse
import pysam
import pandas as pd
from collections import Counter
from pathlib import Path


def extract_barcode(read, ref_start, ref_end):
    """
    Return the query bases spanning ref positions [ref_start, ref_end).
    Returns None if the read doesn't cover the full interval or has only deletions.
    Deletions within the interval are represented as '-'.
    """
    if read.is_unmapped:
        return None

    # get_aligned_pairs returns (query_pos, ref_pos); query_pos is None for deletions
    aligned = read.get_aligned_pairs(matches_only=False, with_seq=False)

    # index by ref_pos
    ref_to_query = {}
    for qpos, rpos in aligned:
        if rpos is not None and ref_start <= rpos < ref_end:
            ref_to_query[rpos] = qpos  # None means deletion

    if len(ref_to_query) < (ref_end - ref_start):
        return None  # read doesn't span the full interval

    seq = read.query_sequence
    bases = []
    for rpos in range(ref_start, ref_end):
        qpos = ref_to_query[rpos]
        if qpos is None:
            bases.append('-')
        else:
            bases.append(seq[qpos])
    return ''.join(bases)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bam', required=True)
    ap.add_argument('--ref-start', type=int, required=True, help='0-based start of barcode region')
    ap.add_argument('--ref-end',   type=int, required=True, help='0-based end (exclusive) of barcode region')
    ap.add_argument('--out-prefix', required=True)
    ap.add_argument('--min-mapq', type=int, default=10, help='Minimum MAPQ (default: 10)')
    ap.add_argument('--require-full', action='store_true',
                    help='Discard reads where any barcode position is a deletion')
    args = ap.parse_args()

    bam = pysam.AlignmentFile(args.bam, 'rb')
    counts = Counter()
    per_read = []
    n_total = n_pass_mapq = n_extracted = 0

    for read in bam.fetch():
        if read.is_secondary or read.is_supplementary:
            continue
        n_total += 1
        if read.mapping_quality < args.min_mapq:
            continue
        n_pass_mapq += 1

        bc = extract_barcode(read, args.ref_start, args.ref_end)
        if bc is None:
            continue
        if args.require_full and '-' in bc:
            continue

        n_extracted += 1
        counts[bc] += 1
        per_read.append({'read_id': read.query_name, 'barcode': bc,
                         'mapq': read.mapping_quality})

    bam.close()

    bc_len = args.ref_end - args.ref_start
    rows = [{'barcode': bc, 'count': cnt,
             'fraction': cnt / n_extracted if n_extracted else 0,
             'has_deletion': '-' in bc,
             'has_N': 'N' in bc}
            for bc, cnt in counts.most_common()]
    df = pd.DataFrame(rows)

    out = Path(args.out_prefix)
    df.to_csv(f'{out}.counts.tsv', sep='\t', index=False)
    pd.DataFrame(per_read).to_csv(f'{out}.per_read.tsv', sep='\t', index=False)

    # summary
    n_unique = len(counts)
    n_clean = int(df.loc[~df['has_deletion'] & ~df['has_N'], 'count'].sum())
    top5 = df.head(5)[['barcode', 'count', 'fraction']].to_string(index=False)

    summary_lines = [
        f'total_primary\t{n_total}',
        f'pass_mapq\t{n_pass_mapq}',
        f'barcode_extracted\t{n_extracted}',
        f'unique_barcodes\t{n_unique}',
        f'clean_barcodes_reads\t{n_clean}',
        f'barcode_region\t{args.ref_start}-{args.ref_end} (0-based, len={bc_len})',
        f'min_mapq\t{args.min_mapq}',
    ]
    Path(f'{out}.summary.txt').write_text('\n'.join(summary_lines) + '\n')

    print('\n'.join(summary_lines))
    print(f'\nTop 5 barcodes:\n{top5}')
    print(f'\nWrote: {out}.counts.tsv, {out}.per_read.tsv, {out}.summary.txt')


if __name__ == '__main__':
    main()

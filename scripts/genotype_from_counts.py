#!/usr/bin/env python3
"""
genotype_from_counts.py

Genotype TE junctions from BBDuk k-mer counting results.
Parses BBDuk stats files to extract matched read counts, then calls
genotypes based on the ratio of absence vs presence allele reads.

Usage:
    python genotype_from_counts.py <results_dir> [-o genotype_results.tsv]
"""

import os
import sys
import glob
import argparse
from collections import defaultdict


def parse_bbduk_stats(stats_file):
    """Parse BBDuk stats file to get #Matched read count."""
    with open(stats_file) as f:
        for line in f:
            if line.startswith('#Matched'):
                parts = line.strip().split('\t')
                return int(parts[1])
    return 0


def call_genotype(abs_count, pre_count, min_depth=5, het_min=0.2, het_max=0.8):
    """
    Call genotype based on allele counts.

    Returns (genotype_string, notes_string).
    """
    total = abs_count + pre_count

    if total < min_depth:
        return "NO_CALL", f"depth={total}"

    abs_ratio = abs_count / total

    if abs_ratio >= het_max:
        return "Abs/Abs", f"abs_ratio={abs_ratio:.3f}"
    elif abs_ratio <= het_min:
        return "Pre/Pre", f"abs_ratio={abs_ratio:.3f}"
    else:
        return "Abs/Pre", f"abs_ratio={abs_ratio:.3f}"


def main():
    parser = argparse.ArgumentParser(
        description='Genotype TE junctions from BBDuk k-mer counts')
    parser.add_argument('results_dir',
        help='Directory containing *_abs_stats.txt and *_pre_stats.txt')
    parser.add_argument('-o', '--output', default=None,
        help='Output file (default: <results_dir>/genotype_results.tsv)')
    parser.add_argument('--min-depth', type=int, default=5,
        help='Minimum total read depth for calling (default: 5)')
    parser.add_argument('--het-min', type=float, default=0.2,
        help='Min abs ratio for heterozygote call (default: 0.2)')
    parser.add_argument('--het-max', type=float, default=0.8,
        help='Max abs ratio for heterozygote call (default: 0.8)')

    args = parser.parse_args()

    abs_files = sorted(glob.glob(
        os.path.join(args.results_dir, '*_abs_stats.txt')))
    pre_files = sorted(glob.glob(
        os.path.join(args.results_dir, '*_pre_stats.txt')))

    if not abs_files and not pre_files:
        print(f"No stats files found in {args.results_dir}", file=sys.stderr)
        sys.exit(1)

    # Collect counts by sample
    data = defaultdict(lambda: {'abs': 0, 'pre': 0})

    for f in abs_files:
        sample = os.path.basename(f).replace('_abs_stats.txt', '')
        data[sample]['abs'] = parse_bbduk_stats(f)

    for f in pre_files:
        sample = os.path.basename(f).replace('_pre_stats.txt', '')
        data[sample]['pre'] = parse_bbduk_stats(f)

    # Output
    output_path = args.output or os.path.join(
        args.results_dir, 'genotype_results.tsv')

    with open(output_path, 'w') as out:
        out.write("sample\tgenotype\tabs_reads\tpre_reads\ttotal\tnotes\n")

        for sample in sorted(data):
            abs_c = data[sample]['abs']
            pre_c = data[sample]['pre']
            total = abs_c + pre_c
            genotype, notes = call_genotype(
                abs_c, pre_c, args.min_depth, args.het_min, args.het_max)
            out.write(f"{sample}\t{genotype}\t{abs_c}\t{pre_c}\t{total}\t{notes}\n")

    print(f"Wrote genotypes for {len(data)} samples to {output_path}")
    print()

    # Print results
    with open(output_path) as f:
        for line in f:
            print(line, end='')


if __name__ == '__main__':
    main()

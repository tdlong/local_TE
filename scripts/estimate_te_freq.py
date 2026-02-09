#!/usr/bin/env python3
"""
estimate_te_freq.py

Estimate TE insertion frequencies from competitive alignment SAM files.

Reads are aligned to a junctions.fa reference containing paired sequences:
  - {id}_Abs: absence allele (wild-type reference, no TE)
  - {id}_Pre: presence allele (junction contig spanning the TE insertion)

A read is counted only if it spans the junction midpoint (position 50 in the
100bp reference) with at least 10bp anchor on each side. This discards reads
that map only to the shared 50bp flank, which would be ambiguous.

Usage:
    python estimate_te_freq.py <sam_file> [sample_name]

Output (TSV to stdout):
    sample  junction_id  abs_count  pre_count  frequency
"""

import sys
import re
from collections import defaultdict


def parse_cigar_reflen(cigar):
    """Calculate reference length consumed by a CIGAR string."""
    total = 0
    for match in re.finditer(r'(\d+)([MIDNSHP=X])', cigar):
        length = int(match.group(1))
        op = match.group(2)
        if op in 'MDN=X':  # ops that consume reference
            total += length
    return total


def get_counts(sam_file):
    """Parse SAM and count spanning reads per Abs/Pre junction pair."""
    counts = defaultdict(lambda: {"Abs": 0, "Pre": 0})

    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue

            cols = line.split('\t')
            if len(cols) < 11:
                continue

            flag = int(cols[1])
            refname = cols[2]
            pos = int(cols[3])
            cigar = cols[5]

            # Skip unmapped reads
            if flag & 4:
                continue
            if refname == '*':
                continue

            # Determine allele (Abs or Pre)
            if refname.endswith('_Abs'):
                junction_id = refname[:-4]
                allele = 'Abs'
            elif refname.endswith('_Pre'):
                junction_id = refname[:-4]
                allele = 'Pre'
            else:
                continue

            # Calculate alignment end on reference using CIGAR
            ref_len = parse_cigar_reflen(cigar)
            end_pos = pos + ref_len

            # SPANNING FILTER:
            # Read must start at/before position 40 and end at/after position 60.
            # This ensures >= 10bp anchor on both sides of the junction midpoint (50).
            # Reads mapping only to the shared 50bp flank are discarded.
            if pos <= 40 and end_pos >= 60:
                counts[junction_id][allele] += 1

    return counts


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    sam_file = sys.argv[1]
    sample_name = sys.argv[2] if len(sys.argv) > 2 else "unknown"

    results = get_counts(sam_file)

    # Output header + results as TSV
    print("sample\tjunction_id\tabs_count\tpre_count\tfrequency")

    for junction_id in sorted(results.keys()):
        abs_c = results[junction_id]["Abs"]
        pre_c = results[junction_id]["Pre"]
        total = abs_c + pre_c
        freq = pre_c / total if total > 0 else 0.0
        print(f"{sample_name}\t{junction_id}\t{abs_c}\t{pre_c}\t{freq:.4f}")

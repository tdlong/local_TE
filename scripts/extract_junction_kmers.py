#!/usr/bin/env python3
"""
extract_junction_kmers.py

Extract allele-specific k-mers spanning junction points from junctions.fa.
For each Abs/Pre pair, extracts k=31 k-mers that cross the junction
(position 50 in 100bp sequences) and keeps only allele-specific ones
(present in one allele but not the other).

Writes abs_kmers.fa and pre_kmers.fa for use as BBDuk references.

Usage:
    python extract_junction_kmers.py <junctions.fa> [output_dir]
"""

import sys
import os
from Bio import SeqIO


def extract_junction_spanning_kmers(abs_seq, pre_seq, junction_pos, k=31):
    """
    Extract k-mers that span the junction point.
    Only returns allele-specific k-mers (set difference).

    junction_pos: 0-based position where junction occurs.
    A k-mer "spans" the junction if it includes bases from both sides.
    """
    abs_junction_kmers = set()
    pre_junction_kmers = set()

    for start in range(max(0, junction_pos - k + 1),
                       min(junction_pos + 1, len(abs_seq) - k + 1)):
        abs_kmer = abs_seq[start:start + k]
        pre_kmer = pre_seq[start:start + k]
        if len(abs_kmer) == k:
            abs_junction_kmers.add(abs_kmer)
        if len(pre_kmer) == k:
            pre_junction_kmers.add(pre_kmer)

    abs_specific = abs_junction_kmers - pre_junction_kmers
    pre_specific = pre_junction_kmers - abs_junction_kmers

    return abs_specific, pre_specific


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    junctions_fa = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "."

    os.makedirs(output_dir, exist_ok=True)

    records = list(SeqIO.parse(junctions_fa, 'fasta'))
    if len(records) == 0:
        print(f"Error: no sequences found in {junctions_fa}")
        sys.exit(1)
    if len(records) % 2 != 0:
        print(f"Error: expected even number of records (Abs/Pre pairs), got {len(records)}")
        sys.exit(1)

    abs_fa = os.path.join(output_dir, "abs_kmers.fa")
    pre_fa = os.path.join(output_dir, "pre_kmers.fa")

    total_abs = 0
    total_pre = 0

    with open(abs_fa, 'w') as af, open(pre_fa, 'w') as pf:
        for i in range(0, len(records), 2):
            abs_rec = records[i]
            pre_rec = records[i + 1]

            if not abs_rec.id.endswith('_Abs') or not pre_rec.id.endswith('_Pre'):
                print(f"Error: expected _Abs/_Pre pair, got {abs_rec.id} and {pre_rec.id}")
                sys.exit(1)

            junction_id = abs_rec.id[:-4]  # strip _Abs
            abs_seq = str(abs_rec.seq)
            pre_seq = str(pre_rec.seq)

            # Junction is at position 50 in these 100bp sequences
            junction_pos = 50

            abs_specific, pre_specific = extract_junction_spanning_kmers(
                abs_seq, pre_seq, junction_pos, k=31)

            print(f"{junction_id}: {len(abs_specific)} abs-specific, "
                  f"{len(pre_specific)} pre-specific k-mers")

            for j, kmer in enumerate(sorted(abs_specific)):
                af.write(f">{junction_id}__kmer_{j:03d}\n{kmer}\n")
                total_abs += 1

            for j, kmer in enumerate(sorted(pre_specific)):
                pf.write(f">{junction_id}__kmer_{j:03d}\n{kmer}\n")
                total_pre += 1

    print(f"\nWrote {total_abs} absence k-mers to {abs_fa}")
    print(f"Wrote {total_pre} presence k-mers to {pre_fa}")


if __name__ == '__main__':
    main()

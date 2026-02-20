#!/usr/bin/env python3
"""
build_junctions_ref.py

Collect junction_*.fasta files from specified directories, extract Abs/Pre
sequence pairs, validate, deduplicate, and write junctions.fa for Phase 2.

Usage:
    python build_junctions_ref.py <output_junctions.fa> <dir1> [dir2 dir3 ...]

Each directory is searched for junction_*.fasta files (non-recursive).
Only directories explicitly listed are searched — no wildcards over
parent directories — so the caller controls exactly which results
are included.

Input:
    junction_*.fasta files with 4 sequences each:
      1. >WT_REF[start:end] insertion=chr3L:8711446 te=FBte0000626 type=right
      2. >REF
      3. >junction_contig_name
      4. >TE_name

Output:
    - junctions.fa:  Competitive alignment reference with Abs/Pre pairs
    - junctions_metadata.tsv: Metadata for each junction
"""

import os
import sys
import glob
from Bio import SeqIO


def parse_wt_ref_header(description):
    """Parse metadata from WT_REF header line.

    Format: WT_REF[start:end] insertion=chr3L:8711446 te=FBte0000626 type=right
    """
    parts = description.split()
    metadata = {}
    for part in parts[1:]:
        if '=' in part:
            key, val = part.split('=', 1)
            metadata[key] = val
    return metadata


def main():
    if len(sys.argv) < 3:
        print("Usage: build_junctions_ref.py <output.fa> <dir1> [dir2 ...]")
        sys.exit(1)

    output_fa = sys.argv[1]
    search_dirs = sys.argv[2:]

    # Derive metadata TSV name
    base, ext = os.path.splitext(output_fa)
    output_tsv = base + "_metadata.tsv"

    # Find junction files in the specified directories only
    junction_files = []
    for d in search_dirs:
        found = sorted(glob.glob(os.path.join(d, 'junction_*.fasta')))
        junction_files.extend(found)

    if not junction_files:
        print(f"No junction files found in: {' '.join(search_dirs)}")
        sys.exit(1)

    print(f"Found {len(junction_files)} junction files in "
          f"{len(search_dirs)} directories")

    # Collect all junctions
    # key: (insertion_pos, type) -> data dict
    junctions = {}

    for jfile in junction_files:
        records = list(SeqIO.parse(jfile, 'fasta'))
        if len(records) < 4:
            print(f"  Skipping {jfile}: expected 4 sequences, got {len(records)}")
            continue

        wt_ref_rec = records[0]   # WT_REF (absence allele)
        # records[1] is REF (duplicate)
        junc_rec = records[2]     # Junction contig (presence allele)
        te_rec = records[3]       # TE sequence

        # Parse metadata from header
        metadata = parse_wt_ref_header(wt_ref_rec.description)

        if 'insertion' not in metadata or 'te' not in metadata:
            print(f"  Skipping {jfile}: missing metadata in WT_REF header")
            print(f"    Header: {wt_ref_rec.description}")
            continue

        insertion = metadata['insertion']
        te_name = metadata['te']
        jtype = metadata.get('type', 'unknown')

        # Strip dashes from sequences to get clean Abs/Pre
        abs_seq = str(wt_ref_rec.seq).replace('-', '')
        pre_seq = str(junc_rec.seq).replace('-', '')

        # Validate: reference half of Pre must match Abs
        half = len(abs_seq) // 2
        if jtype == 'right':
            abs_ref = abs_seq[half:].upper()
            pre_ref = pre_seq[half:].upper()
        else:
            abs_ref = abs_seq[:half].upper()
            pre_ref = pre_seq[:half].upper()

        mismatches = sum(1 for a, b in zip(abs_ref, pre_ref) if a != b)
        if mismatches > 2:
            print(f"  REJECT {jfile}: {insertion} {te_name} ({jtype}) — "
                  f"ref half has {mismatches}/{half} mismatches vs Abs")
            continue

        # Dedup key: insertion position + junction type (NOT TE name).
        # Different TE database entries can match the same physical
        # insertion — they should not produce separate junctions.
        key = (insertion, jtype)

        if key in junctions:
            existing_te = junctions[key]['te_name']
            print(f"  Duplicate: {insertion} {te_name} ({jtype}) — "
                  f"same position as {existing_te}, keeping first")
            continue

        junctions[key] = {
            'insertion': insertion,
            'te_name': te_name,
            'type': jtype,
            'abs_seq': abs_seq,
            'pre_seq': pre_seq,
            'source_file': os.path.basename(jfile),
            'source_dir': os.path.basename(os.path.dirname(jfile)),
            'contig': junc_rec.id,
        }

        print(f"  {insertion} {te_name} ({jtype}) - "
              f"Abs:{len(abs_seq)}bp Pre:{len(pre_seq)}bp")

    if not junctions:
        print("No valid junctions found!")
        sys.exit(1)

    # Write junctions.fa and metadata TSV
    with open(output_fa, 'w') as fa, open(output_tsv, 'w') as tsv:
        tsv.write("junction_id\tinsertion\tte_name\ttype\tcontig\t"
                  "source_dir\tabs_len\tpre_len\n")

        for (insertion, jtype), data in sorted(junctions.items()):
            te_name = data['te_name']
            safe_insertion = insertion.replace(':', '_')
            junction_id = f"{safe_insertion}_{te_name}_{jtype}"

            fa.write(f">{junction_id}_Abs\n{data['abs_seq']}\n")
            fa.write(f">{junction_id}_Pre\n{data['pre_seq']}\n")

            tsv.write(f"{junction_id}\t{data['insertion']}\t{data['te_name']}\t"
                      f"{data['type']}\t{data['contig']}\t{data['source_dir']}\t"
                      f"{len(data['abs_seq'])}\t{len(data['pre_seq'])}\n")

    print(f"\nWrote {len(junctions)} junction pairs to {output_fa}")
    print(f"Metadata: {output_tsv}")


if __name__ == "__main__":
    main()

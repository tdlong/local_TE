#!/usr/bin/env python3
"""
build_junctions_from_reads.py

Read-level TE junction discovery. Replaces SPAdes assembly (Sections 3-7 of
run_te_assembly.sh). BLAST individual reads directly against the TE database
and the reference region, find junction reads (hits to both), cluster by
inferred insertion position, and build consensus junction sequences.

Writes one junction_<type>_<N>.fasta per cluster (compatible with
build_junctions_ref.py 4-record format).

Usage:
    python build_junctions_from_reads.py \\
        --outdir  <outdir>      \\  # contains R1.fq R2.fq singles.fq region.fasta
        --te-fasta <te.fasta>   \\
        --region  <chr:start-end> \\
        [--min-support 3]       \\
        [--snp-min-freq 0.15]
"""

import argparse
import os
import subprocess
import sys
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq


# ---------------------------------------------------------------------------
# IUPAC encoding helpers
# ---------------------------------------------------------------------------
IUPAC_FROM_BASES = {
    frozenset('A'):    'A',
    frozenset('C'):    'C',
    frozenset('G'):    'G',
    frozenset('T'):    'T',
    frozenset('AG'):   'R',
    frozenset('CT'):   'Y',
    frozenset('GC'):   'S',
    frozenset('AT'):   'W',
    frozenset('GT'):   'K',
    frozenset('AC'):   'M',
    frozenset('CGT'):  'B',
    frozenset('AGT'):  'D',
    frozenset('ACT'):  'H',
    frozenset('ACG'):  'V',
    frozenset('ACGT'): 'N',
}


def bases_to_iupac(counts, snp_min_freq):
    """Return IUPAC code for a column given base counts {A:n, C:n, G:n, T:n}.

    Major allele is always included. Any additional allele with frequency >=
    snp_min_freq is also included (IUPAC ambiguity).
    Returns '-' if no coverage.
    """
    total = sum(counts.values())
    if total == 0:
        return 'N'
    present = {b for b, n in counts.items() if n / total >= snp_min_freq}
    key = frozenset(present)
    return IUPAC_FROM_BASES.get(key, 'N')


# ---------------------------------------------------------------------------
# FASTQ → FASTA
# ---------------------------------------------------------------------------
def fastq_to_fasta(fq_paths, out_fa):
    """Concatenate one or more FASTQ files into a single FASTA."""
    with open(out_fa, 'w') as fh:
        for fq in fq_paths:
            if not os.path.exists(fq) or os.path.getsize(fq) == 0:
                continue
            with open(fq) as f:
                while True:
                    header = f.readline().rstrip()
                    seq    = f.readline().rstrip()
                    plus   = f.readline()
                    qual   = f.readline()
                    if not header:
                        break
                    name = header[1:].split()[0]   # strip @, take first token
                    fh.write(f">{name}\n{seq}\n")


# ---------------------------------------------------------------------------
# BLAST helpers
# ---------------------------------------------------------------------------
BLAST_FMT = "6 qseqid sseqid qstart qend sstart send sstrand"


def run_blastn(query, subject, out_tsv, evalue="1e-5"):
    cmd = [
        "blastn",
        "-query", query,
        "-subject", subject,
        "-outfmt", BLAST_FMT,
        "-evalue", evalue,
        "-out", out_tsv,
    ]
    subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL)


def parse_best_hits(tsv_path):
    """Return {read_id: best_hit_dict} keeping hit with longest aligned span."""
    best = {}
    if not os.path.exists(tsv_path):
        return best
    with open(tsv_path) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 7:
                continue
            qid, sid, qstart, qend, sstart, send, sstrand = parts[:7]
            qstart, qend = int(qstart), int(qend)
            sstart, send = int(sstart), int(send)
            span = qend - qstart + 1
            if qid not in best or span > best[qid]['span']:
                best[qid] = {
                    'subject':  sid,
                    'qstart':   qstart,   # 1-based
                    'qend':     qend,
                    'sstart':   sstart,
                    'send':     send,
                    'sstrand':  sstrand,
                    'span':     span,
                }
    return best


# ---------------------------------------------------------------------------
# Reverse complement
# ---------------------------------------------------------------------------
def revcomp(seq):
    return str(Seq(seq).reverse_complement())


# ---------------------------------------------------------------------------
# Consensus building with IUPAC SNP encoding
# ---------------------------------------------------------------------------
def build_consensus(seqs_at_anchor, half=50, snp_min_freq=0.15):
    """
    Build a 2*half bp consensus centred at column `half` (0-based).

    seqs_at_anchor: list of (seq, anchor_col) where anchor_col is the
    position within `seq` that corresponds to the junction (will be placed
    at column `half` of the output).

    Returns None if no reads cover position `half`.
    """
    width = 2 * half
    counts = [defaultdict(int) for _ in range(width)]

    for seq, anchor in seqs_at_anchor:
        # offset so that seq[anchor] → column half
        col_offset = half - anchor
        for i, base in enumerate(seq.upper()):
            col = i + col_offset
            if 0 <= col < width and base in 'ACGT':
                counts[col][base] += 1

    # Check coverage at the junction column itself
    if sum(counts[half].values()) == 0:
        return None

    consensus = ''.join(bases_to_iupac(counts[col], snp_min_freq)
                        for col in range(width))
    return consensus


# ---------------------------------------------------------------------------
# Load reads.fasta into dict
# ---------------------------------------------------------------------------
def load_reads(fa_path):
    reads = {}
    for rec in SeqIO.parse(fa_path, 'fasta'):
        reads[rec.id] = str(rec.seq).upper()
    return reads


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--outdir',      required=True)
    ap.add_argument('--te-fasta',    required=True)
    ap.add_argument('--region',      required=True,
                    help='Region string, e.g. chr3L:8710861-8744900')
    ap.add_argument('--min-support', type=int, default=3)
    ap.add_argument('--snp-min-freq', type=float, default=0.15)
    args = ap.parse_args()

    outdir       = args.outdir
    te_fasta     = args.te_fasta
    region       = args.region
    min_support  = args.min_support
    snp_min_freq = args.snp_min_freq

    # -----------------------------------------------------------------------
    # Parse region coordinates
    # -----------------------------------------------------------------------
    chrom        = region.split(':')[0]
    coords       = region.split(':')[1]
    region_start = int(coords.split('-')[0])   # 1-based genomic start

    # -----------------------------------------------------------------------
    # Step 1: Build reads.fasta (R1 + R2 + singles)
    # -----------------------------------------------------------------------
    reads_fa = os.path.join(outdir, 'reads.fasta')
    fq_files = [
        os.path.join(outdir, 'R1.fq'),
        os.path.join(outdir, 'R2.fq'),
        os.path.join(outdir, 'singles.fq'),
    ]
    print("Step 1: Converting FASTQ to FASTA")
    fastq_to_fasta(fq_files, reads_fa)
    n_reads = sum(1 for _ in open(reads_fa) if _.startswith('>'))
    print(f"  {n_reads} reads in {reads_fa}")

    if n_reads == 0:
        print("  No reads found — no junctions to discover")
        return

    # -----------------------------------------------------------------------
    # Step 2: BLAST reads vs TE database and reference region
    # -----------------------------------------------------------------------
    region_fa   = os.path.join(outdir, 'region.fasta')
    te_blast    = os.path.join(outdir, 'reads_vs_te.tsv')
    ref_blast   = os.path.join(outdir, 'reads_vs_ref.tsv')

    print("Step 2: BLAST reads vs TE database")
    run_blastn(reads_fa, te_fasta, te_blast)
    print(f"  TE hits: {sum(1 for _ in open(te_blast) if _.strip())}")

    print("Step 2: BLAST reads vs reference region")
    run_blastn(reads_fa, region_fa, ref_blast)
    print(f"  Ref hits: {sum(1 for _ in open(ref_blast) if _.strip())}")

    te_hits  = parse_best_hits(te_blast)
    ref_hits = parse_best_hits(ref_blast)

    # -----------------------------------------------------------------------
    # Step 3: Find junction reads (hits to both TE and ref)
    # -----------------------------------------------------------------------
    print("Step 3: Identifying junction reads")
    junction_reads = []   # list of dicts

    for read_id in set(te_hits) & set(ref_hits):
        te_h  = te_hits[read_id]
        ref_h = ref_hits[read_id]

        te_qstart  = te_h['qstart']   # 1-based position in read where TE starts
        ref_qstart = ref_h['qstart']  # 1-based position in read where ref starts

        # Infer insertion point in genomic coordinates
        if ref_h['sstrand'] == 'plus':
            ins_pos = ref_h['sstart'] + (te_qstart - ref_qstart)
            # junction type: ref alignment is to the left of TE hit in read?
            jtype = 'left' if ref_qstart < te_qstart else 'right'
        else:
            ins_pos = ref_h['sstart'] - (te_qstart - ref_qstart)
            jtype = 'left' if ref_qstart < te_qstart else 'right'

        # Convert to absolute genomic coordinate
        # (sstart in BLAST is 1-based offset within the extracted region)
        abs_ins_pos = region_start + ins_pos - 1

        junction_reads.append({
            'read_id':     read_id,
            'te_name':     te_h['subject'],
            'te_qstart':   te_qstart,    # 1-based
            'ref_qstart':  ref_qstart,   # 1-based
            'ref_strand':  ref_h['sstrand'],
            'ins_pos':     abs_ins_pos,
            'jtype':       jtype,
        })

    print(f"  Found {len(junction_reads)} junction reads")
    if not junction_reads:
        print("  No junction reads — no clusters to build")
        return

    # -----------------------------------------------------------------------
    # Step 4: Cluster by inferred insertion position (±20bp)
    # -----------------------------------------------------------------------
    print("Step 4: Clustering by insertion position")
    junction_reads.sort(key=lambda r: r['ins_pos'])

    clusters = []   # list of lists of reads
    current  = [junction_reads[0]]

    for jr in junction_reads[1:]:
        if abs(jr['ins_pos'] - current[0]['ins_pos']) <= 20:
            current.append(jr)
        else:
            clusters.append(current)
            current = [jr]
    clusters.append(current)

    print(f"  {len(clusters)} raw clusters "
          f"(≥{min_support} support required)")

    # -----------------------------------------------------------------------
    # Load reads into memory (needed for consensus building)
    # -----------------------------------------------------------------------
    read_seqs = load_reads(reads_fa)

    # -----------------------------------------------------------------------
    # Load reference region and TE sequences for Abs and TE records
    # -----------------------------------------------------------------------
    region_seq_dict = {rec.id: str(rec.seq).upper()
                       for rec in SeqIO.parse(region_fa, 'fasta')}
    # region_fa has one record; take its sequence
    region_seq = next(iter(region_seq_dict.values()))

    te_seq_dict = {rec.id: str(rec.seq).upper()
                   for rec in SeqIO.parse(te_fasta, 'fasta')}

    # -----------------------------------------------------------------------
    # Step 5–8: Build consensus and write junction FASTA files
    # -----------------------------------------------------------------------
    half = 50   # junction centred at column 50 (0-based) of a 100bp window

    n_written = 0
    for cluster_id, cluster in enumerate(clusters):
        # Separate LEFT and RIGHT junction reads
        left_reads  = [r for r in cluster if r['jtype'] == 'left']
        right_reads = [r for r in cluster if r['jtype'] == 'right']

        # Determine representative insertion position and TE name (majority)
        all_te_names = [r['te_name'] for r in cluster]
        te_name = max(set(all_te_names), key=all_te_names.count)

        med_pos = sorted(r['ins_pos'] for r in cluster)[len(cluster) // 2]

        # Build consensus for whichever side has enough support
        for side, side_reads, min_side in [('left',  left_reads,  2),
                                            ('right', right_reads, 2)]:
            if len(side_reads) < min(min_side, min_support):
                continue

            # Step 5: Anchor reads at the TE boundary
            seqs_at_anchor = []
            for jr in side_reads:
                seq = read_seqs.get(jr['read_id'])
                if seq is None:
                    continue
                # Reverse-complement if ref hit was on minus strand
                if jr['ref_strand'] == 'minus':
                    seq = revcomp(seq)
                if side == 'left':
                    # anchor = te_qstart - 1 (convert to 0-based)
                    anchor = jr['te_qstart'] - 1
                else:  # right: ref begins after TE in the read
                    anchor = jr['ref_qstart'] - 1
                seqs_at_anchor.append((seq, anchor))

            pre_seq = build_consensus(seqs_at_anchor, half=half,
                                      snp_min_freq=snp_min_freq)
            if pre_seq is None:
                print(f"  Cluster {cluster_id} {side}: consensus failed "
                      f"(no coverage at junction)")
                continue

            # Step 6: Abs sequence — reference genome centred at insertion
            # Convert med_pos (1-based genomic) to 0-based index in region_seq
            ins_idx = med_pos - region_start   # 0-based offset into region_seq
            abs_start = ins_idx - half
            abs_end   = ins_idx + half
            if abs_start < 0 or abs_end > len(region_seq):
                print(f"  Cluster {cluster_id}: insertion near region boundary, "
                      f"skipping")
                continue
            abs_seq = region_seq[abs_start:abs_end]

            # Step 7: TE sequence for visualization record
            te_full = te_seq_dict.get(te_name, '')
            if side == 'left':
                te_sub = te_full[:100]
            else:
                te_sub = te_full[-100:] if len(te_full) >= 100 else te_full

            # Step 8: Write junction FASTA (4-record format)
            # Genomic window used for Abs
            abs_gstart = region_start + abs_start        # 1-based
            abs_gend   = abs_gstart + len(abs_seq) - 1

            out_name = os.path.join(
                outdir, f"junction_{side}_{cluster_id:03d}.fasta")

            with open(out_name, 'w') as fh:
                fh.write(
                    f">WT_REF[{abs_gstart}:{abs_gend}] "
                    f"insertion={chrom}:{med_pos} "
                    f"te={te_name} "
                    f"type={side}\n"
                )
                fh.write(abs_seq + '\n')
                fh.write(">REF\n")
                fh.write(abs_seq + '\n')
                fh.write(f">reads_consensus_{side}_{cluster_id}\n")
                fh.write(pre_seq + '\n')
                fh.write(f">{te_name}\n")
                fh.write(te_sub + '\n')

            n_support = len(side_reads)
            print(f"  Wrote {out_name} "
                  f"(ins={chrom}:{med_pos}, te={te_name}, "
                  f"n={n_support})")
            n_written += 1

    print(f"\nPhase 1 (read-level discovery): {n_written} junction files written")


if __name__ == '__main__':
    main()

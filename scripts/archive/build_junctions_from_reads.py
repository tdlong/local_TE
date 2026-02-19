#!/usr/bin/env python3
"""
build_junctions_from_reads.py

Read-level TE junction discovery. Replaces SPAdes assembly (Sections 3-7 of
run_te_assembly.sh). BLAST individual reads directly against the TE database
and the reference region, identify junction reads (hits to both), cluster by
inferred insertion position, and build consensus junction sequences.

Key design: reads are anchored at the GENOMICALLY-DERIVED insertion point
(cluster median position back-projected through the reference alignment), not
at te_qstart.  BLAST can place te_qstart ±3-5bp differently across reads from
the same junction; that jitter makes every column look like a SNP with only
2-3 reads.  The reference alignment is much more precisely constrained and
gives a consistent anchor across all reads in a cluster.

Writes one junction_<type>_<N>.fasta per cluster (compatible with
build_junctions_ref.py 4-record format).

Usage:
    python build_junctions_from_reads.py \\
        --outdir  <outdir>        \\   # contains R1.fq R2.fq singles.fq region.fasta
        --te-fasta <te.fasta>     \\
        --region  <chr:start-end> \\
        [--min-support 3]         \\
        [--snp-min-freq 0.15]     \\
        [--min-flank 10]
"""

import argparse
import os
import subprocess
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq


# ---------------------------------------------------------------------------
# IUPAC encoding
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


def bases_to_iupac(counts, snp_min_freq, snp_min_count=2):
    """
    Return IUPAC code for a base-count column.

    The major allele is always represented.  A minor allele is included only
    when it meets BOTH thresholds:
      - count >= snp_min_count (default 2) — prevents single sequencing errors
      - frequency >= snp_min_freq — prevents low-frequency noise

    With snp_min_count=2 you need at least 3 reads to encode any IUPAC
    ambiguity (2 showing the minor allele + 1 or more showing the major).
    This keeps the consensus clean at low depth.
    """
    total = sum(counts.values())
    if total == 0:
        return 'N'
    major = max(counts, key=counts.get)
    present = {major}
    for b, n in counts.items():
        if b != major and n >= snp_min_count and n / total >= snp_min_freq:
            present.add(b)
    return IUPAC_FROM_BASES.get(frozenset(present), major)


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def revcomp(seq):
    return str(Seq(seq).reverse_complement())


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
                    f.readline()   # +
                    f.readline()   # qual
                    if not header:
                        break
                    name = header[1:].split()[0]
                    fh.write(f">{name}\n{seq}\n")


# ---------------------------------------------------------------------------
# BLAST
# ---------------------------------------------------------------------------
BLAST_FMT = "6 qseqid sseqid qstart qend sstart send sstrand"


def run_blastn(query, subject, out_tsv, evalue="1e-5"):
    cmd = ["blastn", "-query", query, "-subject", subject,
           "-outfmt", BLAST_FMT, "-evalue", evalue, "-out", out_tsv]
    subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL)


def parse_best_hits(tsv_path):
    """Return {read_id: best_hit_dict} keeping hit with longest aligned span."""
    best = {}
    if not os.path.exists(tsv_path):
        return best
    with open(tsv_path) as f:
        for line in f:
            parts = line.rstrip().split('\t')
            if len(parts) < 7:
                continue
            qid, sid, qstart, qend, sstart, send, sstrand = parts[:7]
            qstart, qend = int(qstart), int(qend)
            sstart, send = int(sstart), int(send)
            span = qend - qstart + 1
            if qid not in best or span > best[qid]['span']:
                best[qid] = {
                    'subject': sid,
                    'qstart': qstart, 'qend': qend,
                    'sstart': sstart, 'send': send,
                    'sstrand': sstrand, 'span': span,
                }
    return best


# ---------------------------------------------------------------------------
# Anchor computation
# ---------------------------------------------------------------------------
def compute_anchor(jr, cluster_ins_region, read_seqs, min_flank):
    """
    Back-project the cluster's consensus insertion position (1-based in
    region.fasta) through the read's reference alignment to get the 0-based
    position in the read that corresponds to the insertion point.

    For minus-strand reads the read is reverse-complemented so that the
    returned sequence is in forward genomic orientation.

    Returns (oriented_seq, anchor_0based), or None if the insertion falls
    within min_flank bases of either end of the read.

    Coordinate math:
      plus strand:  subject_pos = ref_sstart + (query_pos - ref_qstart)
                    → query_pos = ref_qstart + (subject_pos - ref_sstart)
      minus strand: subject_pos = ref_sstart - (query_pos - ref_qstart)
                    → query_pos = ref_qstart + (ref_sstart - subject_pos)
      After RC'ing a minus-strand read, position p becomes (len-1-p).
    """
    seq = read_seqs.get(jr['read_id'])
    if not seq:
        return None

    ref_qstart_0 = jr['ref_qstart'] - 1   # 0-based in original read
    ref_sstart   = jr['ref_sstart']        # 1-based in region.fasta

    if jr['ref_strand'] == 'plus':
        anchor_0 = ref_qstart_0 + (cluster_ins_region - ref_sstart)
    else:
        anchor_0 = ref_qstart_0 + (ref_sstart - cluster_ins_region)
        seq      = revcomp(seq)
        anchor_0 = len(seq) - 1 - anchor_0

    anchor_0 = int(round(anchor_0))

    if anchor_0 < min_flank or anchor_0 >= len(seq) - min_flank:
        return None

    return seq, anchor_0


def compute_te_only_anchor(read_id, te_hit, te_junc_coord, canonical_te_strand,
                            read_seqs, side, min_flank):
    """
    Place a TE-only read (no reference BLAST hit) in the junction consensus
    window by projecting the TE junction coordinate into the read.

    The TE junction coordinate (te_junc_coord) is the position in the TE
    database that corresponds to the insertion boundary:
      RIGHT junction → te_send   (end of TE that meets the reference right flank)
      LEFT  junction → te_sstart (start of TE that meets the reference left flank)

    Coordinate math (symmetric to compute_anchor but using TE alignment):
      plus  strand: anchor = qstart + (te_junc_coord - sstart)
      minus strand: anchor = qstart + (sstart - te_junc_coord)
        (for minus strand BLAST hits sstart > send; the alignment runs
         from sstart down to send as read position increases)

    The read is then oriented to match the junction reads (same canonical TE
    strand).  Only the TE side of the anchor is used for the consensus — the
    post-anchor sequence (reference side) is excluded by build_consensus
    naturally if the read ends at the anchor.

    Returns (oriented_seq, anchor_0based), or None if the read cannot
    contribute at least min_flank bases on the TE side of the anchor.
    """
    seq = read_seqs.get(read_id)
    if not seq:
        return None

    qstart  = te_hit['qstart']
    sstart  = te_hit['sstart']
    tstrand = te_hit['sstrand']

    # 1-based anchor in the read
    if tstrand == 'plus':
        anchor_1 = qstart + (te_junc_coord - sstart)
    else:   # minus strand: sstart > te_junc_coord
        anchor_1 = qstart + (sstart - te_junc_coord)

    anchor_0 = anchor_1 - 1   # 0-based

    # Orient consistently with the junction reads
    if tstrand != canonical_te_strand:
        seq      = revcomp(seq)
        anchor_0 = len(seq) - 1 - anchor_0

    if anchor_0 < 0 or anchor_0 >= len(seq):
        return None

    # Require min_flank bases on the TE side of the anchor
    if side == 'right' and anchor_0 < min_flank:
        return None
    if side == 'left'  and (len(seq) - anchor_0) < min_flank:
        return None

    return seq, anchor_0


def canonical_te_half(te_seq, te_junc_coord, canon_te_strand, side, half=50):
    """
    Extract `half` bp of canonical TE sequence to fill the TE half of the
    junction consensus.

    Coordinate logic (te_junc_coord is 1-based from BLAST):
      RIGHT, plus  → te_seq[te_junc_coord-half : te_junc_coord]   (forward)
      RIGHT, minus → RC(te_seq[te_junc_coord-1 : te_junc_coord-1+half])
      LEFT,  plus  → te_seq[te_junc_coord-1    : te_junc_coord-1+half] (forward)
      LEFT,  minus → RC(te_seq[te_junc_coord-half : te_junc_coord])

    Returns None if the window is out of bounds or te_seq is empty.
    """
    if not te_seq:
        return None
    if side == 'right' and canon_te_strand == 'plus':
        s, e, do_rc = te_junc_coord - half, te_junc_coord, False
    elif side == 'right' and canon_te_strand == 'minus':
        s, e, do_rc = te_junc_coord - 1, te_junc_coord - 1 + half, True
    elif side == 'left' and canon_te_strand == 'plus':
        s, e, do_rc = te_junc_coord - 1, te_junc_coord - 1 + half, False
    else:  # left, minus
        s, e, do_rc = te_junc_coord - half, te_junc_coord, True
    if s < 0 or e > len(te_seq):
        return None
    fill = te_seq[s:e]
    return revcomp(fill) if do_rc else fill


# ---------------------------------------------------------------------------
# Consensus building
# ---------------------------------------------------------------------------
def build_consensus(anchored_reads, half=50, snp_min_freq=0.15):
    """
    Build a 2*half bp consensus with column `half` at the insertion point.

    anchored_reads: list of (oriented_seq, anchor_0based)
    Returns None if no reads cover the anchor column itself.
    """
    width  = 2 * half
    counts = [defaultdict(int) for _ in range(width)]

    for seq, anchor in anchored_reads:
        col_offset = half - anchor
        for i, base in enumerate(seq.upper()):
            col = i + col_offset
            if 0 <= col < width and base in 'ACGT':
                counts[col][base] += 1

    if sum(counts[half].values()) == 0:
        return None

    return ''.join(bases_to_iupac(counts[col], snp_min_freq)
                   for col in range(width))


# ---------------------------------------------------------------------------
# Load helpers
# ---------------------------------------------------------------------------
def load_reads(fa_path):
    return {rec.id: str(rec.seq).upper()
            for rec in SeqIO.parse(fa_path, 'fasta')}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--outdir',       required=True)
    ap.add_argument('--te-fasta',     required=True)
    ap.add_argument('--region',       required=True,
                    help='e.g. chr3L:8710861-8744900')
    ap.add_argument('--min-support',  type=int,   default=3,
                    help='Min reads per cluster (default 3)')
    ap.add_argument('--snp-min-freq', type=float, default=0.15,
                    help='Minor-allele frequency threshold for IUPAC (default 0.15)')
    ap.add_argument('--min-flank',    type=int,   default=10,
                    help='Min bases each side of anchor in a read (default 10)')
    ap.add_argument('--max-snps',     type=int,   default=10,
                    help='Max IUPAC codes in Pre consensus before skipping '
                         '(default 10); high counts indicate merged clusters '
                         'from multiple nearby insertions')
    args = ap.parse_args()

    outdir       = args.outdir
    te_fasta     = args.te_fasta
    region       = args.region
    min_support  = args.min_support
    snp_min_freq = args.snp_min_freq
    min_flank    = args.min_flank
    max_snps     = args.max_snps

    chrom        = region.split(':')[0]
    region_start = int(region.split(':')[1].split('-')[0])   # 1-based genomic

    # ------------------------------------------------------------------
    # Step 1: FASTQ → FASTA
    # ------------------------------------------------------------------
    reads_fa = os.path.join(outdir, 'reads.fasta')
    fastq_to_fasta(
        [os.path.join(outdir, f) for f in ('R1.fq', 'R2.fq', 'singles.fq')],
        reads_fa)
    n_reads = sum(1 for line in open(reads_fa) if line.startswith('>'))
    print(f"Step 1: {n_reads} reads → {reads_fa}")
    if n_reads == 0:
        print("  No reads — nothing to discover")
        return

    # ------------------------------------------------------------------
    # Step 2: BLAST reads vs TE database and reference region
    # ------------------------------------------------------------------
    region_fa = os.path.join(outdir, 'region.fasta')
    te_blast  = os.path.join(outdir, 'reads_vs_te.tsv')
    ref_blast = os.path.join(outdir, 'reads_vs_ref.tsv')

    print("Step 2: BLAST reads vs TE database")
    run_blastn(reads_fa, te_fasta, te_blast)
    print(f"  TE hits: {sum(1 for l in open(te_blast) if l.strip())}")

    print("Step 2: BLAST reads vs reference region")
    run_blastn(reads_fa, region_fa, ref_blast)
    print(f"  Ref hits: {sum(1 for l in open(ref_blast) if l.strip())}")

    te_hits  = parse_best_hits(te_blast)
    ref_hits = parse_best_hits(ref_blast)

    # ------------------------------------------------------------------
    # Step 3: Identify junction reads; infer insertion position and type
    #
    # Junction type (LEFT vs RIGHT) describes the read's layout relative
    # to the TE insertion in GENOMIC forward orientation:
    #   LEFT:  [ref_sequence | TE_sequence]  in the read
    #   RIGHT: [TE_sequence  | ref_sequence] in the read
    #
    # For plus-strand ref hits:
    #   ref alignment precedes TE alignment in the read → LEFT
    # For minus-strand ref hits the read is in reverse orientation, so the
    # layout is flipped:
    #   ref alignment precedes TE in original read → RIGHT (genomically)
    #
    # Insertion position: use ref alignment ENDPOINTS, not te_qstart.
    # BLAST places te_qstart ±5-10bp variably; ref alignment endpoints are
    # precisely determined by sequence similarity and stable across reads.
    #
    #   LEFT junction  → insertion = TE-adjacent END of ref alignment
    #                  = max(sstart, send)  [high coord for both strands]
    #   RIGHT junction → insertion = TE-adjacent START of ref alignment
    #                  = min(sstart, send)  [low coord for both strands]
    #
    # Derivation:
    #   plus-strand LEFT:  ref ends at send (high coord) → ins = send ✓
    #   plus-strand RIGHT: ref starts at sstart (low coord) → ins = sstart ✓
    #   minus-strand LEFT: ref is adjacent to TE at sstart (high coord) → ins = sstart ✓
    #   minus-strand RIGHT: ref ends (at TE boundary) at send (low coord) → ins = send ✓
    # ------------------------------------------------------------------
    print("Step 3: Identifying junction reads")
    junction_reads = []

    for read_id in set(te_hits) & set(ref_hits):
        te_h  = te_hits[read_id]
        ref_h = ref_hits[read_id]

        te_qstart  = te_h['qstart']
        ref_qstart = ref_h['qstart']

        # Junction type: determined by relative read positions of TE and ref
        if ref_h['sstrand'] == 'plus':
            jtype = 'left' if ref_qstart < te_qstart else 'right'
        else:
            # Minus-strand: layout is flipped vs plus-strand
            jtype = 'right' if ref_qstart < te_qstart else 'left'

        # Insertion position from ref alignment endpoints (stable, no te_qstart)
        if jtype == 'left':
            ins_region = max(ref_h['sstart'], ref_h['send'])
        else:
            ins_region = min(ref_h['sstart'], ref_h['send'])

        junction_reads.append({
            'read_id':    read_id,
            'te_name':    te_h['subject'],
            'te_sstart':  te_h['sstart'],      # TE alignment start in TE database
            'te_send':    te_h['send'],         # TE alignment end   in TE database
            'te_strand':  te_h['sstrand'],
            'ref_qstart': ref_qstart,
            'ref_sstart': ref_h['sstart'],
            'ref_strand': ref_h['sstrand'],
            'ins_region': ins_region,          # 1-based in region.fasta
            'ins_pos':    region_start + ins_region - 1,  # 1-based genomic
            'jtype':      jtype,
        })

    print(f"  Found {len(junction_reads)} junction reads")
    if not junction_reads:
        return

    # ------------------------------------------------------------------
    # Step 4: Cluster by inferred insertion position (±20bp), then filter
    # ------------------------------------------------------------------
    print("Step 4: Clustering by insertion position")
    junction_reads.sort(key=lambda r: r['ins_pos'])
    clusters = []
    current  = [junction_reads[0]]
    for jr in junction_reads[1:]:
        if abs(jr['ins_pos'] - current[0]['ins_pos']) <= 20:
            current.append(jr)
        else:
            clusters.append(current)
            current = [jr]
    clusters.append(current)

    raw_n = len(clusters)
    clusters = [c for c in clusters if len(c) >= min_support]
    print(f"  {raw_n} raw clusters → {len(clusters)} with ≥{min_support} reads")

    if not clusters:
        return

    # ------------------------------------------------------------------
    # Build TE-only read index
    # Reads with a TE BLAST hit but no reference hit are the paired mates
    # that sit entirely within the TE.  They don't cross the junction so
    # they weren't classified as junction reads, but they do carry sequence
    # from deeper inside the TE — exactly what fills the N's in the TE half
    # of the junction consensus.
    # ------------------------------------------------------------------
    te_only_by_name = defaultdict(list)   # te_name → [(read_id, te_hit), …]
    for read_id, te_hit in te_hits.items():
        if read_id not in ref_hits:
            te_only_by_name[te_hit['subject']].append((read_id, te_hit))
    n_te_only_total = sum(len(v) for v in te_only_by_name.values())
    print(f"  TE-only reads available for TE-half extension: {n_te_only_total}")

    # ------------------------------------------------------------------
    # Load sequences
    # ------------------------------------------------------------------
    read_seqs   = load_reads(reads_fa)
    region_seq  = next(
        str(r.seq).upper() for r in SeqIO.parse(region_fa, 'fasta'))
    te_seq_dict = {r.id: str(r.seq).upper()
                   for r in SeqIO.parse(te_fasta, 'fasta')}

    # ------------------------------------------------------------------
    # Step 5–8: Build consensus and write junction FASTA files
    # ------------------------------------------------------------------
    half      = 50
    n_written = 0

    for cluster_id, cluster in enumerate(clusters):
        left_reads  = [r for r in cluster if r['jtype'] == 'left']
        right_reads = [r for r in cluster if r['jtype'] == 'right']

        # Majority TE name
        all_te  = [r['te_name'] for r in cluster]
        te_name = max(set(all_te), key=all_te.count)

        # Cluster consensus insertion position (median, 1-based in region.fasta)
        ins_list           = sorted(r['ins_region'] for r in cluster)
        cluster_ins_region = ins_list[len(ins_list) // 2]
        abs_ins_pos        = region_start + cluster_ins_region - 1   # genomic

        for side, side_reads in [('left', left_reads), ('right', right_reads)]:
            if len(side_reads) < 2:
                continue

            # Anchor each read at the genomically-derived insertion point
            anchored = []
            for jr in side_reads:
                result = compute_anchor(jr, cluster_ins_region,
                                        read_seqs, min_flank)
                if result is not None:
                    anchored.append(result)

            if len(anchored) < 2:
                continue

            n_junc_anchored = len(anchored)

            # ----------------------------------------------------------------
            # Extend TE half with TE-only reads
            #
            # The TE junction coordinate is the position in the TE database
            # sequence that corresponds to the insertion boundary:
            #   RIGHT → te_send   (TE's right edge that meets the ref flank)
            #   LEFT  → te_sstart (TE's left edge that meets the ref flank)
            #
            # After compute_anchor potentially RC's a junction read (when
            # ref_strand==minus), the effective TE strand is flipped.
            # canonical_te_strand is the TE strand as it appears in the
            # already-oriented junction reads; TE-only reads are oriented
            # to match.
            # ----------------------------------------------------------------
            te_junc_coords  = []
            te_strand_votes = []
            for jr in side_reads:
                te_junc_coords.append(
                    jr['te_send'] if side == 'right' else jr['te_sstart'])
                strand = jr['te_strand']
                if jr['ref_strand'] == 'minus':
                    strand = 'minus' if strand == 'plus' else 'plus'
                te_strand_votes.append(strand)

            te_junc_coord   = sorted(te_junc_coords)[len(te_junc_coords) // 2]
            canon_te_strand = (
                'plus' if te_strand_votes.count('plus')
                          >= te_strand_votes.count('minus') else 'minus')

            n_ext = 0
            for read_id, te_hit in te_only_by_name.get(te_name, []):
                te_lo = min(te_hit['sstart'], te_hit['send'])
                te_hi = max(te_hit['sstart'], te_hit['send'])
                if not (te_lo <= te_junc_coord <= te_hi):
                    continue
                result = compute_te_only_anchor(
                    read_id, te_hit, te_junc_coord, canon_te_strand,
                    read_seqs, side, min_flank)
                if result is not None:
                    anchored.append(result)
                    n_ext += 1

            if n_ext:
                print(f"    extended cluster {cluster_id} {side} with "
                      f"{n_ext} TE-only reads "
                      f"({n_junc_anchored} junction + {n_ext} TE-only = "
                      f"{len(anchored)} total)")

            pre_seq = build_consensus(anchored, half=half,
                                      snp_min_freq=snp_min_freq)
            if pre_seq is None:
                print(f"  Cluster {cluster_id} {side}: no coverage at anchor — skipped")
                continue

            # Abs: reference genome slice centred at insertion (no TE)
            ins_idx   = cluster_ins_region - 1   # 0-based in region_seq
            abs_start = ins_idx - half
            abs_end   = ins_idx + half
            if abs_start < 0 or abs_end > len(region_seq):
                print(f"  Cluster {cluster_id}: too close to region boundary — skipped")
                continue
            abs_seq = region_seq[abs_start:abs_end]

            # Load full TE sequence for canonical fill and visualization
            te_full = te_seq_dict.get(te_name, '')
            te_sub  = te_full[:100] if side == 'left' else te_full[-100:]

            # TE coverage from reads BEFORE canonical fill
            te_half_pre = pre_seq[half:] if side == 'left' else pre_seq[:half]
            read_te_cov = sum(1 for c in te_half_pre if c != 'N')

            # Fill N positions in TE half with canonical TE sequence
            n_canon_fill = 0
            canon_fill = canonical_te_half(
                te_full, te_junc_coord, canon_te_strand, side, half)
            if canon_fill:
                pre_list = list(pre_seq)
                if side == 'right':
                    for j in range(half):
                        if pre_list[j] == 'N' and canon_fill[j] in 'ACGT':
                            pre_list[j] = canon_fill[j]
                            n_canon_fill += 1
                else:
                    for j in range(half):
                        if pre_list[half + j] == 'N' and canon_fill[j] in 'ACGT':
                            pre_list[half + j] = canon_fill[j]
                            n_canon_fill += 1
                pre_seq = ''.join(pre_list)

            abs_gstart = region_start + abs_start        # 1-based genomic
            abs_gend   = abs_gstart + len(abs_seq) - 1

            out_name = os.path.join(
                outdir, f"junction_{side}_{cluster_id:03d}.fasta")

            iupac_count = sum(1 for c in pre_seq if c not in 'ACGTN')
            # TE coverage after canonical fill
            te_half  = pre_seq[half:] if side == 'left' else pre_seq[:half]
            te_cov   = sum(1 for c in te_half if c != 'N')

            canon_str = f"  can={n_canon_fill}" if n_canon_fill else ''
            tag = (f"ins={chrom}:{abs_ins_pos}  te={te_name}  reads={len(anchored)}"
                   f"  SNPs={iupac_count}  TE_cov={read_te_cov}/{half}"
                   f"  filled={te_cov}/{half}{canon_str}")

            if iupac_count > max_snps:
                print(f"  SKIP {os.path.basename(out_name)}  {tag}"
                      f"  ← SNPs>{max_snps}, likely merged clusters")
                continue

            # Case convention: reference bases lowercase, TE bases UPPERCASE.
            # The case boundary at position 50 marks the insertion point.
            abs_out = abs_seq.lower()
            if side == 'left':
                # left half = ref flank (lower), right half = TE start (UPPER)
                pre_out = pre_seq[:half].lower() + pre_seq[half:].upper()
            else:
                # left half = TE end (UPPER), right half = ref flank (lower)
                pre_out = pre_seq[:half].upper() + pre_seq[half:].lower()

            with open(out_name, 'w') as fh:
                fh.write(f">WT_REF[{abs_gstart}:{abs_gend}] "
                         f"insertion={chrom}:{abs_ins_pos} "
                         f"te={te_name} type={side}\n")
                fh.write(abs_out + '\n')
                fh.write(">REF\n" + abs_out + '\n')
                fh.write(f">reads_consensus_{side}_{cluster_id}\n")
                fh.write(pre_out + '\n')
                fh.write(f">{te_name}\n" + te_sub.upper() + '\n')

            print(f"  OK   {os.path.basename(out_name)}  {tag}")
            n_written += 1

    print(f"\nPhase 1 (read-level discovery): {n_written} junction files written")


if __name__ == '__main__':
    main()

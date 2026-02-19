#!/usr/bin/env python3
"""
find_te_junctions.py

SPAdes graph + k-mer walk TE junction discovery. Replaces the read-level
BLAST approach (build_junctions_from_reads.py).

Strategy:
  1. Parse SPAdes assembly graph (FASTG) into node sequences
  2. BLAST graph nodes vs TE database and reference region
  3. Classify nodes: junction (both hits), TE-boundary, TE-interior, ref-only
  4. Build a read k-mer index from R1.fq / R2.fq
  5. K-mer walk from TE-boundary nodes into reference sequence
  6. Cross-validate discovered junctions against gold standard catalog
  7. Write junction files (4-record format, compatible with downstream)

Usage:
    python find_te_junctions.py \\
        --outdir  <outdir>              \\
        --te-fasta <te.fasta>           \\
        --region  <chr:start-end>       \\
        [--gold-standard <gold_standard.tsv>] \\
        [--walk-len 200]                \\
        [--kmer-size 31]
"""

import argparse
import os
import subprocess
import sys
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq


# ---------------------------------------------------------------------------
# IUPAC encoding (same as previous pipeline for downstream compatibility)
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


def bases_to_iupac(counts, min_freq=0.15, min_count=2):
    """Return IUPAC code for a base-count column."""
    total = sum(counts.values())
    if total == 0:
        return 'N'
    major = max(counts, key=counts.get)
    present = {major}
    for b, n in counts.items():
        if b != major and n >= min_count and n / total >= min_freq:
            present.add(b)
    return IUPAC_FROM_BASES.get(frozenset(present), major)


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------
def revcomp(seq):
    return str(Seq(seq).reverse_complement())


def read_fastq_seqs(path):
    """Yield uppercase sequences from a FASTQ file."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return
    with open(path) as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().rstrip()
            f.readline()  # +
            f.readline()  # qual
            yield seq.upper()


# ---------------------------------------------------------------------------
# Step 4a: Parse FASTG
# ---------------------------------------------------------------------------
def parse_fastg(path):
    """Extract node sequences from SPAdes FASTG file.

    FASTG header format:
      >EDGE_1_length_1000_cov_50.5:EDGE_2',...;
    We strip edge info after ':' and the trailing ';'.
    """
    nodes = {}
    current_name = None
    current_seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if current_name:
                    nodes[current_name] = ''.join(current_seq)
                # Strip edge info: ">EDGE_1_length_1000_cov_50.5:EDGE_2,...;"
                name = line[1:].split(':')[0].rstrip(';').strip()
                current_name = name
                current_seq = []
            else:
                current_seq.append(line.strip())
    if current_name:
        nodes[current_name] = ''.join(current_seq)
    return nodes


def write_nodes_fasta(nodes, out_path):
    """Write graph nodes as FASTA for BLAST."""
    with open(out_path, 'w') as f:
        for name, seq in nodes.items():
            f.write(f">{name}\n{seq}\n")


# ---------------------------------------------------------------------------
# Step 4b: BLAST graph nodes
# ---------------------------------------------------------------------------
BLAST_FMT = "6 qseqid sseqid qstart qend sstart send sstrand pident length"


def run_blastn(query, subject, out_tsv, evalue="1e-5"):
    cmd = ["blastn", "-query", query, "-subject", subject,
           "-outfmt", BLAST_FMT, "-evalue", evalue, "-out", out_tsv]
    subprocess.run(cmd, check=True, stderr=subprocess.DEVNULL)


def parse_blast_hits(tsv_path):
    """Return list of hit dicts from BLAST output."""
    hits = []
    if not os.path.exists(tsv_path):
        return hits
    with open(tsv_path) as f:
        for line in f:
            parts = line.rstrip().split('\t')
            if len(parts) < 9:
                continue
            hits.append({
                'qseqid':  parts[0],
                'sseqid':  parts[1],
                'qstart':  int(parts[2]),
                'qend':    int(parts[3]),
                'sstart':  int(parts[4]),
                'send':    int(parts[5]),
                'sstrand': parts[6],
                'pident':  float(parts[7]),
                'length':  int(parts[8]),
            })
    return hits


def get_best_hits_per_query(hits):
    """Group by qseqid, return {qseqid: [hits sorted by length desc]}."""
    by_query = defaultdict(list)
    for h in hits:
        by_query[h['qseqid']].append(h)
    for qid in by_query:
        by_query[qid].sort(key=lambda h: h['length'], reverse=True)
    return dict(by_query)


def classify_nodes(te_hits_by_node, ref_hits_by_node, te_lengths,
                    min_ref_len=50):
    """Classify each node as junction, te_boundary, te_interior, or ref_only.

    A node needs both a TE hit and a substantial ref hit (>=min_ref_len bp)
    to be classified as junction.  A weak/spurious ref hit is ignored so
    TE-boundary nodes aren't misclassified.

    Returns list of dicts with classification and relevant hit info.
    """
    all_nodes = set(te_hits_by_node) | set(ref_hits_by_node)
    classified = []

    for node_id in all_nodes:
        te_hits = te_hits_by_node.get(node_id, [])
        ref_hits = ref_hits_by_node.get(node_id, [])

        # Only count substantial ref hits
        good_ref = [h for h in ref_hits if h['length'] >= min_ref_len]

        if te_hits and good_ref:
            # Junction node: both TE and substantial reference
            classified.append({
                'node_id': node_id,
                'class': 'junction',
                'te_hit': te_hits[0],
                'ref_hit': good_ref[0],
            })
        elif te_hits:
            best = te_hits[0]
            te_name = best['sseqid']
            te_len = te_lengths.get(te_name, 999999)
            sstart = best['sstart']
            send = best['send']
            s_lo = min(sstart, send)
            s_hi = max(sstart, send)

            # Near a TE terminus? (within 50bp of start or end)
            near_start = s_lo <= 50
            near_end = s_hi >= te_len - 50

            if near_start or near_end:
                classified.append({
                    'node_id': node_id,
                    'class': 'te_boundary',
                    'te_hit': best,
                    'near_start': near_start,
                    'near_end': near_end,
                })
            else:
                classified.append({
                    'node_id': node_id,
                    'class': 'te_interior',
                    'te_hit': best,
                })
        elif ref_hits:
            classified.append({
                'node_id': node_id,
                'class': 'ref_only',
                'ref_hit': ref_hits[0],
            })

    return classified


# ---------------------------------------------------------------------------
# Step 4c: Build read k-mer index
# ---------------------------------------------------------------------------
def build_kmer_index(r1_path, r2_path, k=31):
    """Build index: kmer → {A: count, C: count, G: count, T: count}.

    Index both forward and RC of every read so walks work in either direction.
    """
    index = defaultdict(lambda: defaultdict(int))
    n_reads = 0
    for path in [r1_path, r2_path]:
        for seq in read_fastq_seqs(path):
            n_reads += 1
            for s in [seq, revcomp(seq)]:
                for i in range(len(s) - k):
                    kmer = s[i:i+k]
                    next_base = s[i+k]
                    if 'N' not in kmer and next_base in 'ACGT':
                        index[kmer][next_base] += 1
    return dict(index), n_reads


# ---------------------------------------------------------------------------
# Step 4d: K-mer walk
# ---------------------------------------------------------------------------
def kmer_walk_right(seed_seq, kmer_index, max_extend=200, k=31,
                    min_vote_frac=0.6, iupac_threshold=0.15):
    """Extend seed_seq rightward using read k-mer evidence.

    Returns (extended_sequence, n_iupac_positions).

    The walk always appends the majority base (so the k-mer query window
    stays on a concrete sequence), but counts positions where a minor
    allele exceeds the IUPAC threshold.  If the majority base drops
    below min_vote_frac AND there's no clear minor allele, the walk stops.
    """
    extended = list(seed_seq)
    n_iupac = 0
    for _ in range(max_extend):
        query = ''.join(extended[-k:])
        counts = kmer_index.get(query)
        if not counts:
            break
        total = sum(counts.values())
        if total == 0:
            break
        best = max(counts, key=counts.get)
        best_frac = counts[best] / total

        # Check for minor allele worthy of IUPAC encoding
        has_minor = False
        for b, n in counts.items():
            if b != best and n >= 2 and n / total >= iupac_threshold:
                has_minor = True
                break

        if best_frac < min_vote_frac and not has_minor:
            break  # truly ambiguous — stop walking

        if has_minor:
            n_iupac += 1

        # Always extend with majority base to keep k-mer query concrete
        extended.append(best)
    return ''.join(extended), n_iupac


def kmer_walk_left(seed_seq, kmer_index, max_extend=200, k=31,
                   min_vote_frac=0.6):
    """Extend leftward by walking right on the reverse complement."""
    rc_extended, n_iupac = kmer_walk_right(
        revcomp(seed_seq), kmer_index, max_extend, k, min_vote_frac)
    return revcomp(rc_extended), n_iupac


# ---------------------------------------------------------------------------
# Step 4e: Extract 100bp junction from walked sequence
# ---------------------------------------------------------------------------
def extract_junction_100bp(walked_seq, te_boundary_pos, side, half=50):
    """Extract 100bp centered on the junction point.

    te_boundary_pos: 0-based position in walked_seq where TE ends / ref begins.
    side: 'right' (TE on left, ref on right) or 'left' (ref on left, TE on right)

    For a RIGHT junction (TE end → reference):
      walked_seq = [...TE_sequence...][...reference_extension...]
      te_boundary_pos = end of TE portion
      Take half bp of TE (before boundary) + half bp of ref (after boundary)

    For a LEFT junction (reference → TE start):
      walked_seq = [...reference_extension...][...TE_sequence...]
      te_boundary_pos = start of TE portion
      Take half bp of ref (before boundary) + half bp of TE (after boundary)
    """
    if side == 'right':
        # TE is on the left, reference on the right
        te_start = te_boundary_pos - half
        ref_end = te_boundary_pos + half
        if te_start < 0 or ref_end > len(walked_seq):
            return None
        return walked_seq[te_start:ref_end]
    else:
        # Reference on left, TE on right
        ref_start = te_boundary_pos - half
        te_end = te_boundary_pos + half
        if ref_start < 0 or te_end > len(walked_seq):
            return None
        return walked_seq[ref_start:te_end]


# ---------------------------------------------------------------------------
# Step 4f: Cross-validate against gold standard
# ---------------------------------------------------------------------------
def load_gold_standard(path):
    """Load gold_standard.tsv: te_name, approx_pos, total_reads, then per-sample columns."""
    entries = []
    if not path or not os.path.exists(path):
        return entries
    with open(path) as f:
        header = f.readline().rstrip().split('\t')
        for line in f:
            parts = line.rstrip().split('\t')
            if len(parts) < 3:
                continue
            entries.append({
                'te_name': parts[0],
                'approx_pos': int(parts[1]),
                'total_reads': int(parts[2]),
                'samples': dict(zip(header[3:], [int(x) for x in parts[3:]])),
            })
    return entries


def validate_junction(te_name, genomic_pos, gold_entries, tolerance=1000):
    """Check if a discovered junction matches a gold standard entry.

    Returns (status, matching_entry_or_None).
    """
    for entry in gold_entries:
        if entry['te_name'] == te_name:
            if abs(genomic_pos - entry['approx_pos']) <= tolerance:
                return 'VALIDATED', entry
    # Check if TE name matches but position is off
    for entry in gold_entries:
        if entry['te_name'] == te_name:
            return 'POSITION_MISMATCH', entry
    return 'SUSPECT', None


# ---------------------------------------------------------------------------
# BLAST-based genomic coordinate lookup for walked sequence
# ---------------------------------------------------------------------------
def blast_seq_vs_ref(seq, region_fa, outdir, tag):
    """BLAST a sequence against the reference region and return best hit."""
    query_path = os.path.join(outdir, f'_walk_{tag}.fa')
    out_path = os.path.join(outdir, f'_walk_{tag}_vs_ref.tsv')
    with open(query_path, 'w') as f:
        f.write(f">walk_{tag}\n{seq}\n")
    run_blastn(query_path, region_fa, out_path)
    hits = parse_blast_hits(out_path)
    if hits:
        hits.sort(key=lambda h: h['length'], reverse=True)
        return hits[0]
    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--outdir',         required=True)
    ap.add_argument('--te-fasta',       required=True)
    ap.add_argument('--region',         required=True,
                    help='e.g. chr3L:8710861-8744900')
    ap.add_argument('--gold-standard',  default=None,
                    help='Path to gold_standard.tsv')
    ap.add_argument('--walk-len',       type=int, default=200,
                    help='Max bases to extend via k-mer walk (default 200)')
    ap.add_argument('--kmer-size',      type=int, default=31,
                    help='K-mer size for walk index (default 31)')
    ap.add_argument('--min-vote-frac',  type=float, default=0.6,
                    help='Min fraction for majority base during walk (default 0.6)')
    ap.add_argument('--snp-min-freq',   type=float, default=0.15,
                    help='Minor-allele frequency for IUPAC encoding (default 0.15)')
    args = ap.parse_args()

    outdir   = args.outdir
    te_fasta = args.te_fasta
    region   = args.region
    k        = args.kmer_size
    walk_len = args.walk_len

    chrom        = region.split(':')[0]
    region_start = int(region.split(':')[1].split('-')[0])

    region_fa = os.path.join(outdir, 'region.fasta')
    r1_path   = os.path.join(outdir, 'R1.fq')
    r2_path   = os.path.join(outdir, 'R2.fq')

    # Locate assembly graph
    graph_path = os.path.join(outdir, 'assembly', 'assembly_graph.fastg')
    if not os.path.exists(graph_path):
        # Try per-k-value graphs
        for kval in [55, 33, 21]:
            alt = os.path.join(outdir, 'assembly', f'K{kval}', 'assembly_graph.fastg')
            if os.path.exists(alt):
                graph_path = alt
                print(f"Using per-k graph: {alt}")
                break
        else:
            print(f"ERROR: No assembly graph found in {outdir}/assembly/")
            print("  Tried: assembly_graph.fastg, K55/, K33/, K21/")
            sys.exit(1)

    # ------------------------------------------------------------------
    # Step 4a: Parse graph nodes
    # ------------------------------------------------------------------
    print("Step 4a: Parsing assembly graph")
    nodes = parse_fastg(graph_path)
    print(f"  {len(nodes)} graph nodes extracted")

    if not nodes:
        print("  ERROR: No nodes in assembly graph")
        sys.exit(1)

    nodes_fa = os.path.join(outdir, 'graph_nodes.fasta')
    write_nodes_fasta(nodes, nodes_fa)

    # ------------------------------------------------------------------
    # Step 4b: BLAST graph nodes vs TE database and reference
    # ------------------------------------------------------------------
    print("Step 4b: BLAST graph nodes vs TE database and reference")
    te_blast_out  = os.path.join(outdir, 'nodes_vs_te.tsv')
    ref_blast_out = os.path.join(outdir, 'nodes_vs_ref.tsv')

    run_blastn(nodes_fa, te_fasta, te_blast_out)
    run_blastn(nodes_fa, region_fa, ref_blast_out)

    te_hits  = parse_blast_hits(te_blast_out)
    ref_hits = parse_blast_hits(ref_blast_out)
    print(f"  TE hits: {len(te_hits)}, Ref hits: {len(ref_hits)}")

    te_hits_by_node  = get_best_hits_per_query(te_hits)
    ref_hits_by_node = get_best_hits_per_query(ref_hits)

    # Load TE lengths for boundary classification
    te_lengths = {}
    te_seq_dict = {}
    for rec in SeqIO.parse(te_fasta, 'fasta'):
        te_lengths[rec.id] = len(rec.seq)
        te_seq_dict[rec.id] = str(rec.seq).upper()

    # Load reference region
    region_seq = str(next(SeqIO.parse(region_fa, 'fasta')).seq).upper()

    # Classify nodes
    classified = classify_nodes(te_hits_by_node, ref_hits_by_node, te_lengths)
    by_class = defaultdict(list)
    for c in classified:
        by_class[c['class']].append(c)

    print(f"  Classification: "
          f"{len(by_class['junction'])} junction, "
          f"{len(by_class['te_boundary'])} TE-boundary, "
          f"{len(by_class['te_interior'])} TE-interior, "
          f"{len(by_class['ref_only'])} ref-only")

    te_nodes_total = (len(by_class['junction']) + len(by_class['te_boundary'])
                      + len(by_class['te_interior']))
    if te_nodes_total == 0:
        print("  No TE-containing nodes found — cannot discover junctions")
        return

    # ------------------------------------------------------------------
    # Step 4c: Build read k-mer index
    # ------------------------------------------------------------------
    print(f"Step 4c: Building read k-mer index (k={k})")
    kmer_index, n_reads = build_kmer_index(r1_path, r2_path, k)
    print(f"  Indexed {n_reads} reads, {len(kmer_index)} distinct k-mers")

    # ------------------------------------------------------------------
    # Step 4d-e: Process junction and TE-boundary nodes
    # ------------------------------------------------------------------
    print("Step 4d-e: Processing nodes → junction sequences")

    # Load gold standard for validation
    gold_entries = load_gold_standard(args.gold_standard)
    if gold_entries:
        print(f"  Gold standard: {len(gold_entries)} expected TEs")

    half = 50
    n_written = 0
    n_validated = 0
    n_suspect = 0
    missing_tes = set()  # Track gold standard TEs we haven't found
    found_tes = set()

    # Collect all junction candidates from both junction and TE-boundary nodes
    junction_candidates = []

    # --- Process JUNCTION nodes (already contain both TE and ref) ---
    for info in by_class['junction']:
        node_id = info['node_id']
        node_seq = nodes[node_id]
        te_hit = info['te_hit']
        ref_hit = info['ref_hit']
        te_name = te_hit['sseqid']

        # Determine which half of the node is TE vs ref (in node coordinates)
        te_mid = (te_hit['qstart'] + te_hit['qend']) / 2
        ref_mid = (ref_hit['qstart'] + ref_hit['qend']) / 2
        te_is_left_in_node = te_mid < ref_mid

        # Junction point: where TE and ref portions meet in the node
        if te_is_left_in_node:
            junc_pos = te_hit['qend']
            extract_side = 'right'  # TE left, ref right
        else:
            junc_pos = te_hit['qstart']
            extract_side = 'left'   # ref left, TE right

        # Extract 100bp around junction point (orientation not yet determined)
        junc_100 = extract_junction_100bp(node_seq, junc_pos, extract_side, half)
        if junc_100 is None:
            print(f"  {node_id}: junction too close to node boundary — skipped")
            continue

        # Use node-level ref BLAST for orientation (no re-BLAST needed).
        # If ref portion maps to minus strand, node is RC relative to genome.
        if ref_hit['sstrand'] == 'minus':
            junc_100 = revcomp(junc_100)
            te_is_left_in_node = not te_is_left_in_node

        # Determine L/R from which half is TE (in the now-oriented junction)
        # RIGHT: [TE end (0-49) | ref right flank (50-99)]
        # LEFT:  [ref left flank (0-49) | TE start (50-99)]
        if te_is_left_in_node:
            side = 'right'
            ins_region = min(ref_hit['sstart'], ref_hit['send'])
        else:
            side = 'left'
            ins_region = max(ref_hit['sstart'], ref_hit['send'])

        genomic_pos = region_start + ins_region - 1

        junction_candidates.append({
            'node_id': node_id,
            'source': 'junction_node',
            'side': side,
            'te_name': te_name,
            'te_hit': te_hit,
            'ins_region': ins_region,
            'genomic_pos': genomic_pos,
            'junction_seq': junc_100,
        })

    # --- Walk ALL TE-containing nodes (boundary + interior) ---
    # Try walking from BOTH ends of each node. The wrong direction will
    # either fail to extend or won't BLAST to the reference, so it
    # naturally filters. With only ~12 TE nodes the cost is trivial.
    te_walk_nodes = by_class['te_boundary'] + by_class['te_interior']
    print(f"\n  Attempting k-mer walks from {len(te_walk_nodes)} TE nodes "
          f"({len(by_class['te_boundary'])} boundary + "
          f"{len(by_class['te_interior'])} interior)")
    for info in te_walk_nodes:
        node_id = info['node_id']
        node_seq = nodes[node_id]
        te_hit = info['te_hit']
        te_name = te_hit['sseqid']

        for walk_dir in ['right', 'left']:
            if walk_dir == 'right':
                seed = node_seq[-min(len(node_seq), 200):]
                walked, n_iupac = kmer_walk_right(
                    seed, kmer_index, walk_len, k, args.min_vote_frac)
                extension_len = len(walked) - len(seed)

                if extension_len < 30:
                    continue

                # In the walked sequence: TE (seed) on left, extension on right
                junc_pos = len(seed)
                candidate_type = 'right'  # tentative: TE left, ref right

            else:
                seed = node_seq[:min(len(node_seq), 200)]
                walked, n_iupac = kmer_walk_left(
                    seed, kmer_index, walk_len, k, args.min_vote_frac)
                extension_len = len(walked) - len(seed)

                if extension_len < 30:
                    continue

                # In the walked sequence: extension on left, TE (seed) on right
                junc_pos = extension_len
                candidate_type = 'left'  # tentative: ref left, TE right

            junc_100 = extract_junction_100bp(walked, junc_pos, candidate_type, half)
            if junc_100 is None:
                continue

            # BLAST the ref half to verify it hits the reference
            if candidate_type == 'right':
                ref_half = junc_100[half:]
            else:
                ref_half = junc_100[:half]

            ref_hit_loc = blast_seq_vs_ref(ref_half, region_fa, outdir,
                                            f"walk_{node_id}_{walk_dir}")
            if ref_hit_loc is None:
                continue  # wrong direction — extension isn't reference

            # Fix orientation: if ref BLASTs minus, the walked sequence
            # is RC relative to the genome. RC and flip type.
            if ref_hit_loc['sstrand'] == 'minus':
                junc_100 = revcomp(junc_100)
                candidate_type = 'left' if candidate_type == 'right' else 'right'

            if candidate_type == 'right':
                ins_region = min(ref_hit_loc['sstart'], ref_hit_loc['send'])
            else:
                ins_region = max(ref_hit_loc['sstart'], ref_hit_loc['send'])
            genomic_pos = region_start + ins_region - 1

            print(f"  {node_id} walk_{walk_dir}: walked {extension_len}bp, "
                  f"junction {candidate_type} at {chrom}:{genomic_pos}")

            junction_candidates.append({
                'node_id': node_id,
                'source': 'kmer_walk',
                'side': candidate_type,
                'te_name': te_name,
                'te_hit': te_hit,
                'ins_region': ins_region,
                'genomic_pos': genomic_pos,
                'junction_seq': junc_100,
            })

    print(f"\n  Total junction candidates: {len(junction_candidates)}")

    # ------------------------------------------------------------------
    # Deduplicate: if multiple nodes produce the same TE + position (±20bp)
    # + side, keep the one from the better source (junction_node > kmer_walk)
    # ------------------------------------------------------------------
    junction_candidates.sort(key=lambda c: (c['te_name'], c['side'], c['genomic_pos']))
    deduped = []
    for cand in junction_candidates:
        merged = False
        for existing in deduped:
            if (existing['te_name'] == cand['te_name'] and
                existing['side'] == cand['side'] and
                abs(existing['genomic_pos'] - cand['genomic_pos']) <= 20):
                # Keep the one from a junction node if available
                if cand['source'] == 'junction_node' and existing['source'] != 'junction_node':
                    deduped.remove(existing)
                    deduped.append(cand)
                merged = True
                break
        if not merged:
            deduped.append(cand)

    print(f"  After dedup: {len(deduped)} unique junctions")

    # ------------------------------------------------------------------
    # Step 4f-g: Validate and write junction files
    # ------------------------------------------------------------------
    print("\nStep 4f-g: Validate and write junction files")

    for idx, cand in enumerate(deduped):
        te_name     = cand['te_name']
        side        = cand['side']
        genomic_pos = cand['genomic_pos']
        ins_region  = cand['ins_region']
        junc_seq    = cand['junction_seq']
        source      = cand['source']

        found_tes.add(te_name)

        # Cross-validate against gold standard
        if gold_entries:
            status, gs_entry = validate_junction(te_name, genomic_pos, gold_entries)
            gs_reads = gs_entry['total_reads'] if gs_entry else 0
            gs_tag = f"  gold={status} reads={gs_reads}"
            if status == 'VALIDATED':
                n_validated += 1
            elif status == 'SUSPECT':
                n_suspect += 1
        else:
            gs_tag = ""

        # Build Abs sequence (reference genome slice centred at insertion)
        ins_idx   = ins_region - 1  # 0-based in region_seq
        abs_start = ins_idx - half
        abs_end   = ins_idx + half
        if abs_start < 0 or abs_end > len(region_seq):
            print(f"  SKIP junction {idx}: too close to region boundary "
                  f"(ins_region={ins_region})")
            continue
        abs_seq = region_seq[abs_start:abs_end]

        # TE canonical sequence for record 4
        te_full = te_seq_dict.get(te_name, '')
        te_sub = te_full[:100] if side == 'left' else te_full[-100:]
        if len(te_sub) < 100:
            te_sub = te_sub.ljust(100, 'N')

        # Genomic coordinates for header
        abs_gstart = region_start + abs_start
        abs_gend   = abs_gstart + len(abs_seq) - 1

        # Case encoding: reference lowercase, TE UPPERCASE
        abs_out = abs_seq.lower()
        if side == 'left':
            pre_out = junc_seq[:half].lower() + junc_seq[half:].upper()
        else:
            pre_out = junc_seq[:half].upper() + junc_seq[half:].lower()

        # Count IUPAC codes in junction sequence
        iupac_count = sum(1 for c in junc_seq if c.upper() not in 'ACGTN')

        out_name = os.path.join(outdir, f"junction_{side}_{idx:03d}.fasta")

        tag = (f"ins={chrom}:{genomic_pos}  te={te_name}  "
               f"source={source}  SNPs={iupac_count}{gs_tag}")

        with open(out_name, 'w') as fh:
            fh.write(f">WT_REF[{abs_gstart}:{abs_gend}] "
                     f"insertion={chrom}:{genomic_pos} "
                     f"te={te_name} type={side}\n")
            fh.write(abs_out + '\n')
            fh.write(">REF\n" + abs_out + '\n')
            fh.write(f">walked_junction_{side}_{idx}\n")
            fh.write(pre_out + '\n')
            fh.write(f">{te_name}\n" + te_sub.upper() + '\n')

        print(f"  OK   {os.path.basename(out_name)}  {tag}")
        n_written += 1

    # ------------------------------------------------------------------
    # Report missing gold standard TEs
    # ------------------------------------------------------------------
    if gold_entries:
        for entry in gold_entries:
            if entry['te_name'] not in found_tes:
                missing_tes.add(entry['te_name'])
                print(f"  MISSING  {entry['te_name']}  "
                      f"approx_pos={entry['approx_pos']}  "
                      f"gold_reads={entry['total_reads']}")

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    print(f"\nPhase 1 (graph + k-mer walk): {n_written} junction files written")
    if gold_entries:
        print(f"  Validated: {n_validated}  Suspect: {n_suspect}  "
              f"Missing: {len(missing_tes)}")


if __name__ == '__main__':
    main()

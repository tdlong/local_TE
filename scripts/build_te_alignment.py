#!/usr/bin/env python3
"""
build_te_alignment.py

Simple visualization: 50bp ref + 50bp TE, with junction spanning both.

Usage:
    python build_te_alignment.py reference.fa te.fa junctions.fa [region] [outdir]

    region: Optional. Format like "chr3L:8710861-8744900" to show real genomic coordinates.
    outdir: Optional. Directory for PAF files and output. Defaults to current directory.
"""

from Bio import SeqIO
from Bio.Seq import Seq
import sys
import os


def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


def parse_paf(paf_file):
    alignments = []
    if not os.path.exists(paf_file):
        return alignments
    with open(paf_file) as f:
        for line in f:
            if line.strip():
                fields = line.strip().split('\t')
                alignments.append({
                    'query': fields[0],
                    'qlen': int(fields[1]),
                    'qstart': int(fields[2]),
                    'qend': int(fields[3]),
                    'strand': fields[4],
                    'target': fields[5],
                    'tlen': int(fields[6]),
                    'tstart': int(fields[7]),
                    'tend': int(fields[8]),
                })
    return alignments


def build_junction_view(junc_name, junc_seq, ref_aln, te_aln, ref_seq, te_seqs, region_chrom=None, region_start=0):
    """
    Build a 100bp view: 50bp ref + 50bp TE, with junction spanning.
    """
    
    te_name = te_aln['target']
    if te_name not in te_seqs:
        print(f"  WARNING: TE {te_name} not found")
        return None
    te_seq = te_seqs[te_name]
    
    # RC junction if needed to match ref
    need_rc = (ref_aln['strand'] == '-')
    jlen = len(junc_seq)
    
    if need_rc:
        junc = reverse_complement(junc_seq)
        # Flip junction coordinates
        ref_qstart = jlen - ref_aln['qend']
        ref_qend = jlen - ref_aln['qstart']
        te_qstart = jlen - te_aln['qend']
        te_qend = jlen - te_aln['qstart']
        # Also RC the TE
        te = reverse_complement(te_seq)
        te_len = len(te_seq)
        te_tstart = te_len - te_aln['tend']
        te_tend = te_len - te_aln['tstart']
    else:
        junc = junc_seq
        ref_qstart = ref_aln['qstart']
        ref_qend = ref_aln['qend']
        te_qstart = te_aln['qstart']
        te_qend = te_aln['qend']
        te = te_seq
        te_tstart = te_aln['tstart']
        te_tend = te_aln['tend']
    
    ref_tstart = ref_aln['tstart']
    ref_tend = ref_aln['tend']
    
    # Determine insertion point in reference coordinates
    if ref_qend <= te_qstart:
        # Left junction: insertion is at ref_tend (where ref ends, TE begins)
        insertion_point = ref_tend
    else:
        # Right junction: insertion is at ref_tstart (where TE ends, ref begins)
        insertion_point = ref_tstart
    
    # Calculate real genomic coordinate if region provided
    if region_chrom:
        genomic_pos = region_start + insertion_point
        pos_str = f"{region_chrom}:{genomic_pos}"
    else:
        pos_str = f"position {insertion_point}"
    
    print(f"\n{'='*70}")
    print(f"INSERTION of {te_name} at {pos_str}")
    print(f"{'='*70}")
    print(f"  Junction: {junc_name}")
    print(f"  RC={need_rc}")
    print(f"  junc[{ref_qstart}:{ref_qend}] matches ref[{ref_tstart}:{ref_tend}]")
    print(f"  junc[{te_qstart}:{te_qend}] matches TE[{te_tstart}:{te_tend}]")
    
    # Find the transition point in the junction
    # This is where ref alignment ends and TE alignment begins (or vice versa)
    if ref_qend <= te_qstart:
        # Left junction: ref ends at ref_qend, TE starts at te_qstart
        transition = (ref_qend + te_qstart) // 2
        jtype = 'left'
    else:
        # Right junction: TE ends at te_qend, ref starts at ref_qstart
        transition = (te_qend + ref_qstart) // 2
        jtype = 'right'
    
    print(f"    Type: {jtype}, transition at junc position {transition}")
    
    # Extract 50bp on each side of transition from junction
    junc_start = max(0, transition - 50)
    junc_end = min(len(junc), transition + 50)
    junc_portion = junc[junc_start:junc_end]
    
    # Calculate padding needed
    left_pad = 0
    right_pad = 0
    if transition - 50 < 0:
        left_pad = 50 - transition
    if transition + 50 > len(junc):
        right_pad = transition + 50 - len(junc)
    
    junc_portion = '-' * left_pad + junc_portion + '-' * right_pad
    
    # For ref: extract portion that matches the junction portion we're showing
    # Junction's ref portion: junc[ref_qstart:ref_qend] <-> ref[ref_tstart:ref_tend]
    # We're showing junc[junc_start:junc_end], so find overlap with ref portion
    
    ref_junc_overlap_start = max(junc_start, ref_qstart)
    ref_junc_overlap_end = min(junc_end, ref_qend)
    
    # Calculate which portion of junction view overlaps with ref
    # junc[ref_qstart:ref_qend] <-> ref[ref_tstart:ref_tend]
    # Note: lengths may differ due to indels, so we anchor at boundaries
    view_ref_start = max(junc_start, ref_qstart)
    view_ref_end = min(junc_end, ref_qend)
    
    if view_ref_start < view_ref_end:
        # How many bases of junction's ref portion are visible
        junc_ref_visible = view_ref_end - view_ref_start
        
        if jtype == 'right':
            # Right junction: anchor at START of ref alignment (ref_qstart <-> ref_tstart)
            # Work forwards from ref_tstart
            ref_overlap_start = ref_tstart
            ref_overlap_end = min(ref_tend, ref_tstart + junc_ref_visible)
        else:
            # Left junction: anchor at END of ref alignment (ref_qend <-> ref_tend)
            # Work backwards from ref_tend
            ref_overlap_end = ref_tend
            ref_overlap_start = max(ref_tstart, ref_tend - junc_ref_visible)
        
        ref_portion = ref_seq[ref_overlap_start:ref_overlap_end]
        # Position in our 100bp view
        ref_pos_in_view = (view_ref_start - junc_start) + left_pad
    else:
        ref_portion = ""
        ref_pos_in_view = 50
    
    # For TE: extract portion that matches the junction portion we're showing
    # Junction's TE portion: junc[te_qstart:te_qend] <-> te[te_tstart:te_tend]
    
    te_junc_overlap_start = max(junc_start, te_qstart)
    te_junc_overlap_end = min(junc_end, te_qend)
    
    # Calculate which portion of junction view overlaps with TE
    # junc[te_qstart:te_qend] <-> te[te_tstart:te_tend]
    # Note: lengths may differ due to indels, so we anchor at boundaries
    view_te_start = max(junc_start, te_qstart)
    view_te_end = min(junc_end, te_qend)
    
    if view_te_start < view_te_end:
        # How many bases of junction's TE portion are visible
        junc_te_visible = view_te_end - view_te_start
        
        if jtype == 'right':
            # Right junction: anchor at END of alignment (te_qend <-> te_tend)
            # Work backwards from te_tend
            te_overlap_end = te_tend
            te_overlap_start = max(te_tstart, te_tend - junc_te_visible)
        else:
            # Left junction: anchor at START of alignment (te_qstart <-> te_tstart)
            # Work forwards from te_tstart
            te_overlap_start = te_tstart
            te_overlap_end = min(te_tend, te_tstart + junc_te_visible)
        
        te_portion = te[te_overlap_start:te_overlap_end]
        # Position in our 100bp view
        te_pos_in_view = (view_te_start - junc_start) + left_pad
    else:
        te_portion = ""
        te_pos_in_view = 0
    
    # Build the 100bp alignment
    # Reference and TE positioned where they overlap with junction
    
    ref_aligned = '-' * ref_pos_in_view + ref_portion + '-' * (100 - ref_pos_in_view - len(ref_portion))
    te_aligned = '-' * te_pos_in_view + te_portion + '-' * (100 - te_pos_in_view - len(te_portion))
    junc_aligned = junc_portion[:100]
    
    # Ensure all are exactly 100bp
    ref_aligned = ref_aligned[:100].ljust(100, '-')
    te_aligned = te_aligned[:100].ljust(100, '-')
    junc_aligned = junc_aligned[:100].ljust(100, '-')
    
    # WT_REF: What the reference looks like without the TE insertion
    # This is a continuous stretch of reference around the insertion point
    wt_start = max(0, insertion_point - 50)
    wt_end = min(len(ref_seq), insertion_point + 50)
    wt_ref = ref_seq[wt_start:wt_end]
    # Pad to 100bp
    left_wt_pad = 50 - (insertion_point - wt_start)
    wt_ref_aligned = '-' * left_wt_pad + wt_ref + '-' * (100 - left_wt_pad - len(wt_ref))
    
    print(f"\n    100bp view (insertion at position 50):")
    print(f"    WT_REF: {wt_ref_aligned}")
    print(f"    REF:    {ref_aligned}")
    print(f"    JUNC:   {junc_aligned}")
    print(f"    TE:     {te_aligned}")
    
    return {
        'name': junc_name,
        'type': jtype,
        'insertion_point': insertion_point,
        'te_name': te_name,
        'wt_ref_label': f"WT_REF[{wt_start}:{wt_end}] insertion={pos_str} te={te_name} type={jtype}",
        'wt_ref_seq': wt_ref_aligned,
        'ref_label': f"REF",
        'ref_seq': ref_aligned,
        'junc_label': junc_name + ("_RC" if need_rc else ""),
        'junc_seq': junc_aligned,
        'te_label': te_name + ("_RC" if need_rc else ""),
        'te_seq': te_aligned
    }


def parse_region(region_str):
    """Parse region string like 'chr3L:8710861-8744900' into (chrom, start, end)"""
    if not region_str:
        return None, 0, 0
    try:
        chrom, coords = region_str.split(':')
        start, end = coords.split('-')
        return chrom, int(start), int(end)
    except:
        return None, 0, 0


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)
    
    ref_file, te_file, junc_file = sys.argv[1:4]
    
    # Optional region argument for real genomic coordinates
    region_str = sys.argv[4] if len(sys.argv) > 4 else None
    region_chrom, region_start, region_end = parse_region(region_str)
    
    # Optional output directory (defaults to current directory)
    outdir = sys.argv[5] if len(sys.argv) > 5 else "."
    
    # PAF files are in the output directory
    junc_to_ref_paf = os.path.join(outdir, "junctions_to_ref.paf")
    junc_to_te_paf = os.path.join(outdir, "junctions_to_te.paf")
    
    for f in [ref_file, te_file, junc_file, junc_to_ref_paf, junc_to_te_paf]:
        if not os.path.exists(f):
            print(f"Error: {f} not found")
            sys.exit(1)
    
    print("="*60)
    print("TE Junction Alignment (100bp view)")
    print("="*60)
    
    ref_rec = SeqIO.read(ref_file, "fasta")
    ref_seq = str(ref_rec.seq)
    te_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(te_file, "fasta")}
    junc_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(junc_file, "fasta")}
    
    # Keep the best (longest aligned span) alignment per query
    def best_per_query(alignments):
        best = {}
        for a in alignments:
            q = a['query']
            span = a['qend'] - a['qstart']
            if q not in best or span > (best[q]['qend'] - best[q]['qstart']):
                best[q] = a
        return best
    
    ref_alns = best_per_query(parse_paf(junc_to_ref_paf))
    te_alns = best_per_query(parse_paf(junc_to_te_paf))
    
    print(f"\nJunction contigs: {len(junc_seqs)}")
    print(f"  with ref alignment: {len(ref_alns)}")
    print(f"  with TE alignment:  {len(te_alns)}")
    
    n_processed = 0
    for jname in junc_seqs:
        if jname not in ref_alns:
            print(f"\n  SKIPPED: {jname} -- no alignment to reference")
            print(f"    (contig length: {len(junc_seqs[jname])}bp, "
                  f"likely too little reference flank -- assembly issue)")
            if jname in te_alns:
                te_span = te_alns[jname]['qend'] - te_alns[jname]['qstart']
                ref_flank = len(junc_seqs[jname]) - te_span
                print(f"    (TE alignment to {te_alns[jname]['target']}: "
                      f"{te_span}bp TE + ~{ref_flank}bp flank)")
            continue
        if jname not in te_alns:
            print(f"\n  SKIPPED: {jname} -- no alignment to TE")
            print(f"    (contig length: {len(junc_seqs[jname])}bp)")
            if jname in ref_alns:
                print(f"    (has ref alignment: "
                      f"qpos {ref_alns[jname]['qstart']}-{ref_alns[jname]['qend']})")
            continue
        
        result = build_junction_view(
            jname, junc_seqs[jname],
            ref_alns[jname], te_alns[jname],
            ref_seq, te_seqs,
            region_chrom, region_start
        )
        if result:
            outfile = os.path.join(outdir, f"junction_{result['type']}_{jname[:20]}.fasta")
            with open(outfile, 'w') as f:
                f.write(f">{result['wt_ref_label']}\n{result['wt_ref_seq']}\n")
                f.write(f">{result['ref_label']}\n{result['ref_seq']}\n")
                f.write(f">{result['junc_label']}\n{result['junc_seq']}\n")
                f.write(f">{result['te_label']}\n{result['te_seq']}\n")
            print(f"  Written: {outfile}")
            n_processed += 1
    
    print(f"\nDone! {n_processed} junctions visualized out of {len(junc_seqs)} contigs.")


if __name__ == "__main__":
    main()

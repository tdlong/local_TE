#!/usr/bin/env bash
# Phase 1: TE junction assembly pipeline (pooled discovery)
# SPAdes graph + k-mer walk approach
set -e

#=============================================================================
# USAGE
#=============================================================================
if [ $# -lt 3 ]; then
    echo "Usage: $0 <region> <outdir> <bam1> [bam2 bam3 ...]"
    echo ""
    echo "  region: Genomic region, e.g., chr3L:8710861-8744900"
    echo "  outdir: Directory for intermediate and output files"
    echo "  bam:    One or more BAM files to pool reads from"
    echo ""
    echo "Environment variables (optional):"
    echo "  REF:     Reference genome FASTA (default: dm6_conTE.fa)"
    echo "  TEFASTA: TE sequence database FASTA"
    exit 1
fi

REGION="$1"
OUTDIR="$2"
shift 2
BAMS=("$@")

# Reference files (can be overridden via environment)
REF="${REF:-/dfs7/adl/sruckman/XQTL/XQTL2/ref/dm6_conTE.fa}"
TEFASTA="${TEFASTA:-/dfs7/adl/sruckman/XQTL/XQTL2/ref/transposon_sequence_set.fa}"

# Create output directory
mkdir -p "$OUTDIR"

# Redirect all output to log file
exec > "$OUTDIR/pipeline.log" 2>&1

# Parse region coordinates
REGION_CHR="${REGION%%:*}"
REGION_COORDS="${REGION#*:}"
REGION_START="${REGION_COORDS%-*}"
REGION_END="${REGION_COORDS#*-}"

# Extend region by 1kb on each side (we exclude reads near ends pointing out)
REGION_EXT="${REGION_CHR}:$((REGION_START - 1000))-$((REGION_END + 1000))"
END_MARGIN=700

echo "=== TE Junction Assembly Pipeline (SPAdes Graph + K-mer Walk) ==="
echo "REGION: $REGION (extended: $REGION_EXT)"
echo "BAMs (${#BAMS[@]}): ${BAMS[*]}"
echo "OUTDIR: $OUTDIR"

#=============================================================================
# SECTIONS 1+2: Extract candidate reads from BAMs + build gold standard
#   Skip if R1.fq and region.fasta already exist (re-extraction takes ~15 min).
#   Delete R1.fq to force a fresh extraction.
#=============================================================================
if [ -s "$OUTDIR/R1.fq" ] && [ -s "$OUTDIR/region.fasta" ]; then
    n_cached=$(grep -c '^@' "$OUTDIR/R1.fq" 2>/dev/null || echo 0)
    echo ""
    echo "=== Sections 1+2: Skipped (reads already extracted) ==="
    echo "  R1.fq: $n_cached reads  |  delete $OUTDIR/R1.fq to re-extract"
else

#-----------------------------------------------------------------------------
# SECTION 1: Extract candidate reads from ALL BAMs (with catalog)
#   Output candidates_catalog.tsv preserving sample identity and mate contig.
#   Criterion: maps to region (MAPQ>=20, primary) AND mate is NOT on the
#   same chromosome ($7 != "="). This captures in one pass:
#     - Gold: mate maps to FBte reference
#     - Junction: mate unmapped
#     - Discordant: mate maps to different chromosome (TE copy elsewhere)
#   Exclude: reads near region ends pointing outward (boundary, not TE)
#-----------------------------------------------------------------------------
echo ""
echo "=== Section 1: Extract candidate reads (pooled from ${#BAMS[@]} BAMs) ==="
module load samtools/1.15.1

> "$OUTDIR/candidates_catalog.tsv"

for bam in "${BAMS[@]}"; do
    bam_name=$(basename "$bam" .bam)
    echo "  Processing BAM: $bam_name"

    # Output: read_name  sample  mate_contig  read_pos  sam_flag
    samtools view -q 20 "$bam" "$REGION_EXT" | \
      awk -v sample="$bam_name" '
        int($2/256)%2==0 && int($2/2048)%2==0 {
          len = length($10)
          pos_end = $4 + len - 1
          reverse = int($2/16)%2
          if (reverse && pos_end >= '$REGION_START' && pos_end <= '$REGION_START' + '$END_MARGIN') next
          if (!reverse && $4 >= '$REGION_END' - '$END_MARGIN' && $4 <= '$REGION_END') next
          if ($7 != "=") print $1 "\t" sample "\t" $7 "\t" $4 "\t" $2
        }' >> "$OUTDIR/candidates_catalog.tsv"

    echo "    done"
done

# Derive deduplicated read name list from catalog
awk '{print $1}' "$OUTDIR/candidates_catalog.tsv" | sort -u > "$OUTDIR/candidates.txt"
echo "Total candidates: $(wc -l < "$OUTDIR/candidates.txt")"
echo "Catalog entries:  $(wc -l < "$OUTDIR/candidates_catalog.tsv")"

#-----------------------------------------------------------------------------
# SECTION 2: Extract paired-end reads + build gold standard catalog
#   Extract all alignments for candidate read names, name-sort, then
#   convert to paired-end FASTQ.
#   Then build gold standard from FBte* mate-mapped reads.
#-----------------------------------------------------------------------------
echo ""
echo "=== Section 2: Extract paired-end reads (pooled) ==="

# Extract primary alignments for candidate reads from each BAM, merge into one name-sorted BAM
> "$OUTDIR/cand_bam_list.txt"
for bam in "${BAMS[@]}"; do
    bam_name=$(basename "$bam" .bam)
    samtools view -b -F 0x900 "$bam" -N "$OUTDIR/candidates.txt" \
        > "$OUTDIR/cand_${bam_name}.bam"
    echo "$OUTDIR/cand_${bam_name}.bam" >> "$OUTDIR/cand_bam_list.txt"
    echo "  $bam_name: extracted"
done

# Merge (handles single or multiple BAMs) and name-sort
if [ ${#BAMS[@]} -eq 1 ]; then
    samtools sort -n -o "$OUTDIR/candidates_namesorted.bam" \
        "$OUTDIR/cand_$(basename "${BAMS[0]}" .bam).bam"
else
    samtools merge -o "$OUTDIR/candidates_merged.bam" -b "$OUTDIR/cand_bam_list.txt"
    samtools sort -n -o "$OUTDIR/candidates_namesorted.bam" "$OUTDIR/candidates_merged.bam"
    rm -f "$OUTDIR/candidates_merged.bam"
fi

# Convert to paired-end FASTQ
samtools fastq \
    -1 "$OUTDIR/R1.fq" \
    -2 "$OUTDIR/R2.fq" \
    -s "$OUTDIR/singles.fq" \
    -0 "$OUTDIR/other.fq" \
    "$OUTDIR/candidates_namesorted.bam"

# Count reads
n_r1=$(grep -c '^@' "$OUTDIR/R1.fq" 2>/dev/null || echo 0)
n_r2=$(grep -c '^@' "$OUTDIR/R2.fq" 2>/dev/null || echo 0)
n_singles=$(grep -c '^@' "$OUTDIR/singles.fq" 2>/dev/null || echo 0)
echo "Paired: $n_r1 R1 + $n_r2 R2, Singles: $n_singles"

# Get reference region
samtools faidx "$REF" "$REGION" > "$OUTDIR/region.fasta"

# Clean up temp BAMs
rm -f "$OUTDIR"/cand_*.bam "$OUTDIR/cand_bam_list.txt" "$OUTDIR/candidates_namesorted.bam"

#-----------------------------------------------------------------------------
# Build gold standard catalog from FBte* mate-mapped reads
#   Filter catalog for mates mapping to FBte* contigs, then summarize
#   per TE name: approximate position (median), total reads, per-sample counts.
#-----------------------------------------------------------------------------
echo ""
echo "=== Section 2b: Build gold standard catalog ==="

# Filter for FBte* mates
awk '$3 ~ /^FBte/' "$OUTDIR/candidates_catalog.tsv" > "$OUTDIR/gold_standard_reads.tsv"
n_gold=$(wc -l < "$OUTDIR/gold_standard_reads.tsv")
echo "Gold standard reads (FBte* mates): $n_gold"

if [ "$n_gold" -gt 0 ]; then
    # Build summary: te_name  approx_pos  total_reads  sample1_count  sample2_count ...
    # Get sample names in order
    sample_names=$(awk '{print $2}' "$OUTDIR/candidates_catalog.tsv" | sort -u | tr '\n' '\t' | sed 's/\t$//')

    python3 -c "
import sys
from collections import defaultdict

reads_file = sys.argv[1]
out_file = sys.argv[2]

# Parse gold standard reads: read_name  sample  mate_contig(=te_name)  read_pos  flag
te_data = defaultdict(lambda: {'positions': [], 'samples': defaultdict(int)})
all_samples = set()

with open(reads_file) as f:
    for line in f:
        parts = line.rstrip().split('\t')
        if len(parts) < 5:
            continue
        read_name, sample, te_name, pos, flag = parts[:5]
        te_data[te_name]['positions'].append(int(pos))
        te_data[te_name]['samples'][sample] += 1
        all_samples.add(sample)

samples = sorted(all_samples)

with open(out_file, 'w') as f:
    f.write('te_name\tapprox_pos\ttotal_reads\t' + '\t'.join(samples) + '\n')
    for te_name in sorted(te_data.keys()):
        d = te_data[te_name]
        positions = sorted(d['positions'])
        median_pos = positions[len(positions) // 2]
        total = len(positions)
        sample_counts = '\t'.join(str(d['samples'].get(s, 0)) for s in samples)
        f.write(f'{te_name}\t{median_pos}\t{total}\t{sample_counts}\n')
" "$OUTDIR/gold_standard_reads.tsv" "$OUTDIR/gold_standard.tsv"

    echo "Gold standard catalog:"
    cat "$OUTDIR/gold_standard.tsv"
else
    echo "  No FBte* mates found — no gold standard catalog"
    > "$OUTDIR/gold_standard.tsv"
fi

fi  # end of skip block

#=============================================================================
# SECTION 3: SPAdes assembly
#   Use metagenome mode to preserve low-coverage TE junction paths.
#   We use the assembly graph (FASTG), not contigs.fasta.
#=============================================================================
echo ""
echo "=== Section 3: SPAdes assembly (--meta, paired-end) ==="

if [ -s "$OUTDIR/assembly/assembly_graph.fastg" ]; then
    n_nodes=$(grep -c '^>' "$OUTDIR/assembly/assembly_graph.fastg" 2>/dev/null || echo 0)
    echo "  Assembly graph already exists ($n_nodes nodes) — skipping SPAdes"
    echo "  Delete $OUTDIR/assembly/ to force reassembly"
else
    module load SPAdes/3.15.4

    rm -rf "$OUTDIR/assembly"
    echo "Running: spades.py -1 R1.fq -2 R2.fq -o assembly -k 21,33,55 --meta"
    spades.py \
        -1 "$OUTDIR/R1.fq" \
        -2 "$OUTDIR/R2.fq" \
        -o "$OUTDIR/assembly" \
        -k 21,33,55 --meta > "$OUTDIR/spades.log" 2>&1

    if [ -s "$OUTDIR/assembly/assembly_graph.fastg" ]; then
        n_nodes=$(grep -c '^>' "$OUTDIR/assembly/assembly_graph.fastg" 2>/dev/null || echo 0)
        echo "Assembly graph: $n_nodes nodes"
    else
        echo "WARNING: No assembly_graph.fastg produced. Check $OUTDIR/spades.log"
        echo "  Trying per-k-value graphs..."
        for kval in 55 33 21; do
            if [ -s "$OUTDIR/assembly/K${kval}/assembly_graph.fastg" ]; then
                echo "  Found K${kval}/assembly_graph.fastg"
            fi
        done
    fi

    if [ -s "$OUTDIR/assembly/contigs.fasta" ]; then
        echo "Contigs: $(grep -c '^>' "$OUTDIR/assembly/contigs.fasta") sequences (for reference)"
    fi
fi

#=============================================================================
# SECTION 4: Find TE junctions (SPAdes graph + k-mer walk)
#   Parse graph nodes, BLAST vs TE database, classify as junction/boundary,
#   k-mer walk from TE boundaries into reference, cross-validate against
#   gold standard catalog, write junction files.
#=============================================================================
echo ""
echo "=== Section 4: Find TE junctions (graph + k-mer walk) ==="
# Remove stale junction files so re-runs don't mix old and new results
rm -f "$OUTDIR"/junction_*.fasta
module load ncbi-blast/2.13.0
module load python/3.10.2

# Build gold standard arg if file exists and is non-empty
GS_ARG=""
if [ -s "$OUTDIR/gold_standard.tsv" ]; then
    GS_ARG="--gold-standard $OUTDIR/gold_standard.tsv"
fi

python scripts/find_te_junctions.py \
    --outdir   "$OUTDIR" \
    --te-fasta "$TEFASTA" \
    --region   "$REGION" \
    $GS_ARG

echo ""
echo "=== Phase 1 Complete ==="
echo "Output files in: $OUTDIR"
ls -la "$OUTDIR"/junction_*.fasta 2>/dev/null || echo "No junction files produced"

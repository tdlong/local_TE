#!/usr/bin/env bash
# Phase 1: TE junction assembly pipeline (pooled discovery)
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

echo "=== TE Junction Assembly Pipeline (Pooled Discovery) ==="
echo "REGION: $REGION (extended: $REGION_EXT)"
echo "BAMs (${#BAMS[@]}): ${BAMS[*]}"
echo "OUTDIR: $OUTDIR"

#=============================================================================
# SECTIONS 1+2: Extract candidate reads from BAMs
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
# SECTION 1: Extract candidate reads from ALL BAMs
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

> "$OUTDIR/candidates_all.txt"

for bam in "${BAMS[@]}"; do
    bam_name=$(basename "$bam" .bam)
    echo "  Processing BAM: $bam_name"

    # Reads in region, MAPQ>=20, primary, not near-end pointing out,
    # mate NOT on same chromosome
    samtools view -q 20 "$bam" "$REGION_EXT" | \
      awk '
        int($2/256)%2==0 && int($2/2048)%2==0 {
          len = length($10)
          pos_end = $4 + len - 1
          reverse = int($2/16)%2
          if (reverse && pos_end >= '$REGION_START' && pos_end <= '$REGION_START' + '$END_MARGIN') next
          if (!reverse && $4 >= '$REGION_END' - '$END_MARGIN' && $4 <= '$REGION_END') next
          if ($7 != "=") print $1
        }' | sort -u >> "$OUTDIR/candidates_all.txt"

    echo "    done"
done

# Deduplicate across BAMs
sort -u "$OUTDIR/candidates_all.txt" > "$OUTDIR/candidates.txt"
echo "Total candidates: $(wc -l < "$OUTDIR/candidates.txt")"

#-----------------------------------------------------------------------------
# SECTION 2: Extract paired-end reads
#   Extract all alignments for candidate read names, name-sort, then
#   convert to paired-end FASTQ.
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

# Get reference region (needed by build_junctions_from_reads.py)
samtools faidx "$REF" "$REGION" > "$OUTDIR/region.fasta"

# Clean up temp BAMs
rm -f "$OUTDIR"/cand_*.bam "$OUTDIR/cand_bam_list.txt" "$OUTDIR/candidates_namesorted.bam"

fi  # end of skip block

#=============================================================================
# SECTION 3: Read-level junction discovery
#   BLAST individual reads against TE DB and reference region to identify
#   junction reads (hits to both), cluster by inferred insertion position,
#   and build consensus junction sequences with IUPAC SNP encoding.
#   No assembly required — works at any coverage level.
#=============================================================================
echo ""
echo "=== Section 3: Read-level junction discovery ==="
# Remove stale junction files so re-runs don't mix old and new results
rm -f "$OUTDIR"/junction_*.fasta
module load ncbi-blast/2.13.0
module load python/3.10.2

python scripts/build_junctions_from_reads.py \
    --outdir   "$OUTDIR" \
    --te-fasta "$TEFASTA" \
    --region   "$REGION" \
    --min-support 3

echo ""
echo "=== Phase 1 Complete ==="
echo "Output files in: $OUTDIR"
ls -la "$OUTDIR"/junction_*.fasta 2>/dev/null || echo "No junction files produced"

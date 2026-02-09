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
# SECTION 1: Extract candidate reads from ALL BAMs
#   - Gold: maps uniquely to region AND has alignment to FBte (canonical TE)
#   - Junction: maps uniquely to region AND mate unmapped or maps elsewhere
#   - Exclude: reads near region ends pointing outward (they span region boundary, not TE)
#   - Pool candidates across all BAMs for maximum assembly coverage
#=============================================================================
echo ""
echo "=== Section 1: Extract candidate reads (pooled from ${#BAMS[@]} BAMs) ==="
module load samtools/1.15.1

# Filter: primary alignments, MAPQ>=20, exclude near-end pointing out
awk_filter='
  int($2/256)%2==0 && int($2/2048)%2==0 {
    len = length($10)
    pos_end = $4 + len - 1
    reverse = int($2/16)%2
    # Near left end pointing out (reverse strand, ends within margin of region start)
    if (reverse && pos_end >= '$REGION_START' && pos_end <= '$REGION_START' + '$END_MARGIN') next
    # Near right end pointing out (forward strand, starts within margin of region end)
    if (!reverse && $4 >= '$REGION_END' - '$END_MARGIN' && $4 <= '$REGION_END') next
    # If requiring mate unmapped, check FLAG bit 8
    if (require_mate_unmapped && int($2/8)%2 != 1) next
    print $1
  }'

# Initialize pooled files
> "$OUTDIR/gold_all.txt"
> "$OUTDIR/junction_all.txt"

for bam in "${BAMS[@]}"; do
    bam_name=$(basename "$bam" .bam)
    echo "  Processing BAM: $bam_name"

    # Get region-unique read names (MAPQ>=20, primary, not near-end-out)
    samtools view -q 20 "$bam" "$REGION_EXT" | \
      awk -v require_mate_unmapped=0 "$awk_filter" | sort -u > "$OUTDIR/region_unique_${bam_name}.txt"

    # Gold: region-unique AND has any alignment to FBte
    samtools view "$bam" -N "$OUTDIR/region_unique_${bam_name}.txt" | \
      awk '$3 ~ /^FBte/ || $7 ~ /^FBte/ {print $1}' | sort -u >> "$OUTDIR/gold_all.txt"

    # Junction: region-unique AND mate unmapped
    samtools view -q 20 "$bam" "$REGION_EXT" | \
      awk -v require_mate_unmapped=1 "$awk_filter" | sort -u >> "$OUTDIR/junction_all.txt"

    echo "    region_unique: $(wc -l < "$OUTDIR/region_unique_${bam_name}.txt")"
done

# Deduplicate across all BAMs
sort -u "$OUTDIR/gold_all.txt" > "$OUTDIR/gold.txt"
sort -u "$OUTDIR/junction_all.txt" > "$OUTDIR/junction.txt"
cat "$OUTDIR/gold.txt" "$OUTDIR/junction.txt" | sort -u > "$OUTDIR/candidates.txt"

echo "Gold (region+FBte): $(wc -l < "$OUTDIR/gold.txt")"
echo "Junction (mate unmapped): $(wc -l < "$OUTDIR/junction.txt")"
echo "Total candidates: $(wc -l < "$OUTDIR/candidates.txt")"

#=============================================================================
# SECTION 2: Build FASTA of candidate reads from ALL BAMs
#=============================================================================
echo ""
echo "=== Section 2: Build reads FASTA (pooled) ==="

> "$OUTDIR/reads.fasta"
for bam in "${BAMS[@]}"; do
    bam_name=$(basename "$bam" .bam)
    samtools view "$bam" -N "$OUTDIR/candidates.txt" | \
      awk 'int($2/256)%2==0 && int($2/2048)%2==0 {
        print ">" $1 "/" (int($2/64)%2==1 ? "1" : "2")
        print $10
      }' >> "$OUTDIR/reads.fasta"
    echo "  $bam_name: done"
done

echo "Reads FASTA: $OUTDIR/reads.fasta ($(grep -c '^>' "$OUTDIR/reads.fasta") sequences)"

# Also get reference region
samtools faidx "$REF" "$REGION" > "$OUTDIR/region.fasta"

#=============================================================================
# SECTION 3: Assemble with SPAdes
#=============================================================================
echo ""
echo "=== Section 3: SPAdes assembly ==="
module load SPAdes/3.15.4

# SPAdes options for low-depth TE junction assembly
# -k 21,33,55: smaller k-mers (skip 77 which often has no coverage)
# --only-assembler: skip error correction (can drop low-coverage TE reads)
SPADES_OPTS="-k 21,33,55 --only-assembler"

rm -rf "$OUTDIR/assembly"
echo "Running: spades.py -s $OUTDIR/reads.fasta -o $OUTDIR/assembly $SPADES_OPTS"
spades.py -s "$OUTDIR/reads.fasta" -o "$OUTDIR/assembly" $SPADES_OPTS > "$OUTDIR/spades.log" 2>&1

if [ ! -s "$OUTDIR/assembly/contigs.fasta" ]; then
  echo "ERROR: No contigs produced. Check $OUTDIR/spades.log" >&2
  exit 1
fi
echo "Contigs: $(grep -c '^>' "$OUTDIR/assembly/contigs.fasta") sequences"

#=============================================================================
# SECTION 4: BLAST (reality check)
#   - Reads should hit TE DB (they contain TE sequence)
#   - Contigs should also hit TE DB (assembly retained TE signal)
#   - If reads hit but contigs don't, assembly failed to capture junctions
#=============================================================================
echo ""
echo "=== Section 4: BLAST reality check ==="
module purge
module load ncbi-blast/2.13.0

makeblastdb -in "$TEFASTA" -dbtype nucl -out "$OUTDIR/te_db" -parse_seqids 2>/dev/null

blastn -query "$OUTDIR/reads.fasta" -db "$OUTDIR/te_db" -evalue 1e-5 \
  -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" \
  -out "$OUTDIR/reads_vs_te.tsv"
n_reads=$(wc -l < "$OUTDIR/reads_vs_te.tsv")
echo "Reads vs TE DB: $n_reads hits"

blastn -query "$OUTDIR/assembly/contigs.fasta" -db "$OUTDIR/te_db" -evalue 1e-5 \
  -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore" \
  -out "$OUTDIR/contigs_vs_te.tsv"
n_contigs=$(wc -l < "$OUTDIR/contigs_vs_te.tsv")
echo "Contigs vs TE DB: $n_contigs hits"

if [ "$n_reads" -gt 0 ] && [ "$n_contigs" -eq 0 ]; then
  echo "WARNING: Reads hit TEs but contigs do not. Assembly may have lost TE junctions."
fi

#=============================================================================
# SECTION 5: Collect sequences for alignment
#   - Region (reference)
#   - Contigs that hit TE (the junction-spanning fragments)
#   - The TE(s) that were hit
#=============================================================================
echo ""
echo "=== Section 5: Collect sequences for alignment ==="

# Contigs with good TE hits (>90% identity, >100bp)
awk '$3 > 90 && $4 > 100 {print $1}' "$OUTDIR/contigs_vs_te.tsv" | sort -u > "$OUTDIR/hit_contigs.txt"
awk '$3 > 90 && $4 > 100 {print $2}' "$OUTDIR/contigs_vs_te.tsv" | sort -u > "$OUTDIR/hit_tes.txt"

# Extract full FASTA records (multi-line aware)
awk 'NR==FNR {ids[$1]; next} 
     /^>/ {p=0; for(id in ids) if(index($0,">"id)>0){p=1;break}} 
     p' "$OUTDIR/hit_contigs.txt" "$OUTDIR/assembly/contigs.fasta" > "$OUTDIR/te_contigs.fasta" 2>/dev/null || true

awk 'NR==FNR {ids[$1]; next} 
     /^>/ {p=0; for(id in ids) if(index($0,">"id)>0){p=1;break}} 
     p' "$OUTDIR/hit_tes.txt" "$TEFASTA" > "$OUTDIR/te_seqs.fasta" 2>/dev/null || true

# Combine: region + contigs + TE
cat "$OUTDIR/region.fasta" "$OUTDIR/te_contigs.fasta" "$OUTDIR/te_seqs.fasta" > "$OUTDIR/for_alignment.fasta"

echo "Region: $OUTDIR/region.fasta"
echo "TE-hitting contigs: $OUTDIR/te_contigs.fasta ($(grep -c '^>' "$OUTDIR/te_contigs.fasta" 2>/dev/null || echo 0))"
echo "Matching TEs: $OUTDIR/te_seqs.fasta ($(grep -c '^>' "$OUTDIR/te_seqs.fasta" 2>/dev/null || echo 0))"
echo "Combined: $OUTDIR/for_alignment.fasta"

#=============================================================================
# SECTION 6: Generate alignments with minimap2
#=============================================================================
echo ""
echo "=== Section 6: Minimap2 alignment ==="
module load minimap2/2.28 2>/dev/null || true

# Map junction contigs to reference
minimap2 -c "$OUTDIR/region.fasta" "$OUTDIR/te_contigs.fasta" > "$OUTDIR/junctions_to_ref.paf" 2>/dev/null

# Map junction contigs to TE database
minimap2 -c "$OUTDIR/te_seqs.fasta" "$OUTDIR/te_contigs.fasta" > "$OUTDIR/junctions_to_te.paf" 2>/dev/null

echo "Generated: $OUTDIR/junctions_to_ref.paf, $OUTDIR/junctions_to_te.paf"

#=============================================================================
# SECTION 7: Build junction alignment visualization
#=============================================================================
echo ""
echo "=== Section 7: Junction visualization ==="
module load python/3.10.2

python scripts/build_te_alignment.py \
    "$OUTDIR/region.fasta" \
    "$OUTDIR/te_seqs.fasta" \
    "$OUTDIR/te_contigs.fasta" \
    "$REGION" \
    "$OUTDIR"

echo ""
echo "=== Phase 1 Complete ==="
echo "Output files in: $OUTDIR"
ls -la "$OUTDIR"/junction_*.fasta 2>/dev/null || echo "No junction files produced"

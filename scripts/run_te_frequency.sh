#!/usr/bin/env bash
# Phase 2: TE frequency estimation from raw FASTQs
set -e

#=============================================================================
# USAGE
#=============================================================================
if [ $# -lt 2 ]; then
    echo "Usage: $0 <junctions.fa> <samples.tsv> [outdir]"
    echo ""
    echo "  junctions.fa: Combined competitive-alignment reference from Phase 1"
    echo "  samples.tsv:  Tab-delimited file with columns:"
    echo "                  sample_name  /path/to/R1.fq.gz  /path/to/R2.fq.gz"
    echo "  outdir:       Working directory (default: freq_work)"
    echo ""
    echo "Output: <outdir>/frequency_results.tsv"
    exit 1
fi

JUNCTIONS_FA="$1"
SAMPLES_TSV="$2"
OUTDIR="${3:-freq_work}"

mkdir -p "$OUTDIR"

echo "=== Phase 2: TE Frequency Estimation ==="
echo "Junctions ref: $JUNCTIONS_FA"
echo "Samples file:  $SAMPLES_TSV"
echo "Working dir:   $OUTDIR"

#=============================================================================
# Index junctions.fa with bwa (once)
#=============================================================================
module load bwa/0.7.17

if [ ! -f "${JUNCTIONS_FA}.bwt" ]; then
    echo ""
    echo "=== Indexing junctions.fa with bwa ==="
    bwa index "$JUNCTIONS_FA"
fi

#=============================================================================
# Process each sample
#=============================================================================
module load python/3.10.2

# Write header to results file
echo -e "sample\tjunction_id\tabs_count\tpre_count\tfrequency" > "$OUTDIR/frequency_results.tsv"

n_samples=0
while IFS=$'\t' read -r sample_name r1 r2; do
    # Skip empty lines and comments
    [[ -z "$sample_name" || "$sample_name" == \#* ]] && continue

    n_samples=$((n_samples + 1))
    echo ""
    echo "=== [$n_samples] Processing: $sample_name ==="
    echo "  R1: $r1"
    echo "  R2: $r2"

    # Verify FASTQ files exist
    if [ ! -f "$r1" ]; then
        echo "  WARNING: R1 not found: $r1 -- skipping"
        continue
    fi
    if [ ! -f "$r2" ]; then
        echo "  WARNING: R2 not found: $r2 -- skipping"
        continue
    fi

    sam_file="$OUTDIR/${sample_name}.sam"

    # Align raw FASTQs to junction reference
    # -W 15: minimum seed length for sensitivity near junction edges
    # -a: output all alignments (let competitive scoring resolve ambiguity)
    echo "  Aligning with bwa mem..."
    bwa mem -t 8 -W 15 -a "$JUNCTIONS_FA" "$r1" "$r2" \
        > "$sam_file" 2>"$OUTDIR/${sample_name}_bwa.log"

    # Count spanning reads per junction
    echo "  Counting spanning reads..."
    python scripts/estimate_te_freq.py "$sam_file" "$sample_name" | \
        tail -n +2 >> "$OUTDIR/frequency_results.tsv"

    echo "  Done."

    # Optional: remove SAM to save disk (uncomment if space is tight)
    # rm -f "$sam_file"

done < "$SAMPLES_TSV"

echo ""
echo "=== Phase 2 Complete ==="
echo "Processed $n_samples samples"
echo "Results: $OUTDIR/frequency_results.tsv"
echo ""
echo "--- Results ---"
column -t "$OUTDIR/frequency_results.tsv"

#!/bin/bash
#SBATCH --job-name=te_discovery
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G

# Phase 1: TE Junction Discovery - SLURM wrapper
# Run from: /dfs7/adl/tdlong/Sarah/local_TE

source config.sh

echo "=== Phase 1: TE Junction Discovery (Pooled) ==="
echo "BAMs: ${BAMS[*]}"
echo "Regions: ${REGIONS[*]}"
echo ""

# Build full BAM paths
BAM_PATHS=()
for bam in "${BAMS[@]}"; do
    BAM_PATHS+=("$BAM_DIR/$bam")
done

# Process each region, pooling reads from ALL BAMs
for region in "${REGIONS[@]}"; do
    region_safe="${region//:/_}"
    outdir="temp_work/${region_safe}"

    echo "Processing region: $region"
    echo "  Pooling ${#BAMS[@]} BAMs"
    echo "  Output: $outdir"

    REF="$REF" TEFASTA="$TEFASTA" \
        bash scripts/run_te_assembly.sh "$region" "$outdir" "${BAM_PATHS[@]}"

    echo "  Log: $outdir/pipeline.log"
    echo ""
done

echo "=== Phase 1 Complete ==="
echo "Results in temp_work/"
ls -la temp_work/*/junction_*.fasta 2>/dev/null || echo "No junction files produced"
echo ""
echo "Next steps:"
echo "  1. python scripts/build_junctions_ref.py temp_work junctions.fa"
echo "  2. Fill in SAMPLES in config.sh with FASTQ paths"
echo "  3. sbatch scripts/submit_te_frequency.sh"

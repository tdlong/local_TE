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

echo "=== Building competitive reference from all junctions ==="
module load python/3.10.2
python scripts/build_junctions_ref.py temp_work junctions.fa

echo ""
echo "=== Phase 1 Complete ==="
echo ""
echo "Review junctions.fa before proceeding to Phase 2:"
echo "  cat junctions.fa"
echo ""
echo "When satisfied, run Phase 2:"
echo "  sbatch scripts/submit_te_kmer_count.sh"

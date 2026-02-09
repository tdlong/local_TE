#!/bin/bash
#SBATCH --job-name=te_assembly
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G

# TE Junction Assembly - SLURM wrapper
# Run from: /dfs7/adl/tdlong/Sarah/local_TE

# Configuration arrays (expand later for multiple BAMs/regions)
BAM_DIR="/dfs7/adl/sruckman/XQTL/XQTL2/data/bam/colorTEcon"
BAMS=("HOULE_L2F.bam")
REGIONS=("chr3L:8710861-8744900")

echo "=== TE Junction Assembly ==="
echo "BAMs: ${BAMS[*]}"
echo "Regions: ${REGIONS[*]}"
echo ""

# Iterate over all BAM x region combinations
for bam in "${BAMS[@]}"; do
    for region in "${REGIONS[@]}"; do
        # Create output directory name from bam and region
        bam_base="${bam%.bam}"
        region_safe="${region//:/_}"
        outdir="temp_work/${bam_base}_${region_safe}"
        
        echo "Processing: $bam x $region"
        echo "  Output: $outdir"
        
        mkdir -p "$outdir"
        
        bash scripts/run_te_assembly.sh "$region" "$BAM_DIR/$bam" "$outdir"
        
        echo "  Log: $outdir/pipeline.log"
        echo ""
    done
done

echo "=== Complete ==="
echo "Results in temp_work/"
ls -la temp_work/*/junction_*.fasta 2>/dev/null || echo "No junction files produced"

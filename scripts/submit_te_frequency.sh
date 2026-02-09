#!/bin/bash
#SBATCH --job-name=te_frequency
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G

# Phase 2: TE Frequency Estimation - SLURM wrapper
# Run from: /dfs7/adl/tdlong/Sarah/local_TE

source config.sh

# Step 1: Build combined junctions.fa from discovery results (if not done)
if [ ! -f "junctions.fa" ]; then
    echo "Building junctions.fa from discovery results..."
    module load python/3.10.2
    python scripts/build_junctions_ref.py temp_work junctions.fa
    echo ""
fi

# Step 2: Write samples.tsv from config
> samples.tsv
for entry in "${SAMPLES[@]}"; do
    echo "$entry" >> samples.tsv
done

echo "Wrote $(wc -l < samples.tsv) samples to samples.tsv"

# Step 3: Run frequency estimation
bash scripts/run_te_frequency.sh junctions.fa samples.tsv freq_work

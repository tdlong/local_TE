#!/bin/bash
#SBATCH --job-name=te_kmer_setup
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

# Phase 2: TE Genotyping via K-mer Counting - SLURM orchestrator
# Builds k-mer databases, submits per-sample BBDuk array job,
# then submits genotyping job with dependency.
#
# Run from: /dfs7/adl/tdlong/Sarah/local_TE
set -e

source config.sh

OUTDIR="kmer_work"
mkdir -p "$OUTDIR"

echo "=== Phase 2: TE Genotyping (K-mer Counting) ==="

# Step 1: Build junctions.fa from discovery results (if not done)
if [ ! -f "junctions.fa" ]; then
    echo "Building junctions.fa from discovery results..."
    module load python/3.10.2
    python scripts/build_junctions_ref.py temp_work junctions.fa
    echo ""
fi

# Step 2: Extract diagnostic k-mers
echo "Extracting diagnostic k-mers from junctions.fa..."
module load python/3.10.2
python scripts/extract_junction_kmers.py junctions.fa "$OUTDIR"
echo ""

# Step 3: Write samples.tsv from config
> samples.tsv
for entry in "${SAMPLES[@]}"; do
    echo "$entry" >> samples.tsv
done
N_SAMPLES=$(wc -l < samples.tsv)
echo "Wrote $N_SAMPLES samples to samples.tsv"

# Step 4: Submit BBDuk array job (one task per sample)
echo "Submitting BBDuk array job..."
ARRAY_JOB=$(sbatch --parsable \
    --array=1-${N_SAMPLES} \
    scripts/run_te_kmer_count.sh \
    samples.tsv "$OUTDIR" "$OUTDIR/abs_kmers.fa" "$OUTDIR/pre_kmers.fa")
echo "  Array job: $ARRAY_JOB"

# Step 5: Submit genotyping job after array completes
echo "Submitting genotyping job (depends on $ARRAY_JOB)..."
GENO_JOB=$(sbatch --parsable \
    --dependency=afterok:${ARRAY_JOB} \
    --job-name=te_genotype \
    -A tdlong_lab -p standard \
    --cpus-per-task=1 --mem-per-cpu=4G \
    --wrap="module load python/3.10.2; python scripts/genotype_from_counts.py $OUTDIR")
echo "  Genotype job: $GENO_JOB"

echo ""
echo "Pipeline submitted. Monitor with: squeue -u \$USER"
echo "Results will be in: $OUTDIR/genotype_results.tsv"

#!/bin/bash
#SBATCH --job-name=te_kmer
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G

# Per-sample BBDuk k-mer counting (designed for SLURM array job)
# Reads sample info from samples.tsv using SLURM_ARRAY_TASK_ID
set -e

SAMPLES_TSV="$1"
OUTDIR="$2"
ABS_KMERS="$3"
PRE_KMERS="$4"

if [ -z "$SAMPLES_TSV" ] || [ -z "$OUTDIR" ] || [ -z "$ABS_KMERS" ] || [ -z "$PRE_KMERS" ]; then
    echo "Usage: run_te_kmer_count.sh <samples.tsv> <outdir> <abs_kmers.fa> <pre_kmers.fa>"
    echo "  Normally called by submit_te_kmer_count.sh as a SLURM array job"
    exit 1
fi

# Read sample line for this array task
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_TSV")
SAMPLE=$(echo "$LINE" | cut -f1)
R1=$(echo "$LINE" | cut -f2)
R2=$(echo "$LINE" | cut -f3)

echo "=== BBDuk k-mer counting: $SAMPLE ==="
echo "  R1: $R1"
echo "  R2: $R2"

mkdir -p "$OUTDIR"

# Verify inputs
for f in "$R1" "$R2" "$ABS_KMERS" "$PRE_KMERS"; do
    if [ ! -f "$f" ]; then
        echo "Error: file not found: $f"
        exit 1
    fi
done

# BBDuk is Java-based; load Java before activating conda
module load java/17

# Activate conda/BBMap
source ~/miniforge3/etc/profile.d/conda.sh
source ~/miniforge3/etc/profile.d/mamba.sh
mamba activate bbmap

# Count absence allele k-mers
# -Xmx5g: explicit heap limit so BBDuk doesn't auto-detect node RAM
echo "  Counting absence k-mers..."
bbduk.sh -Xmx5g in="$R1" in2="$R2" \
    ref="$ABS_KMERS" \
    k=31 maskmiddle=f rcomp=t \
    stats="$OUTDIR/${SAMPLE}_abs_stats.txt"

# Count presence allele k-mers
echo "  Counting presence k-mers..."
bbduk.sh -Xmx5g in="$R1" in2="$R2" \
    ref="$PRE_KMERS" \
    k=31 maskmiddle=f rcomp=t \
    stats="$OUTDIR/${SAMPLE}_pre_stats.txt"

echo "  Done: $SAMPLE"

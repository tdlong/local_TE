#!/usr/bin/env bash
# test_phase1.sh -- Interactive end-to-end Phase 1 test using config.sh values.
#
# Sources config.sh, picks the first region and all configured BAMs, runs
# run_te_assembly.sh with output to stdout (not a log file) so you can watch
# progress in real time.
#
# Run from the project root on a compute node:
#   srun --pty bash
#   bash scripts/test_phase1.sh
#
# Optionally override the output directory:
#   TEST_OUTDIR=my_test_dir bash scripts/test_phase1.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

cd "$PROJECT_ROOT"

# Load config
source config.sh

# Pick first region
TEST_REGION="${REGIONS[0]}"
if [ -z "$TEST_REGION" ]; then
    echo "ERROR: No REGIONS defined in config.sh" >&2
    exit 1
fi

# Build full BAM paths
BAM_PATHS=()
for bam in "${BAMS[@]}"; do
    BAM_PATHS+=("$BAM_DIR/$bam")
done

if [ ${#BAM_PATHS[@]} -eq 0 ]; then
    echo "ERROR: No BAMS defined in config.sh" >&2
    exit 1
fi

# Output directory
region_safe="${TEST_REGION//:/_}"
TEST_OUTDIR="${TEST_OUTDIR:-temp_work/test_${region_safe}}"

echo "=== Phase 1 Interactive Test ==="
echo "Region:  $TEST_REGION"
echo "BAMs:    ${BAM_PATHS[*]}"
echo "Outdir:  $TEST_OUTDIR"
echo "REF:     $REF"
echo "TEFASTA: $TEFASTA"
if [ -s "$TEST_OUTDIR/R1.fq" ]; then
    echo ""
    echo "NOTE: $TEST_OUTDIR/R1.fq exists — BAM extraction will be skipped."
    echo "      Delete it to force re-extraction."
fi
echo ""

# pipeline.log is written inside OUTDIR; print it at the end.

REF="$REF" TEFASTA="$TEFASTA" \
    bash scripts/run_te_assembly.sh \
        "$TEST_REGION" \
        "$TEST_OUTDIR" \
        "${BAM_PATHS[@]}" || true

echo ""
echo "=== Pipeline log ==="
cat "$TEST_OUTDIR/pipeline.log" 2>/dev/null || echo "(no log file)"

echo ""
echo "=== Results ==="
echo "Junction files:"
ls -lh "$TEST_OUTDIR"/junction_*.fasta 2>/dev/null || echo "  (none)"

echo ""
echo "BLAST hit counts:"
echo "  reads_vs_te.tsv:  $(wc -l < "$TEST_OUTDIR/reads_vs_te.tsv" 2>/dev/null || echo 0) hits"
echo "  reads_vs_ref.tsv: $(wc -l < "$TEST_OUTDIR/reads_vs_ref.tsv" 2>/dev/null || echo 0) hits"

echo ""
echo "=== First junction file (if any) ==="
first=$(ls "$TEST_OUTDIR"/junction_*.fasta 2>/dev/null | head -1)
if [ -n "$first" ]; then
    cat "$first"
else
    echo "  (none produced)"
fi

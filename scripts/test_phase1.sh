#!/usr/bin/env bash
# test_phase1.sh -- Interactive end-to-end Phase 1 test using config.sh values.
#
# Sources config.sh, picks the first region and all configured BAMs, and runs
# run_te_assembly.sh. Prints a clean diagnostic summary on completion.
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

source config.sh

TEST_REGION="${REGIONS[0]}"
[ -n "$TEST_REGION" ] || { echo "ERROR: No REGIONS in config.sh" >&2; exit 1; }

BAM_PATHS=()
for bam in "${BAMS[@]}"; do BAM_PATHS+=("$BAM_DIR/$bam"); done
[ ${#BAM_PATHS[@]} -gt 0 ] || { echo "ERROR: No BAMS in config.sh" >&2; exit 1; }

region_safe="${TEST_REGION//:/_}"
TEST_OUTDIR="${TEST_OUTDIR:-temp_work/test_${region_safe}}"

echo "=== Phase 1 Interactive Test ==="
echo "Region:  $TEST_REGION"
echo "BAMs:    ${#BAM_PATHS[@]} files"
echo "Outdir:  $TEST_OUTDIR"
if [ -s "$TEST_OUTDIR/R1.fq" ]; then
    echo "Reads:   cached (delete R1.fq to re-extract)"
fi
echo ""

REF="$REF" TEFASTA="$TEFASTA" \
    bash scripts/run_te_assembly.sh \
        "$TEST_REGION" "$TEST_OUTDIR" "${BAM_PATHS[@]}" || true

LOG="$TEST_OUTDIR/pipeline.log"

# ------------------------------------------------------------------
# Summary: key metrics extracted from the pipeline log
# ------------------------------------------------------------------
echo ""
echo "=== Summary ==="

# Read/BLAST/cluster counts
grep -E "^  (R1\.fq|Paired:|Total candidates:|Found [0-9]+ junction|raw clusters)" \
    "$LOG" 2>/dev/null || true
grep -E "^Step [12]:" "$LOG" 2>/dev/null || true
grep -E "raw clusters" "$LOG" 2>/dev/null || true

echo ""
echo "Junction results (OK=written  SKIP=too many SNPs):"
echo ""
printf "  %-6s  %-32s  %-20s  %-14s  %5s  %5s  %8s  %8s  %4s\n" \
    STATUS FILE INS_POS TE READS SNPs TE_COV FILLED CAN
printf "  %-6s  %-32s  %-20s  %-14s  %5s  %5s  %8s  %8s  %4s\n" \
    "------" "------" "-------" "--" "-----" "----" "------" "------" "---"

grep -E "^  (OK|SKIP)" "$LOG" 2>/dev/null | \
while IFS= read -r line; do
    status=$(echo "$line"  | awk '{print $1}')
    file=$(echo "$line"    | awk '{print $2}')
    ins=$(echo "$line"     | grep -oP 'ins=\S+')
    te=$(echo "$line"      | grep -oP 'te=\S+')
    reads=$(echo "$line"   | grep -oP 'reads=\K[0-9]+')
    snps=$(echo "$line"    | grep -oP 'SNPs=\K[0-9]+')
    tecov=$(echo "$line"   | grep -oP 'TE_cov=\S+')
    filled=$(echo "$line"  | grep -oP 'filled=\S+')
    can=$(echo "$line"     | grep -oP 'can=\K[0-9]+')
    printf "  %-6s  %-32s  %-20s  %-14s  %5s  %5s  %8s  %8s  %4s\n" \
        "$status" "$file" "${ins#ins=}" "${te#te=}" "$reads" "$snps" \
        "${tecov#TE_cov=}" "${filled#filled=}" "${can:-0}"
done || echo "  (none)"

echo ""
n_ok=$(grep -c "^  OK" "$LOG" 2>/dev/null || echo 0)
n_skip=$(grep -c "^  SKIP" "$LOG" 2>/dev/null || echo 0)
echo "  Written: $n_ok    Skipped (SNPs>max): $n_skip"

# ------------------------------------------------------------------
# Show first written junction file as a spot-check
# ------------------------------------------------------------------
first=$(grep "^  OK" "$LOG" 2>/dev/null | head -1 | awk '{print $2}')
if [ -n "$first" ]; then
    echo ""
    echo "=== Spot-check: $first ==="
    cat "$TEST_OUTDIR/$first" 2>/dev/null || cat "$first" 2>/dev/null || true
fi

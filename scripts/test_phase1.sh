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

# Read extraction counts
grep -E "^  (R1\.fq|Paired:|Total candidates:|Catalog entries:)" \
    "$LOG" 2>/dev/null || true

# Gold standard catalog
echo ""
echo "Gold standard:"
grep -E "^Gold standard reads" "$LOG" 2>/dev/null || true
grep -A 100 "^Gold standard catalog:" "$LOG" 2>/dev/null | \
    head -20 || echo "  (none)"

# SPAdes assembly
echo ""
echo "Assembly:"
grep -E "^(Assembly graph:|Contigs:|  Assembly graph already|WARNING:)" \
    "$LOG" 2>/dev/null || true

# Graph node classification
grep -E "Classification:" "$LOG" 2>/dev/null || true

# K-mer walk stats
echo ""
echo "K-mer walk:"
grep -E "^  (Indexed|.*walked.*into reference)" "$LOG" 2>/dev/null || true
grep -E "Total junction candidates" "$LOG" 2>/dev/null || true
grep -E "After dedup" "$LOG" 2>/dev/null || true

echo ""
echo "Junction results (OK=written  SKIP=skipped  MISSING=expected but not found):"
echo ""
printf "  %-9s  %-32s  %-20s  %-14s  %s\n" \
    STATUS FILE INS_POS TE SOURCE
printf "  %-9s  %-32s  %-20s  %-14s  %s\n" \
    "--------" "------" "-------" "--" "------"

grep -E "^  (OK|SKIP|MISSING|VALIDATED|SUSPECT)" "$LOG" 2>/dev/null | \
while IFS= read -r line; do
    status=$(echo "$line"  | awk '{print $1}')
    file=$(echo "$line"    | awk '{print $2}')
    ins=$(echo "$line"     | grep -oP 'ins=\S+')
    te=$(echo "$line"      | grep -oP 'te=\S+')
    source=$(echo "$line"  | grep -oP 'source=\S+')
    gold=$(echo "$line"    | grep -oP 'gold=\S+')
    printf "  %-9s  %-32s  %-20s  %-14s  %s  %s\n" \
        "$status" "$file" "${ins#ins=}" "${te#te=}" "${source#source=}" "${gold}"
done || echo "  (none)"

echo ""
n_ok=$(grep -c "^  OK" "$LOG" 2>/dev/null || echo 0)
n_missing=$(grep -c "^  MISSING" "$LOG" 2>/dev/null || echo 0)
echo "  Written: $n_ok    Missing (gold standard): $n_missing"

# Final summary line from find_te_junctions.py
grep -E "^Phase 1 \(graph" "$LOG" 2>/dev/null || true
grep -E "^  Validated:" "$LOG" 2>/dev/null || true

# ------------------------------------------------------------------
# Show first written junction file as a spot-check
# ------------------------------------------------------------------
first=$(grep "^  OK" "$LOG" 2>/dev/null | head -1 | awk '{print $2}')
if [ -n "$first" ]; then
    echo ""
    echo "=== Spot-check: $first ==="
    cat "$TEST_OUTDIR/$first" 2>/dev/null || cat "$first" 2>/dev/null || true
fi

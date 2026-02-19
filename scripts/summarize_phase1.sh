#!/usr/bin/env bash
# summarize_phase1.sh -- Dump all key Phase 1 diagnostics into one text file.
#
# Run after test_phase1.sh completes. Paste the output file back to Claude
# for analysis.
#
# Usage:
#   bash scripts/summarize_phase1.sh [outdir]
#
# Default outdir: temp_work/test_chr3L_8710861-8744900

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
cd "$PROJECT_ROOT"

source config.sh

TEST_REGION="${REGIONS[0]}"
region_safe="${TEST_REGION//:/_}"
OUTDIR="${1:-temp_work/test_${region_safe}}"

OUT="$OUTDIR/phase1_summary.txt"

{
echo "========================================"
echo "Phase 1 Summary: $TEST_REGION"
echo "Generated: $(date)"
echo "Outdir: $OUTDIR"
echo "========================================"

# --- Gold standard catalog ---
echo ""
echo "=== GOLD STANDARD CATALOG ==="
if [ -s "$OUTDIR/gold_standard.tsv" ]; then
    cat "$OUTDIR/gold_standard.tsv"
else
    echo "(not found or empty)"
fi

echo ""
echo "=== GOLD STANDARD RAW COUNTS ==="
if [ -s "$OUTDIR/gold_standard_reads.tsv" ]; then
    echo "Total gold standard reads: $(wc -l < "$OUTDIR/gold_standard_reads.tsv")"
    echo ""
    echo "Reads per TE:"
    awk '{print $3}' "$OUTDIR/gold_standard_reads.tsv" | sort | uniq -c | sort -rn
    echo ""
    echo "Reads per sample:"
    awk '{print $2}' "$OUTDIR/gold_standard_reads.tsv" | sort | uniq -c | sort -rn
else
    echo "(not found)"
fi

# --- Candidates catalog ---
echo ""
echo "=== CANDIDATES CATALOG ==="
if [ -s "$OUTDIR/candidates_catalog.tsv" ]; then
    echo "Total catalog entries: $(wc -l < "$OUTDIR/candidates_catalog.tsv")"
    echo ""
    echo "Entries per sample:"
    awk '{print $2}' "$OUTDIR/candidates_catalog.tsv" | sort | uniq -c | sort -rn
    echo ""
    echo "Mate contigs (top 20):"
    awk '{print $3}' "$OUTDIR/candidates_catalog.tsv" | sort | uniq -c | sort -rn | head -20
else
    echo "(not found)"
fi

# --- Read counts ---
echo ""
echo "=== READ COUNTS ==="
for f in R1.fq R2.fq singles.fq; do
    if [ -s "$OUTDIR/$f" ]; then
        n=$(grep -c '^@' "$OUTDIR/$f" 2>/dev/null || echo 0)
        echo "$f: $n reads"
    else
        echo "$f: (not found)"
    fi
done

# --- SPAdes assembly ---
echo ""
echo "=== SPADES ASSEMBLY ==="
if [ -s "$OUTDIR/assembly/assembly_graph.fastg" ]; then
    n_nodes=$(grep -c '^>' "$OUTDIR/assembly/assembly_graph.fastg" 2>/dev/null || echo 0)
    echo "assembly_graph.fastg: $n_nodes node headers"
    # Node length distribution
    echo ""
    echo "Node length distribution:"
    awk '/^>/{if(name && len>0) print len; name=$0; len=0; next} {len+=length($0)} END{if(len>0) print len}' \
        "$OUTDIR/assembly/assembly_graph.fastg" | \
    awk '{
        if($1<100) b="<100"
        else if($1<500) b="100-499"
        else if($1<1000) b="500-999"
        else if($1<5000) b="1000-4999"
        else b="5000+"
        counts[b]++
    } END {
        split("<100 100-499 500-999 1000-4999 5000+", order, " ")
        for(i=1;i<=5;i++) printf "  %10s: %d\n", order[i], counts[order[i]]+0
    }'
else
    echo "assembly_graph.fastg: NOT FOUND"
    # Check per-k graphs
    for kval in 21 33 55; do
        f="$OUTDIR/assembly/K${kval}/assembly_graph.fastg"
        if [ -s "$f" ]; then
            echo "  K${kval}/assembly_graph.fastg: $(grep -c '^>' "$f") nodes"
        fi
    done
fi
if [ -s "$OUTDIR/assembly/contigs.fasta" ]; then
    echo ""
    echo "contigs.fasta: $(grep -c '^>' "$OUTDIR/assembly/contigs.fasta") contigs"
fi

# --- BLAST results ---
echo ""
echo "=== BLAST: GRAPH NODES VS TE DATABASE ==="
if [ -s "$OUTDIR/nodes_vs_te.tsv" ]; then
    echo "Total hits: $(wc -l < "$OUTDIR/nodes_vs_te.tsv")"
    echo ""
    echo "Hits per TE (top 20):"
    awk '{print $2}' "$OUTDIR/nodes_vs_te.tsv" | sort | uniq -c | sort -rn | head -20
    echo ""
    echo "Top hits (by alignment length):"
    sort -t$'\t' -k9 -rn "$OUTDIR/nodes_vs_te.tsv" | head -10
else
    echo "(not found)"
fi

echo ""
echo "=== BLAST: GRAPH NODES VS REFERENCE ==="
if [ -s "$OUTDIR/nodes_vs_ref.tsv" ]; then
    echo "Total hits: $(wc -l < "$OUTDIR/nodes_vs_ref.tsv")"
    echo ""
    echo "Top hits (by alignment length):"
    sort -t$'\t' -k9 -rn "$OUTDIR/nodes_vs_ref.tsv" | head -10
else
    echo "(not found)"
fi

# --- Pipeline log (key sections) ---
echo ""
echo "=== PIPELINE LOG (key lines) ==="
if [ -s "$OUTDIR/pipeline.log" ]; then
    grep -E "^(Step |  |Phase 1|=== |Total |Gold |Assembly |Contigs:|WARNING)" \
        "$OUTDIR/pipeline.log" 2>/dev/null || echo "(no matching lines)"
else
    echo "(not found)"
fi

# --- Junction files ---
echo ""
echo "=== JUNCTION FILES ==="
jfiles=$(ls "$OUTDIR"/junction_*.fasta 2>/dev/null)
if [ -n "$jfiles" ]; then
    echo "Count: $(echo "$jfiles" | wc -l | tr -d ' ')"
    echo ""
    for jf in $jfiles; do
        echo "--- $(basename "$jf") ---"
        cat "$jf"
        echo ""
    done
else
    echo "No junction files found"
fi

# --- SPAdes log (last 30 lines, for errors) ---
echo ""
echo "=== SPADES LOG (last 30 lines) ==="
if [ -s "$OUTDIR/spades.log" ]; then
    tail -30 "$OUTDIR/spades.log"
else
    echo "(not found)"
fi

} > "$OUT" 2>&1

echo "Summary written to: $OUT"
echo "Size: $(wc -c < "$OUT") bytes, $(wc -l < "$OUT") lines"
echo ""
echo "Paste it back with:  cat $OUT | pbcopy   (macOS)"
echo "Or on the cluster:   cat $OUT"

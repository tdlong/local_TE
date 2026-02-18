#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Usage: $0 <directory_path> <TE_database_fasta>"
    exit 1
fi

WORK_DIR="$1"
TE_DATABASE="$2"
SCRATCH_DIR="tempexperiment"

# Check files
if [ ! -f "${WORK_DIR}/R1.fq" ] || [ ! -f "${WORK_DIR}/R2.fq" ] || [ ! -f "${WORK_DIR}/region.fasta" ]; then
    echo "Error: Missing required files"
    exit 1
fi

rm -rf "$SCRATCH_DIR"
mkdir -p "$SCRATCH_DIR"

echo "=== TE Junction Reconstruction ==="
echo ""

# Convert to FASTA
cat "${WORK_DIR}/R1.fq" | paste - - - - | awk '{print ">"substr($1,2)"\n"$2}' > "${SCRATCH_DIR}/R1.fa"
cat "${WORK_DIR}/R2.fq" | paste - - - - | awk '{print ">"substr($1,2)"\n"$2}' > "${SCRATCH_DIR}/R2.fa"

# BLAST
blastn -query "${SCRATCH_DIR}/R1.fa" -subject "$TE_DATABASE" -outfmt "6 qseqid sseqid qstart qend sstart send" -evalue 1e-5 > "${SCRATCH_DIR}/R1_vs_TE.txt" 2>/dev/null
blastn -query "${SCRATCH_DIR}/R2.fa" -subject "$TE_DATABASE" -outfmt "6 qseqid sseqid qstart qend sstart send" -evalue 1e-5 > "${SCRATCH_DIR}/R2_vs_TE.txt" 2>/dev/null
blastn -query "${SCRATCH_DIR}/R1.fa" -subject "${WORK_DIR}/region.fasta" -outfmt "6 qseqid qstart qend sstart send" -evalue 1e-5 > "${SCRATCH_DIR}/R1_vs_ref.txt" 2>/dev/null
blastn -query "${SCRATCH_DIR}/R2.fa" -subject "${WORK_DIR}/region.fasta" -outfmt "6 qseqid qstart qend sstart send" -evalue 1e-5 > "${SCRATCH_DIR}/R2_vs_ref.txt" 2>/dev/null

# Find junction reads
cat "${SCRATCH_DIR}"/R*_vs_TE.txt | awk '{print $1}' | sort -u > "${SCRATCH_DIR}/te_reads.txt"
cat "${SCRATCH_DIR}"/R*_vs_ref.txt | awk '{print $1}' | sort -u > "${SCRATCH_DIR}/ref_reads.txt"
join "${SCRATCH_DIR}/te_reads.txt" "${SCRATCH_DIR}/ref_reads.txt" > "${SCRATCH_DIR}/junction_reads.txt"

N_JUNCTION=$(wc -l < "${SCRATCH_DIR}/junction_reads.txt")
echo "Found $N_JUNCTION junction reads"
echo ""

# Load reference
REF_SEQ=$(grep -v "^>" "${WORK_DIR}/region.fasta" | tr -d '\n')

# Analyze each junction read
echo "Analyzing junctions..."
> "${SCRATCH_DIR}/junctions.txt"

while read READ_ID; do
    # Get BLAST hits
    REF=$(grep "^${READ_ID}" "${SCRATCH_DIR}"/R*_vs_ref.txt | head -1)
    TE=$(grep "^${READ_ID}" "${SCRATCH_DIR}"/R*_vs_TE.txt | head -1)
    
    if [ -z "$REF" ] || [ -z "$TE" ]; then
        continue
    fi
    
    # Parse alignments (qstart qend sstart send)
    REF_QSTART=$(echo "$REF" | awk '{print $2}')
    REF_QEND=$(echo "$REF" | awk '{print $3}')
    REF_SSTART=$(echo "$REF" | awk '{print $4}')
    REF_SEND=$(echo "$REF" | awk '{print $5}')
    
    TE_QSTART=$(echo "$TE" | awk '{print $3}')
    TE_QEND=$(echo "$TE" | awk '{print $4}')
    TE_SSTART=$(echo "$TE" | awk '{print $5}')
    TE_SEND=$(echo "$TE" | awk '{print $6}')
    TE_NAME=$(echo "$TE" | awk '{print $2}')
    
    # Determine junction type based on read coordinates
    # LEFT junction: reference comes first (lower q), then TE
    # RIGHT junction: TE comes first (lower q), then reference
    
    if [ $REF_QSTART -lt $TE_QSTART ]; then
        # LEFT junction: ref...TE
        # Overlap bases go to TE, so junction is where TE starts
        # Map TE_QSTART back to reference position
        SIDE="LEFT"
        
        # Calculate reference position at TE_QSTART
        # If ref is forward: ref_pos = REF_SSTART + (TE_QSTART - REF_QSTART - 1)
        # If ref is reverse: ref_pos = REF_SSTART - (TE_QSTART - REF_QSTART - 1)
        if [ $REF_SSTART -lt $REF_SEND ]; then
            # Forward
            POS=$((REF_SSTART + TE_QSTART - REF_QSTART - 1))
        else
            # Reverse
            POS=$((REF_SSTART - TE_QSTART + REF_QSTART + 1))
        fi
        
    else
        # RIGHT junction: TE...ref
        # Overlap bases go to reference, so junction is where reference starts
        SIDE="RIGHT"
        
        if [ $REF_SSTART -lt $REF_SEND ]; then
            POS=$REF_SSTART
        else
            POS=$REF_SEND
        fi
    fi
    
    echo "$POS $SIDE $TE_NAME $TE_SSTART $TE_SEND" >> "${SCRATCH_DIR}/junctions.txt"
    
done < "${SCRATCH_DIR}/junction_reads.txt"

# Cluster by position (within 20bp)
sort -n "${SCRATCH_DIR}/junctions.txt" | awk '
BEGIN { cluster=0; prev_pos=-1000; }
{
    if ($1 - prev_pos > 20) {
        cluster++;
    }
    print cluster, $0;
    prev_pos = $1;
}' > "${SCRATCH_DIR}/clustered.txt"

N_CLUSTERS=$(awk '{print $1}' "${SCRATCH_DIR}/clustered.txt" | sort -u | wc -l)
echo "Found $N_CLUSTERS distinct insertion sites"
echo ""
echo "========================================"
echo ""

# Process each cluster
for CLUSTER in $(awk '{print $1}' "${SCRATCH_DIR}/clustered.txt" | sort -u); do
    
    # Get cluster data
    awk -v c=$CLUSTER '$1==c' "${SCRATCH_DIR}/clustered.txt" > "${SCRATCH_DIR}/cluster_${CLUSTER}.txt"
    
    # Get most common TE
    TE=$(awk '{print $4}' "${SCRATCH_DIR}/cluster_${CLUSTER}.txt" | sort | uniq -c | sort -rn | head -1 | awk '{print $2}')
    TE_COUNT=$(awk '{print $4}' "${SCRATCH_DIR}/cluster_${CLUSTER}.txt" | sort | uniq -c | sort -rn | head -1 | awk '{print $1}')
    
    # Filter to this TE
    grep "$TE" "${SCRATCH_DIR}/cluster_${CLUSTER}.txt" > "${SCRATCH_DIR}/cluster_${CLUSTER}_filt.txt"
    
    # Count left/right
    N_LEFT=$(grep "LEFT" "${SCRATCH_DIR}/cluster_${CLUSTER}_filt.txt" | wc -l)
    N_RIGHT=$(grep "RIGHT" "${SCRATCH_DIR}/cluster_${CLUSTER}_filt.txt" | wc -l)
    
    # Skip if not enough support
    if [ $TE_COUNT -lt 3 ] || [ $(($N_LEFT + $N_RIGHT)) -lt 3 ]; then
        continue
    fi
    
    # Get TE sequence
    TE_SEQ=$(grep -A1 "^>$TE" "$TE_DATABASE" | grep -v "^>" | tr -d '\n')
    TE_LEN=${#TE_SEQ}
    
    echo "INSERTION: $TE (length ${TE_LEN}bp)"
    echo "Supporting reads: $N_LEFT left, $N_RIGHT right"
    echo ""
    
    # Left junction
    if [ $N_LEFT -gt 0 ]; then
        LEFT_POS=$(grep "LEFT" "${SCRATCH_DIR}/cluster_${CLUSTER}_filt.txt" | awk '{sum+=$2; n++} END {print int(sum/n)}')
        echo "  Left junction at position $LEFT_POS:"
        
        # Extract 100bp reference before junction
        UPSTREAM=${REF_SEQ:$((LEFT_POS-101)):100}
        
        # TE starts at beginning
        TE_START=${TE_SEQ:0:100}
        
        echo "    Ref: ...${UPSTREAM}|"
        echo "    TE:  |${TE_START}..."
        echo ""
    fi
    
    # Right junction  
    if [ $N_RIGHT -gt 0 ]; then
        RIGHT_POS=$(grep "RIGHT" "${SCRATCH_DIR}/cluster_${CLUSTER}_filt.txt" | awk '{sum+=$2; n++} END {print int(sum/n)}')
        echo "  Right junction at position $RIGHT_POS:"
        
        # TE ends at end
        TE_END=${TE_SEQ:$((TE_LEN-100)):100}
        
        # Extract 100bp reference after junction
        DOWNSTREAM=${REF_SEQ:$((RIGHT_POS-1)):100}
        
        echo "    TE:  ...${TE_END}|"
        echo "    Ref: |${DOWNSTREAM}..."
        echo ""
    fi
    
    # Show position difference
    if [ $N_LEFT -gt 0 ] && [ $N_RIGHT -gt 0 ]; then
        DIFF=$((RIGHT_POS - LEFT_POS))
        echo "  Position difference: ${DIFF}bp (likely TSD or sequence overlap)"
        echo ""
    fi
    
    echo "----------------------------------------"
    echo ""
done

echo "Done!"

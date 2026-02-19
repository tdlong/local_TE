#!/usr/bin/env bash
# config.sh -- Edit this file to configure both pipeline phases
#
# This is the only file you need to edit to run the pipeline.
# Both Phase 1 (discovery) and Phase 2 (frequency) read from here.

#=============================================================================
# Reference files
#=============================================================================
REF="/dfs7/adl/sruckman/XQTL/XQTL2/ref/dm6_conTE.fa"
TEFASTA="/dfs7/adl/sruckman/XQTL/XQTL2/ref/transposon_sequence_set.fa"

#=============================================================================
# Phase 1: Discovery -- BAMs and regions
#=============================================================================
BAM_DIR="/dfs7/adl/sruckman/XQTL/XQTL2/data/bam/colorTEcon"
BAMS=("HOULE_L2F.bam" "HOULE_L3F.bam" "HOUSTON_L2F.bam" "HOUSTON_L3F.bam" "HOULE_L1F.bam" "HOUSTON_L1F.bam")
REGIONS=("chr3L:8710861-8744900" "chrX:8450000-8500000")

#=============================================================================
# Phase 2: Frequency -- raw FASTQs (tab-delimited: sample R1 R2)
#=============================================================================
SAMPLES=(
"HOULE_L2F	/dfs7/adl/sruckman/XQTL/XQTL2/data/raw/Dmel_color/xR065-L1-G1-P02-GTCCTAGA-GGTACTTC-READ1-Sequences.txt.gz	/dfs7/adl/sruckman/XQTL/XQTL2/data/raw/Dmel_color/xR065-L1-G1-P02-GTCCTAGA-GGTACTTC-READ2-Sequences.txt.gz"
"HOULE_L3F	/dfs7/adl/sruckman/XQTL/XQTL2/data/raw/Dmel_color/xR065-L1-G1-P03-ACTTGCCA-GGTACTTC-READ1-Sequences.txt.gz	/dfs7/adl/sruckman/XQTL/XQTL2/data/raw/Dmel_color/xR065-L1-G1-P03-ACTTGCCA-GGTACTTC-READ2-Sequences.txt.gz"
"HOUSTON_L2F	/dfs7/adl/sruckman/XQTL/XQTL2/data/raw/Dmel_color/xR065-L1-G1-P05-CTGTGCTT-GGTACTTC-READ1-Sequences.txt.gz	/dfs7/adl/sruckman/XQTL/XQTL2/data/raw/Dmel_color/xR065-L1-G1-P05-CTGTGCTT-GGTACTTC-READ2-Sequences.txt.gz"
"HOUSTON_L3F	/dfs7/adl/sruckman/XQTL/XQTL2/data/raw/Dmel_color/xR065-L1-G1-P06-TCGAACTC-GGTACTTC-READ1-Sequences.txt.gz	/dfs7/adl/sruckman/XQTL/XQTL2/data/raw/Dmel_color/xR065-L1-G1-P06-TCGAACTC-GGTACTTC-READ2-Sequences.txt.gz"
"HOULE_L1	/dfs7/adl/sruckman/XQTL/XQTL2/data/raw/Dmel_color/xR065-L1-G1-P01-TGGCTATG-GGTACTTC-READ1-Sequences.txt.gz	/dfs7/adl/sruckman/XQTL/XQTL2/data/raw/Dmel_color/xR065-L1-G1-P01-TGGCTATG-GGTACTTC-READ2-Sequences.txt.gz"
"HOUSTON_L1F	/dfs7/adl/sruckman/XQTL/XQTL2/data/raw/Dmel_color/xR065-L1-G1-P04-TCCACTCA-GGTACTTC-READ1-Sequences.txt.gz	/dfs7/adl/sruckman/XQTL/XQTL2/data/raw/Dmel_color/xR065-L1-G1-P04-TCCACTCA-GGTACTTC-READ2-Sequences.txt.gz"
)

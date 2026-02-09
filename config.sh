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
BAMS=(
"HOULE_L2F.bam"
)

REGIONS=(
"chr3L:8710861-8744900"
)

#=============================================================================
# Phase 2: Frequency -- raw FASTQs
# Tab-delimited: sample_name	R1_path	R2_path
#=============================================================================
SAMPLES=(
"HOULE_L2F	/path/to/HOULE_L2F.F.fq.gz	/path/to/HOULE_L2F.R.fq.gz"
)

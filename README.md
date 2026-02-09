# TE Junction Assembly & Frequency Estimation Pipeline

A two-phase pipeline for detecting transposable element (TE) insertions from paired-end sequencing data and estimating their population frequency.

## Overview

**Phase 1 (Discovery):** Pool reads from all BAMs for each genomic region, assemble with SPAdes, identify TE junction contigs via BLAST/minimap2, and generate visual alignments.

**Phase 2 (Frequency):** Build a competitive alignment reference from all discovered junctions, align each sample's raw FASTQs, and count reads spanning each junction to estimate TE insertion frequency per sample.

## Requirements

HPC modules:
- `samtools/1.15.1`
- `SPAdes/3.15.4`
- `ncbi-blast/2.13.0`
- `minimap2/2.28`
- `bwa/0.7.17`
- `python/3.10.2` with Biopython

## Configuration

All user-editable settings live in **`config.sh`** in the project root. This is the only file you need to edit.

```bash
# Reference files
REF="/dfs7/adl/sruckman/XQTL/XQTL2/ref/dm6_conTE.fa"
TEFASTA="/dfs7/adl/sruckman/XQTL/XQTL2/ref/transposon_sequence_set.fa"

# Phase 1: Discovery -- BAMs and regions
BAM_DIR="/dfs7/adl/sruckman/XQTL/XQTL2/data/bam/colorTEcon"
BAMS=("HOULE_L2F.bam" "sample2.bam" "sample3.bam")
REGIONS=("chr3L:8710861-8744900" "chr2R:1000000-1050000")

# Phase 2: Frequency -- raw FASTQs (tab-delimited: sample  R1  R2)
SAMPLES=(
"HOULE_L2F	/path/to/HOULE_L2F.F.fq.gz	/path/to/HOULE_L2F.R.fq.gz"
"sample2	/path/to/sample2.F.fq.gz	/path/to/sample2.R.fq.gz"
)
```

## Quick Start

```bash
cd /dfs7/adl/tdlong/Sarah/local_TE

# 1. Edit config.sh with your BAMs, regions, and FASTQ paths

# 2. Phase 1: Discover TE insertions
sbatch scripts/submit_te_analysis.sh

# 3. Build competitive reference from discovered junctions
module load python/3.10.2
python scripts/build_junctions_ref.py temp_work junctions.fa

# 4. Phase 2: Estimate frequencies from raw FASTQs
sbatch scripts/submit_te_frequency.sh
```

## Phase 1: Discovery

All BAMs are pooled per region to maximize assembly coverage. Output is one directory per region (not per BAM).

### Interactive Mode

```bash
srun --pty bash
bash scripts/run_te_assembly.sh \
    "chr3L:8710861-8744900" \
    "temp_work/chr3L_8710861-8744900" \
    "/path/to/sample1.bam" \
    "/path/to/sample2.bam"
```

### Arguments

| Argument | Description | Example |
|----------|-------------|---------|
| `region` | Genomic region to analyze | `chr3L:8710861-8744900` |
| `outdir` | Output directory | `temp_work/chr3L_8710861-8744900` |
| `bam...` | One or more BAM files to pool | `/path/to/sample1.bam /path/to/sample2.bam` |

### Environment Variables (Optional)

| Variable | Description | Default |
|----------|-------------|---------|
| `REF` | Reference genome FASTA | `/dfs7/adl/sruckman/XQTL/XQTL2/ref/dm6_conTE.fa` |
| `TEFASTA` | TE sequence database | `/dfs7/adl/sruckman/XQTL/XQTL2/ref/transposon_sequence_set.fa` |

### Discovery Output

```
temp_work/
└── chr3L_8710861-8744900/
    ├── pipeline.log
    ├── region.fasta
    ├── reads.fasta              # Pooled reads from all BAMs
    ├── assembly/
    ├── te_contigs.fasta
    ├── te_seqs.fasta
    ├── junctions_to_ref.paf
    ├── junctions_to_te.paf
    └── junction_*.fasta         # Visual alignments (one per TE insertion)
```

### Output Visualization

Each junction file shows a 100bp view of the insertion:

```
INSERTION of FBte0000626 at chr3L:8711446
======================================================================
  Junction: NODE_21_length_563_cov_4.021654
    Type: right, transition at junc position 237

    100bp view (insertion at position 50):
    WT_REF: AGTGCCGAAAGTACAAGTTAAGTACATACATCGTGCCACTATTAACGCTCCACTGACAGCGGCAAAACACGCATCAAAAACACACATACAAATCGGCAGA
    REF:    --------------------------------------------------CACTGACAGCGGCAAAACACGCATCAAAAACACACATACAAATCGGCAGA
    JUNC:   GGTCATCATTTCGAATTTCTGCCAAAAAAAACGCATAAAAAACCACTGTGCACTGACAGCGGCAAAACACGCATCAAAAACACACATACAAATCGGCAGA
    TE:     GTCATCATTTCGAATTTCTGCCAAAAAAAAACACATAAAAAACCACTGTG--------------------------------------------------
```

- **WT_REF**: Wild-type reference (continuous, no TE)
- **REF**: Reference portion aligned to the junction (gap where TE is)
- **JUNC**: Assembled junction contig spanning the insertion
- **TE**: Transposable element sequence

## Building the Competitive Reference

After Phase 1, build `junctions.fa` from all discovered junctions:

```bash
python scripts/build_junctions_ref.py temp_work junctions.fa
```

This scans all `junction_*.fasta` files, extracts Abs (absence/WT_REF) and Pre (presence/JUNC) sequences, deduplicates by genomic position + TE name, and writes:

- `junctions.fa` -- competitive alignment reference with paired entries:
  ```
  >chr3L_8711446_FBte0000626_Abs
  AGTGCCGAAA...reference...AATCGGCAGA
  >chr3L_8711446_FBte0000626_Pre
  GGTCATCATTT...TE junction...AATCGGCAGA
  ```
- `junctions_metadata.tsv` -- metadata for each junction pair

## Phase 2: Frequency Estimation

### Why Raw FASTQs?

Extracting reads from BAMs introduces bias: the whole-genome aligner has already made decisions about multimappers and unmapped reads. A read that truly spans a TE junction might have been discarded as a multimapper. Going back to raw FASTQs and aligning to the small `junctions.fa` reference lets bwa make a clean, unbiased competitive alignment.

### Sample FASTQs

Fill in the `SAMPLES` array in `config.sh` (tab-delimited: sample_name, R1, R2). The submit script writes these to `samples.tsv` automatically.

### Running

```bash
# SLURM (builds junctions.fa if needed, writes samples.tsv from config, runs estimation)
sbatch scripts/submit_te_frequency.sh

# Or manually:
bash scripts/run_te_frequency.sh junctions.fa samples.tsv freq_work
```

### How It Works

1. `bwa index junctions.fa` (once)
2. For each sample: `bwa mem -t 8 -W 15 -a junctions.fa R1.fq.gz R2.fq.gz`
3. Count reads spanning each junction midpoint (>= 10bp anchor on both sides)
4. Report frequency = pre_count / (abs_count + pre_count)

The `-W 15` flag ensures bwa doesn't ignore reads with short overlaps near the junction edge. The spanning filter (`pos <= 40 AND end >= 60`) discards reads that map only to the shared 50bp flank, which would otherwise be ambiguous.

### Frequency Output

```
freq_work/
├── frequency_results.tsv     # Combined results for all samples
├── HOULE_L2F.sam             # Per-sample alignment (can delete to save space)
├── HOULE_L2F_bwa.log
├── sample2.sam
└── ...
```

`frequency_results.tsv` format:

```
sample      junction_id                     abs_count  pre_count  frequency
HOULE_L2F   chr3L_8711446_FBte0000626       42         18         0.3000
sample2     chr3L_8711446_FBte0000626       55         5          0.0833
```

## Pipeline Sections (Phase 1)

1. **Read Extraction**: Pools "gold" reads (region + FBte) and junction candidates (mate unmapped) from all BAMs
2. **Assembly**: SPAdes assembly of pooled junction-spanning reads
3. **BLAST**: Identifies contigs with TE similarity
4. **TE Contig Extraction**: Extracts full sequences of TE-hitting contigs
5. **TE Sequence Extraction**: Gets canonical TE sequences from database
6. **Minimap2 Alignment**: Aligns contigs to reference and TE sequences
7. **Visualization**: Generates visual alignments with `build_te_alignment.py`

## Files

```
local_TE/
├── README.md
├── Notes.txt
├── config.sh                            # Edit this: BAMs, regions, FASTQs, refs
├── junctions.fa                         # Built by build_junctions_ref.py
├── junctions_metadata.tsv               # Metadata for each junction
└── scripts/
    ├── submit_te_analysis.sh            # Phase 1 SLURM wrapper
    ├── run_te_assembly.sh               # Phase 1 main pipeline
    ├── build_te_alignment.py            # Phase 1 visualization
    ├── build_junctions_ref.py           # Bridge: build competitive reference
    ├── submit_te_frequency.sh           # Phase 2 SLURM wrapper
    ├── run_te_frequency.sh              # Phase 2 main pipeline
    └── estimate_te_freq.py              # Phase 2 frequency counter
```

## Troubleshooting

**No contigs produced (Phase 1):**
- Check that input BAMs have reads in the specified region
- Verify the region contains a TE insertion

**No BLAST hits:**
- The assembled contigs may not span the TE junction
- Try a different/larger region

**No spanning reads (Phase 2):**
- Check that FASTQ files are correct for the samples
- Verify `junctions.fa` contains the expected Abs/Pre pairs
- Check bwa log files in `freq_work/`

**Missing metadata in junction headers:**
- Ensure Phase 1 was run with the updated `build_te_alignment.py` that embeds `insertion=` and `te=` in WT_REF headers

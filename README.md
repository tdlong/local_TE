# TE Junction Assembly & Genotyping Pipeline

A two-phase pipeline for detecting transposable element (TE) insertions from paired-end sequencing data and genotyping them across samples.

## Overview

**Phase 1 (Discovery):** Pool reads from all BAMs for each genomic region, assemble with SPAdes, identify TE junction contigs via BLAST/minimap2, and generate visual alignments.

**Phase 2 (Genotyping):** Extract diagnostic k-mers from junction sequences, scan raw FASTQs with BBDuk, and call genotypes from k-mer counts.

## Requirements

HPC modules:
- `samtools/1.15.1`
- `SPAdes/3.15.4`
- `ncbi-blast/2.13.0`
- `minimap2/2.28`
- `python/3.10.2` with Biopython

Conda:
- `bbmap` (BBDuk for k-mer counting): `conda install -c bioconda bbmap`

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

# Phase 2: Genotyping -- raw FASTQs (tab-delimited: sample  R1  R2)
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

# 3. Phase 2: Genotype via k-mer counting
sbatch scripts/submit_te_kmer_count.sh
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

- `junctions.fa` -- paired entries:
  ```
  >chr3L_8711446_FBte0000626_Abs
  AGTGCCGAAA...reference...AATCGGCAGA
  >chr3L_8711446_FBte0000626_Pre
  GGTCATCATTT...TE junction...AATCGGCAGA
  ```
- `junctions_metadata.tsv` -- metadata for each junction pair

## Phase 2: Genotyping via K-mer Counting

### Why K-mers Instead of Alignment?

Aligning hundreds of millions of paired-end reads against a tiny `junctions.fa` reference with bwa is wasteful and introduces alignment artifacts. Instead, we extract short diagnostic k-mers (31bp) that span the junction point and are unique to each allele, then use BBDuk to efficiently scan raw FASTQs for those k-mers.

### How It Works

1. **Extract diagnostic k-mers** (`extract_junction_kmers.py`):
   - From each Abs/Pre pair in `junctions.fa`, extract k=31 k-mers spanning the junction at position 50
   - Keep only allele-specific k-mers (set difference: abs - pre, pre - abs)
   - Write `abs_kmers.fa` and `pre_kmers.fa`

2. **Count k-mers per sample** (`run_te_kmer_count.sh`):
   - BBDuk scans each sample's paired FASTQs against the k-mer databases
   - Two runs per sample: once for absence k-mers, once for presence k-mers
   - Flags: `k=31 maskmiddle=f rcomp=t`
   - Output: `<sample>_abs_stats.txt` and `<sample>_pre_stats.txt`

3. **Call genotypes** (`genotype_from_counts.py`):
   - Parse `#Matched` read counts from BBDuk stats files
   - Compute abs_ratio = abs_count / (abs_count + pre_count)
   - Call: Abs/Abs (ratio >= 0.8), Pre/Pre (ratio <= 0.2), Abs/Pre (between), or NO_CALL (depth < 5)

### Running

```bash
# Automated (builds k-mers, submits array job, submits genotyping with dependency)
sbatch scripts/submit_te_kmer_count.sh

# Or step by step:
python scripts/extract_junction_kmers.py junctions.fa kmer_work
# then run BBDuk per sample (see run_te_kmer_count.sh)
python scripts/genotype_from_counts.py kmer_work
```

### Genotyping Output

```
kmer_work/
├── abs_kmers.fa              # Absence allele diagnostic k-mers
├── pre_kmers.fa              # Presence allele diagnostic k-mers
├── HOULE_L2F_abs_stats.txt   # BBDuk stats per sample
├── HOULE_L2F_pre_stats.txt
├── sample2_abs_stats.txt
├── sample2_pre_stats.txt
└── genotype_results.tsv      # Final genotype calls
```

`genotype_results.tsv` format:

```
sample      genotype  abs_reads  pre_reads  total  notes
HOULE_L2F   Abs/Pre   420        180        600    abs_ratio=0.700
sample2     Abs/Abs   550        50         600    abs_ratio=0.917
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
    ├── submit_te_kmer_count.sh          # Phase 2 SLURM orchestrator
    ├── run_te_kmer_count.sh             # Phase 2 BBDuk per-sample
    ├── extract_junction_kmers.py        # Phase 2 k-mer extraction
    └── genotype_from_counts.py          # Phase 2 genotyping
```

## Troubleshooting

**No contigs produced (Phase 1):**
- Check that input BAMs have reads in the specified region
- Verify the region contains a TE insertion

**No BLAST hits:**
- The assembled contigs may not span the TE junction
- Try a different/larger region

**No k-mer matches (Phase 2):**
- Check that FASTQ files are correct for the samples
- Verify `junctions.fa` contains the expected Abs/Pre pairs
- Inspect `abs_kmers.fa` / `pre_kmers.fa` to confirm diagnostic k-mers were extracted
- Check BBDuk stats files for total read counts

**Missing metadata in junction headers:**
- Ensure Phase 1 was run with the updated `build_te_alignment.py` that embeds `insertion=` and `te=` in WT_REF headers

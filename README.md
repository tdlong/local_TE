# TE Junction Assembly Pipeline

A pipeline for detecting and visualizing transposable element (TE) insertions from paired-end sequencing data. The pipeline extracts reads spanning TE junctions, performs local assembly, and generates visual alignments showing the insertion points.

## Overview

The pipeline identifies TE insertions by:
1. Extracting read pairs where one read maps to a genomic region of interest and the other maps to a TE (or has an unmapped mate suggesting a junction)
2. Assembling these reads with SPAdes
3. Identifying contigs that hit the TE database via BLAST
4. Aligning junction contigs to both reference and TE sequences with minimap2
5. Generating a visual alignment showing the insertion point

## Requirements

The following modules are expected to be available (HPC environment):
- `samtools/1.17`
- `SPAdes/3.15.4`
- `ncbi-blast/2.13.0`
- `minimap2/2.28`
- `python/3.10.2` with Biopython

## Usage

### SLURM Batch Mode (Recommended)

```bash
cd /dfs7/adl/tdlong/Sarah/local_TE
sbatch scripts/submit_te_analysis.sh
```

Edit `scripts/submit_te_analysis.sh` to configure BAMs and regions:

```bash
BAMS=("HOULE_L2F.bam" "sample2.bam" "sample3.bam" "sample4.bam")
REGIONS=("chr3L:8710861-8744900" "chr2R:1000000-1050000")
```

### Interactive Mode

```bash
cd /dfs7/adl/tdlong/Sarah/local_TE
srun --pty bash
bash scripts/run_te_assembly.sh \
    "chr3L:8710861-8744900" \
    "/dfs7/adl/sruckman/XQTL/XQTL2/data/bam/colorTEcon/HOULE_L2F.bam" \
    "temp_work/test_run"
```

### Arguments

| Argument | Description | Example |
|----------|-------------|---------|
| `region` | Genomic region to analyze | `chr3L:8710861-8744900` |
| `bam` | Path to input BAM file | `/path/to/alignments.bam` |
| `outdir` | Output directory for all files | `temp_work/sample_region` |

### Environment Variables (Optional)

| Variable | Description | Default |
|----------|-------------|---------|
| `REF` | Reference genome FASTA | `/dfs7/adl/sruckman/XQTL/XQTL2/ref/dm6_conTE.fa` |
| `TEFASTA` | TE sequence database | `/dfs7/adl/sruckman/XQTL/XQTL2/ref/transposon_sequence_set.fa` |

## Output Files

All output files are created in the specified `scratch_dir`:

| File | Description |
|------|-------------|
| `region.fasta` | Reference sequence for the region |
| `candidate_reads.bam` | Extracted junction-spanning reads |
| `te_contigs.fasta` | Assembled contigs with TE BLAST hits |
| `junctions_to_ref.paf` | Minimap2 alignment of contigs to reference |
| `junctions_to_te.paf` | Minimap2 alignment of contigs to TEs |
| `junction_*.fasta` | Visual alignment files for each insertion |

## Output Visualization

The pipeline generates a 100bp visualization for each detected TE insertion:

```
INSERTION of FBte0000626 at chr3L:8711446
======================================================================
  Junction: NODE_21_length_563_cov_4.021654
  RC=False
  junc[237:563] matches ref[585:911]
  junc[0:237] matches TE[1923:2167]
    Type: right, transition at junc position 237

    100bp view (insertion at position 50):
    WT_REF: AGTGCCGAAAGTACAAGTTAAGTACATACATCGTGCCACTATTAACGCTCCACTGACAGCGGCAAAACACGCATCAAAAACACACATACAAATCGGCAGA
    REF:    --------------------------------------------------CACTGACAGCGGCAAAACACGCATCAAAAACACACATACAAATCGGCAGA
    JUNC:   GGTCATCATTTCGAATTTCTGCCAAAAAAAACGCATAAAAAACCACTGTGCACTGACAGCGGCAAAACACGCATCAAAAACACACATACAAATCGGCAGA
    TE:     GTCATCATTTCGAATTTCTGCCAAAAAAAAACACATAAAAAACCACTGTG--------------------------------------------------
```

**Reading the visualization:**
- **WT_REF**: Wild-type reference (continuous sequence without insertion)
- **REF**: Reference portion that aligns to the junction (with gap for TE)
- **JUNC**: The assembled junction contig spanning the insertion
- **TE**: The transposable element sequence

The transition point (position 50) marks where the reference ends and the TE begins (or vice versa for right junctions).

## Pipeline Sections

1. **Read Extraction**: Extracts "gold standard" reads (region + TE) and junction candidates (mate unmapped)
2. **Assembly**: SPAdes assembly of junction-spanning reads
3. **BLAST**: Identifies contigs with TE similarity
4. **TE Contig Extraction**: Extracts full sequences of TE-hitting contigs
5. **TE Sequence Extraction**: Gets canonical TE sequences from database
6. **Minimap2 Alignment**: Aligns contigs to reference and TE sequences
7. **Visualization**: Generates visual alignments with `build_te_alignment.py`

## Troubleshooting

**No contigs produced:**
- Check that input BAM has reads in the specified region
- Verify the region contains a TE insertion

**No BLAST hits:**
- The assembled contigs may not span the TE junction
- Try a different/larger region

**Empty junction files:**
- The minimap2 alignments may not meet quality thresholds
- Check `junctions_to_ref.paf` and `junctions_to_te.paf` manually

## Files

```
assemble_TE/
├── README.md
├── Notes.txt                    # Original planning notes
└── scripts/
    ├── submit_te_analysis.sh    # SLURM wrapper (edit BAMs/regions here)
    ├── run_te_assembly.sh       # Main pipeline script
    └── build_te_alignment.py    # Visualization script
```

## Output Structure

After running, results are in `temp_work/`:

```
temp_work/
└── HOULE_L2F_chr3L_8710861-8744900/
    ├── pipeline.log             # Full pipeline output
    ├── region.fasta
    ├── reads.fasta
    ├── assembly/
    ├── te_contigs.fasta
    ├── te_seqs.fasta
    ├── junctions_to_ref.paf
    ├── junctions_to_te.paf
    └── junction_*.fasta         # Visual alignments
```

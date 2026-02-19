# TE Junction Assembly & Genotyping Pipeline

A two-phase pipeline for detecting transposable element (TE) insertions from paired-end sequencing data and genotyping them across samples.

## Overview

**Phase 1 (Discovery):** Pool reads from all BAMs for each genomic region, assemble with SPAdes (metagenome mode), parse the assembly graph for TE-containing nodes, extend TE boundaries into reference sequence via k-mer walking, cross-validate against a gold standard catalog, and write junction files.

**Phase 2 (Genotyping):** Extract diagnostic k-mers from junction sequences, scan raw FASTQs with BBDuk, and call genotypes from k-mer counts.

## Requirements

HPC modules:
- `samtools/1.15.1`
- `ncbi-blast/2.13.0`
- `SPAdes/3.15.4`
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

# 2. Phase 1: Discover TE insertions (pool BAMs, assemble, walk junctions)
sbatch scripts/submit_te_analysis.sh
# → temp_work/<region>/junction_*.fasta   (one file per TE found, per region)
# → junctions.fa                          (Abs/Pre sequence pairs for all TEs)

# 3. Phase 2: Genotype each TE across all samples via k-mer counting
sbatch scripts/submit_te_kmer_count.sh
# → kmer_work/genotype_results.tsv        (one row per sample × TE junction)
```

## Phase 1: Discovery

### Pipeline Sections

1. **Section 1 — Extract Candidate Reads:** Pools junction-candidate reads (mate not on same chromosome) from all BAMs. Outputs `candidates_catalog.tsv` preserving sample identity and mate contig for downstream gold standard analysis.

2. **Section 2 — Paired Reads + Gold Standard:** Extracts paired-end FASTQs (R1.fq, R2.fq) and reference region. Builds a **gold standard catalog** (`gold_standard.tsv`) from reads whose mates map to FBte* contigs — this tells us before any assembly which TEs are expected, approximate positions, and per-sample read counts.

3. **Section 3 — SPAdes Assembly:** Assembles reads with `spades.py --meta -k 21,33,55`. We use the **assembly graph** (`assembly_graph.fastg`), not `contigs.fasta`, because the graph preserves all edges including low-coverage TE junction paths.

4. **Section 4 — Find TE Junctions** (`find_te_junctions.py`):
   - Parse graph nodes from FASTG
   - BLAST graph nodes vs TE database and reference region
   - Classify nodes: junction (both hits), TE-boundary (near TE terminus), TE-interior, ref-only
   - K-mer walk from TE-boundary nodes into reference (~200bp extension)
   - Cross-validate discovered junctions against gold standard catalog
   - Write `junction_{type}_{N}.fasta` per junction

### Gold Standard Catalog

The gold standard (`gold_standard.tsv`) is built automatically as a pipeline step from Section 1 output. It summarizes reads whose mates map to FBte* contigs:

```
te_name        approx_pos  total_reads  HOULE_L1  HOULE_L2F  HOUSTON_L1F  ...
FBte0000559    8720500     45           8          12          25
FBte0000301    8735100     12           0          12          0
```

This provides ground truth for validating junction discovery results: which TEs should be found, approximate positions, and which populations contain them.

### Interactive Mode

```bash
srun --pty bash
bash scripts/test_phase1.sh
```

### Junction File Format

Each `junction_{type}_{N}.fasta` contains 4 records, each exactly 100 bp, with the
insertion point at position 50 (0-indexed).

**Case convention:** reference-derived bases are written **lowercase**; TE-derived
bases are written **UPPERCASE**. The case boundary in the junction record marks the
insertion point exactly.

#### The 4 Records

**Record 1 — WT_REF (absence allele)**

```
>WT_REF[8733808:8733907] insertion=chr3L:8733858 te=FBte0000559 type=right
aagcggcgcacacgggtggtggtctgctgggagacaccctcctgctcggacagctggcggcggtagatgttgatcttggcagtggacttgtccagcgcgt
```

100 bp of the reference genome centered on the insertion site. All lowercase
(entirely reference sequence). The genomic coordinates, insertion coordinate, and
TE name are encoded in the header.

**Record 2 — REF (duplicate)**

Identical to WT_REF. Included for alignment visualization tools.

**Record 3 — walked_junction (presence allele)**

```
>walked_junction_right_0
CCCTCCTGCTCGGCCATGCTACGTACGTACGTACGTACGTACGTACGTAcagctggcggcggtagatgttgatcttggcagtggacttgtccagcgcgt
```

Junction sequence discovered from SPAdes graph node + k-mer walk. UPPERCASE = TE-derived
bases, lowercase = reference-derived bases.

**Record 4 — TE canonical sequence**

```
>FBte0000559
CGAGCGGAAAGACAGCAATTTTGGCCGTCACCAAAAAAGTGGCTGCATAGTGCCAAACCAATGTATGGCCGTTACGCATCTTGTTATTCTAGTGTCTTTG
```

100 bp of the canonical TE sequence from the database for comparison.

#### LEFT vs. RIGHT Junction Types

```
RIGHT junction (type=right):
  Sequence:   [--- TE end (50 bp) ---][--- ref right flank (50 bp) ---]
  Position:    0                     50                              99
  Case:        UPPERCASE              lowercase

LEFT junction (type=left):
  Sequence:   [--- ref left flank (50 bp) ---][--- TE start (50 bp) ---]
  Position:    0                             50                        99
  Case:        lowercase                      UPPERCASE
```

## Building the Competitive Reference

After Phase 1, build `junctions.fa` from all discovered junctions:

```bash
python scripts/build_junctions_ref.py temp_work junctions.fa
```

This scans all `junction_*.fasta` files, extracts Abs (absence/WT_REF) and Pre (presence/junction) sequences, deduplicates by genomic position + TE name, and writes:

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
    ├── find_te_junctions.py             # Phase 1 graph + k-mer walk discovery
    ├── build_junctions_ref.py           # Bridge: build competitive reference
    ├── submit_te_kmer_count.sh          # Phase 2 SLURM orchestrator
    ├── run_te_kmer_count.sh             # Phase 2 BBDuk per-sample
    ├── extract_junction_kmers.py        # Phase 2 k-mer extraction
    ├── genotype_from_counts.py          # Phase 2 genotyping
    └── archive/
        ├── build_junctions_from_reads.py  # Previous read-level approach
        └── run_te_assembly_spades.sh      # Previous SPAdes-only approach
```

## Troubleshooting

**No junction files produced (Phase 1):**
- Check `gold_standard.tsv` — if empty, no candidate reads with FBte* mates were found
- Inspect `assembly_graph.fastg` — if empty/missing, SPAdes assembly failed
- Check `nodes_vs_te.tsv` — if empty, no graph nodes matched TEs
- If TE-boundary nodes exist but walks fail, try lower `--min-vote-frac`
- Try per-k-value graphs: `K21/`, `K33/`, `K55/` may have junctions the main graph lost

**SPAdes fails or produces empty graph:**
- Check `spades.log` for errors
- Ensure R1.fq and R2.fq are properly paired
- Try with fewer reads or a smaller k range

**No k-mer matches (Phase 2):**
- Check that FASTQ files are correct for the samples
- Verify `junctions.fa` contains the expected Abs/Pre pairs
- Inspect `abs_kmers.fa` / `pre_kmers.fa` to confirm diagnostic k-mers were extracted
- Check BBDuk stats files for total read counts

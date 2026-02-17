# TE Junction Assembly & Genotyping Pipeline

Two-phase pipeline: discover TE insertions from BAMs, then genotype them across samples using k-mer counting.

## Phases

**Phase 1 (Discovery):** Pool reads from BAMs per region, assemble with SPAdes, identify TE junction contigs (BLAST + minimap2), generate 100bp visual alignments, build `junctions.fa`.

**Phase 2 (Genotyping):** Extract diagnostic k-mers from junction sequences, scan raw FASTQs with BBDuk, call genotypes from k-mer counts.

## Requirements

HPC modules:
- `samtools/1.15.1`
- `SPAdes/3.15.4`
- `ncbi-blast/2.13.0`
- `minimap2/2.28`
- `python/3.10.2` with Biopython

Conda:
- `bbmap` (BBDuk for k-mer counting) -- activate via `mamba activate bbmap`

## Configuration

All settings in **`config.sh`** (the only file you edit).

## Quick Start

```bash
cd /dfs7/adl/tdlong/Sarah/local_TE

# 1. Edit config.sh with BAMs, regions, and FASTQ paths

# 2. Phase 1: Discover TE insertions
sbatch scripts/submit_te_analysis.sh

# 3. Phase 2: Genotype via k-mer counting
sbatch scripts/submit_te_kmer_count.sh
```

## Files

```
local_TE/
├── config.sh                            # BAMs, regions, FASTQs, refs
├── junctions.fa                         # Built by Phase 1
├── junctions_metadata.tsv
└── scripts/
    ├── submit_te_analysis.sh            # Phase 1 SLURM wrapper
    ├── run_te_assembly.sh               # Phase 1 pipeline
    ├── build_te_alignment.py            # Phase 1 visualization
    ├── build_junctions_ref.py           # Build competitive reference
    ├── submit_te_kmer_count.sh          # Phase 2 SLURM orchestrator
    ├── run_te_kmer_count.sh             # Phase 2 BBDuk per-sample
    ├── extract_junction_kmers.py        # Phase 2 k-mer extraction
    └── genotype_from_counts.py          # Phase 2 genotyping
```

## Cluster

- Account: `tdlong_lab`
- See `slurm.md` for partition details and templates

## Git

Commit scripts and docs. Never commit data files (>10MB), results, or logs.

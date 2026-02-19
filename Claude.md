# TE Junction Assembly & Genotyping Pipeline

Two-phase pipeline: discover TE insertions from BAMs, then genotype them across samples using k-mer counting.

## Phases

**Phase 1 (Discovery):** Pool reads from BAMs per region, assemble with SPAdes (metagenome mode), parse the assembly graph for TE-containing nodes, k-mer walk from TE boundaries into reference sequence, cross-validate against a gold standard catalog built from BAM mate-pair info, and write `junction_*.fasta` files → `junctions.fa`.

**Phase 2 (Genotyping):** Extract diagnostic k-mers from junction sequences, scan raw FASTQs with BBDuk, call genotypes from k-mer counts.

## Important: Data files live on HPC

Code is written locally and pushed via git. All data files (BAMs, FASTQs, results, assembly output) live on the HPC cluster at `/dfs7/adl/tdlong/Sarah/local_TE/`. To examine output files, check on the cluster.

## Requirements

HPC modules:
- `samtools/1.15.1`
- `ncbi-blast/2.13.0`
- `SPAdes/3.15.4`
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
    ├── run_te_assembly.sh               # Phase 1 pipeline (SPAdes + k-mer walk)
    ├── find_te_junctions.py             # Phase 1 graph + k-mer walk junction discovery
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

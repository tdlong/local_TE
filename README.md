# TE Junction Assembly & Genotyping Pipeline

A two-phase pipeline for detecting transposable element (TE) insertions from paired-end sequencing data and genotyping them across samples.

## Overview

**Phase 1 (Discovery):** Pool reads from all BAMs for each genomic region, assemble with SPAdes (metagenome mode), parse the assembly graph for TE-containing nodes, and write junction files. Three approaches are used to extract junctions:
1. **Junction nodes** — graph nodes that already contain both TE and reference sequence
2. **K-mer walking** — extend TE nodes into reference using a read k-mer index
3. **Graph edge tracing** — follow FASTG connectivity from TE nodes to reference neighbors

Results are cross-validated against a gold standard catalog built from BAM mate-pair info.

**Phase 2 (Genotyping):** Extract diagnostic k-mers from junction sequences, scan raw FASTQs with BBDuk, and call genotypes from k-mer counts.

## Requirements

HPC modules:
- `samtools/1.15.1`
- `ncbi-blast/2.13.0`
- `SPAdes/3.15.4`
- `python/3.10.2` with Biopython

Conda:
- `bbmap` (BBDuk for k-mer counting): `mamba activate bbmap`

## Configuration

All user-editable settings live in **`config.sh`** in the project root.

## Quick Start

```bash
cd /dfs7/adl/tdlong/Sarah/local_TE

# 1. Edit config.sh with your BAMs, regions, and FASTQ paths

# 2. Phase 1: Discover TE insertions → junctions.fa
sbatch scripts/submit_te_analysis.sh
# → junctions.fa + junctions_metadata.tsv

# 3. Review junctions.fa, then Phase 2: Genotype via k-mer counting
cat junctions.fa
sbatch scripts/submit_te_kmer_count.sh
# → kmer_work/genotype_results.tsv        (one row per sample × TE junction)
```

## Phase 1: Discovery

### Pipeline Sections

1. **Section 1 — Extract Candidate Reads:** Pools junction-candidate reads (mate not on same chromosome) from all BAMs. Outputs `candidates_catalog.tsv` preserving sample identity and mate contig.

2. **Section 2 — Paired Reads + Gold Standard:** Extracts paired-end FASTQs (R1.fq, R2.fq) and reference region. Builds a **gold standard catalog** (`gold_standard.tsv`) from reads whose mates map to FBte* contigs — this tells us before any assembly which TEs are expected, approximate positions, and per-sample read counts.

3. **Section 3 — SPAdes Assembly:** Assembles reads with `spades.py --meta -k 21,33,55`. We use the **assembly graph** (`assembly_graph.fastg`), not `contigs.fasta`, because the graph preserves all edges including low-coverage TE junction paths.

4. **Section 4 — Find TE Junctions** (`find_te_junctions.py`):
   - Parse graph nodes and edge connectivity from FASTG
   - BLAST graph nodes vs TE database and reference region
   - Classify nodes: junction (both hits), TE-boundary, TE-interior, ref-only
   - **Junction nodes**: extract 100bp junction directly from nodes with both TE and ref hits
   - **K-mer walk**: extend TE boundary/interior nodes into reference via read k-mer index
   - **Graph edges**: concatenate TE nodes with adjacent reference neighbors from FASTG connectivity
   - Cross-validate discovered junctions against gold standard catalog
   - Write `junction_{type}_{N}.fasta` per junction

### Gold Standard Catalog

The gold standard (`gold_standard.tsv`) is built automatically from BAM mate-pair info. It summarizes reads whose mates map to FBte* contigs:

```
te_name        approx_pos  total_reads  HOULE_L1F  HOULE_L2F  HOUSTON_L1F  ...
FBte0000626    8711446     14           0          14          0
FBte0001399    8738348     2            1          1           0
```

This provides ground truth for validating junction discovery. TEs with only 1-2 supporting reads are often not assemblable. The gold standard only counts FBte* mates; TEs supported by discordant mates to other chromosomes may also be discovered as "SUSPECT" junctions.

### Interactive Testing

```bash
srun --pty bash
bash scripts/test_phase1.sh          # runs first region from config.sh
bash scripts/summarize_phase1.sh     # dumps diagnostics to phase1_summary.txt
```

### Junction File Format

Each `junction_{type}_{N}.fasta` contains 4 records, each ~100 bp, with the
insertion point at position 50 (0-indexed).

**Case convention:** reference-derived bases are written **lowercase**; TE-derived
bases are written **UPPERCASE**. The case boundary marks the insertion point.

#### The 4 Records

1. **WT_REF** — 100bp reference (absence allele), all lowercase
2. **REF** — duplicate of WT_REF
3. **walked_junction** — 100bp junction (presence allele), case-encoded
4. **TE canonical** — 100bp of canonical TE sequence for comparison

#### LEFT vs. RIGHT Junction Types

```
RIGHT junction (type=right):
  Sequence:   [--- TE end (50 bp) ---][--- ref right flank (50 bp) ---]
  Case:        UPPERCASE              lowercase

LEFT junction (type=left):
  Sequence:   [--- ref left flank (50 bp) ---][--- TE start (50 bp) ---]
  Case:        lowercase                      UPPERCASE
```

Not all TEs produce both left and right junctions. A single junction is sufficient for genotyping. When both are found, their positions should differ by ~4-8bp (target site duplication).

## Building the Competitive Reference

After Phase 1, build `junctions.fa` from all discovered junctions:

```bash
python scripts/build_junctions_ref.py temp_work junctions.fa
```

This scans all `junction_*.fasta` files across regions, extracts Abs/Pre pairs, deduplicates, and writes `junctions.fa` + `junctions_metadata.tsv`.

## Phase 2: Genotyping via K-mer Counting

### How It Works

1. **Extract diagnostic k-mers** (`extract_junction_kmers.py`):
   - From each Abs/Pre pair, extract k=31 k-mers spanning the junction at position 50
   - Keep only allele-specific k-mers (set difference: abs - pre, pre - abs)

2. **Count k-mers per sample** (`run_te_kmer_count.sh`):
   - BBDuk scans each sample's paired FASTQs against the k-mer databases
   - Flags: `k=31 maskmiddle=f rcomp=t`

3. **Call genotypes** (`genotype_from_counts.py`):
   - Compute abs_ratio = abs_count / (abs_count + pre_count)
   - Call: Abs/Abs (ratio >= 0.8), Pre/Pre (ratio <= 0.2), Abs/Pre (between), or NO_CALL (depth < 5)

### Running

```bash
sbatch scripts/submit_te_kmer_count.sh
```

## Files

```
local_TE/
├── README.md
├── config.sh                            # Edit this: BAMs, regions, FASTQs, refs
├── junctions.fa                         # Built by build_junctions_ref.py
├── junctions_metadata.tsv
└── scripts/
    ├── submit_te_analysis.sh            # Phase 1 SLURM wrapper
    ├── run_te_assembly.sh               # Phase 1 pipeline (SPAdes + junction finding)
    ├── find_te_junctions.py             # Phase 1 junction discovery (graph + walk + edges)
    ├── build_junctions_ref.py           # Build competitive reference from junctions
    ├── test_phase1.sh                   # Interactive Phase 1 test
    ├── summarize_phase1.sh              # Diagnostic summary dump
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
- TEs with only 1-2 gold standard reads are often not assemblable

**SPAdes fails or produces empty graph:**
- Check `spades.log` for errors
- Ensure R1.fq and R2.fq are properly paired

**No k-mer matches (Phase 2):**
- Verify `junctions.fa` contains the expected Abs/Pre pairs
- Inspect `abs_kmers.fa` / `pre_kmers.fa` to confirm diagnostic k-mers were extracted
- Check BBDuk stats files for total read counts

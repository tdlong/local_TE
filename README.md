# TE Junction Assembly & Genotyping Pipeline

A two-phase pipeline for detecting transposable element (TE) insertions from paired-end sequencing data and genotyping them across samples.

## Overview

**Phase 1 (Discovery):** Pool reads from all BAMs for each genomic region, assemble with SPAdes (metagenome mode), parse the assembly graph for TE-containing nodes, extract junction sequences, validate them, and build `junctions.fa`.

**Phase 2 (Genotyping):** Extract diagnostic k-mers from junction sequences, scan raw FASTQs with BBDuk, and call genotypes from k-mer counts.

## Test Run Results

End-to-end test on 2 regions (chr3L:8710861-8744900, chrX:8450000-8500000) across 6 samples (3 HOULE, 3 HOUSTON):

```
junction                              HOULE_L1  HOULE_L2F  HOULE_L3F  HOUSTON_L1F  HOUSTON_L2F  HOUSTON_L3F
chr3L_8711446_FBte0000626_right       0.00      0.05       0.00       0.00         0.00         0.00
chr3L_8734174_FBte0000559_left        NA        0.00       0.00       0.43         0.00         0.00
chrX_8490354_FBte0000021_right        0.60      0.71       0.80       0.02         0.03         0.08
```

3 TE insertions discovered and genotyped, all showing population-specific frequency patterns:

- **FBte0000626** (chr3L:8711446): Rare insertion at ~5% in HOULE_L2F only. Consistent with the gold standard (14 reads, all from HOULE_L2F).
- **FBte0000559** (chr3L:8734174): Present at ~43% in HOUSTON_L1F, absent elsewhere. Only the left junction was recovered; the right junction (at 8734167, 7bp TSD) was found in earlier test runs but not in the production assembly — likely due to SPAdes graph variability between runs.
- **FBte0000021** (chrX:8490354): Strong population differentiation — 60-80% in HOULE, 2-8% in HOUSTON.

### Known Limitations

- **Not all junctions are recovered.** SPAdes assembly is stochastic; low-coverage junctions may appear in some runs but not others. K-mer walks rarely succeed because SPAdes error-corrects reads, breaking the k-mer chain between assembled nodes and the original read index.
- **Single junctions are sufficient.** Not all TEs produce both left and right junctions. A single junction provides enough diagnostic k-mers for genotyping.
- **NA = insufficient depth.** Phase 2 reports NA when fewer than 5 total k-mers are found for a junction in a sample.

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

# 3. Review junctions.fa before proceeding
cat junctions.fa

# 4. Phase 2: Genotype via k-mer counting (~1 hour)
sbatch scripts/submit_te_kmer_count.sh
# → kmer_work/genotype_results_freq_matrix.tsv
```

## Phase 1: Discovery

### How It Works

1. **Extract candidate reads** from all BAMs — reads whose mates map to a different chromosome. Outputs `candidates_catalog.tsv` preserving sample identity and mate contig.

2. **Build gold standard catalog** from reads whose mates map to FBte* contigs. This tells us before any assembly which TEs to expect, their approximate positions, and which samples contain them.

3. **SPAdes assembly** (`--meta -k 21,33,55`). We use the **assembly graph** (`assembly_graph.fastg`), not `contigs.fasta`, because the graph preserves all paths including low-coverage TE junctions.

4. **Find TE junctions** (`find_te_junctions.py`):
   - BLAST graph nodes vs TE database and reference region
   - Classify nodes: junction (both hits), TE-boundary, TE-interior, ref-only
   - Extract junctions from junction nodes; attempt k-mer walks and graph edge tracing for the rest
   - Validate each junction (ref half must match Abs) and deduplicate by position
   - Write `junction_{type}_{N}.fasta` per junction

5. **Build competitive reference** (`build_junctions_ref.py`): collect junction files from the directories just processed (not a wildcard), validate, deduplicate by position, write `junctions.fa`.

### Junction Discovery

**Junction nodes** are the primary source. These are graph nodes that BLAST to both a TE and the reference region. The junction point is determined from the **reference BLAST boundary** — where the reference alignment begins (right junction) or ends+1 (left junction) in the node. This ensures the junction extraction point and the genomic position are derived from the same source.

**K-mer walking** attempts to extend TE-only nodes into reference using a read k-mer index. In practice, SPAdes error-correction often prevents walks from succeeding because the corrected node sequences don't match the original read k-mers.

**Graph edge tracing** follows FASTG connectivity from TE nodes to adjacent reference nodes. Useful when the junction spans two nodes connected by a graph edge.

### Built-in Correctness Check

Each junction's reference half (the lowercase portion) must match the corresponding half of the Abs (absence allele) — both are the same stretch of reference genome. If they don't match (>2 mismatches out of 50bp), the junction point is wrong and the junction is **rejected**. This is checked in both `find_te_junctions.py` and `build_junctions_ref.py`.

This validation caught two bugs during development: an inconsistency between TE-BLAST-derived and ref-BLAST-derived junction points, and an off-by-one in left junction coordinate calculation.

### Deduplication

Related TE families can share sequence, causing a single physical insertion to BLAST to multiple TE database entries. Junctions are deduplicated by **(position, side)** regardless of TE name — same position means same insertion. The best BLAST hit is kept. This is enforced in both `find_te_junctions.py` (during discovery) and `build_junctions_ref.py` (during collection).

### Gold Standard Catalog

The gold standard (`gold_standard.tsv`) is built from BAM mate-pair info:

```
te_name        approx_pos  total_reads  HOULE_L1F  HOULE_L2F  HOUSTON_L1F  ...
FBte0000626    8711446     14           0          14          0
FBte0001399    8738348     2            1          1           0
```

This provides ground truth for validating junction discovery. Notes:
- TEs with only 1-2 supporting reads are often not assemblable
- The gold standard only counts FBte* mates; other discordant mates may produce "SUSPECT" junctions
- TE insertions can be truncated/incomplete — a mate mapping deep into a canonical TE doesn't mean the full TE is inserted

### Junction File Format

Each `junction_{type}_{N}.fasta` contains 4 records, each 100bp, with the insertion point at position 50 (0-indexed).

**Case convention:** reference bases are **lowercase**; TE bases are **UPPERCASE**. The case boundary marks the insertion point.

```
RIGHT junction:  [UPPERCASE TE end 50bp][lowercase ref right flank 50bp]
LEFT  junction:  [lowercase ref left flank 50bp][UPPERCASE TE start 50bp]
```

Not all TEs produce both left and right junctions. A single junction is sufficient for genotyping. When both are found, their positions should differ by ~4-8bp (target site duplication).

### Interactive Testing

```bash
srun --pty bash
bash scripts/test_phase1.sh          # runs first region from config.sh
bash scripts/summarize_phase1.sh     # dumps diagnostics to phase1_summary.txt
```

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
├── config.sh                            # Edit this: BAMs, regions, FASTQs, refs
├── junctions.fa                         # Built by Phase 1
├── junctions_metadata.tsv
└── scripts/
    ├── submit_te_analysis.sh            # Phase 1 SLURM wrapper (includes build_junctions)
    ├── run_te_assembly.sh               # Phase 1 per-region pipeline
    ├── find_te_junctions.py             # Junction discovery (graph + walk + edges)
    ├── build_junctions_ref.py           # Collect, validate, and deduplicate junctions
    ├── test_phase1.sh                   # Interactive Phase 1 test
    ├── summarize_phase1.sh              # Diagnostic summary
    ├── submit_te_kmer_count.sh          # Phase 2 SLURM orchestrator
    ├── run_te_kmer_count.sh             # Phase 2 BBDuk per-sample
    ├── extract_junction_kmers.py        # Phase 2 k-mer extraction
    └── genotype_from_counts.py          # Phase 2 genotyping
```

## Troubleshooting

**No junction files produced:**
- Check `gold_standard.tsv` — if empty, no candidate reads with FBte* mates were found
- Inspect `assembly_graph.fastg` — if empty/missing, SPAdes assembly failed
- Check `nodes_vs_te.tsv` — if empty, no graph nodes matched TEs
- TEs with only 1-2 gold standard reads are often not assemblable

**Junctions rejected by correctness check:**
- The ref half of the Pre didn't match the Abs — junction point determination is wrong
- Check the BLAST hits in `nodes_vs_te.tsv` and `nodes_vs_ref.tsv` for the offending node
- Common causes: TE and ref BLAST hits don't cleanly abut in the node, strand handling issues

**No k-mer matches (Phase 2):**
- Verify `junctions.fa` contains expected Abs/Pre pairs
- Inspect `abs_kmers.fa` / `pre_kmers.fa` to confirm diagnostic k-mers were extracted
- Check BBDuk stats files for total read counts

## Future Work

- **Scale to more regions.** The test run used 2 regions; the pipeline is ready for genome-wide application by adding regions to `config.sh`.
- **Recover more junctions.** The FBte0000559 right junction and FBte0000626 left junction exist but weren't recovered in the production assembly. Possible improvements: use per-k-value SPAdes graphs (K21, K33) in addition to the final merged graph, or use SPAdes error-corrected reads for the k-mer walk index.
- **Handle truncated TE insertions.** Some TE insertions are incomplete — mate pairs mapping deep into a canonical TE don't mean the full TE is present. The gold standard catalog flags these but the pipeline doesn't yet distinguish them.

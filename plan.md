# Plan: Replace SPAdes Assembly with Read-Level Junction Discovery

## What Changed

Phase 1 previously used SPAdes to assemble junction-spanning contigs, then
minimap2 to align them. SPAdes fails systematically when coverage is low,
reads are diverged from the reference, or the TE region is complex — but
individual junction reads are still present. The new approach BLASTs
individual reads directly against the TE database and the reference region,
identifies junction reads (hits to both), clusters by inferred insertion
position, and builds a consensus from those reads. No assembly, no Java, no
graph traversal.

## Why

- Works at any coverage level (even 2-3 reads spanning a junction)
- No SPAdes graph-collapse failures for repetitive or diverged TEs
- Consensus captures TSDs naturally (they appear in actual reads)
- IUPAC SNP encoding in Pre consensus captures all TE-bearing haplotypes
  in a pooled population without extra Phase 2 complexity

## New Script: `scripts/build_junctions_from_reads.py`

Replaces Sections 3–7 of `run_te_assembly.sh`. Called after read extraction.

**Inputs:**
- `outdir/R1.fq`, `outdir/R2.fq`, `outdir/singles.fq` — junction reads
- `outdir/region.fasta` — reference region sequence
- `$TEFASTA` — TE sequence database
- `$REGION` — region string (e.g. `chr3L:8710861-8744900`)
- `--min-support` (default 3) — minimum reads per cluster
- `--snp-min-freq` (default 0.15) — minor allele threshold for IUPAC encoding

**Algorithm:**
1. Convert FASTQ → FASTA
2. BLAST reads vs TE database (`blastn -subject $TEFASTA`)
3. BLAST reads vs reference region (`blastn -subject region.fasta`)
4. Junction reads = reads with hits to both; infer insertion position
5. Cluster reads within ±20bp of the same inferred insertion point
6. For each cluster: anchor reads at TE boundary (LEFT) or ref boundary (RIGHT),
   build base-count consensus, encode IUPAC at SNP positions ≥ `--snp-min-freq`
7. Abs sequence = reference genome slice centred at insertion
8. Write `junction_{type}_{N}.fasta` (4-record format for `build_junctions_ref.py`)

**Output format** (unchanged; compatible with downstream scripts):
```
>WT_REF[start:end] insertion=chr3L:8711446 te=FBte0000626 type=left
[100bp reference slice — Abs allele]
>REF
[same 100bp reference slice]
>reads_consensus_left_N
[100bp read consensus with IUPAC at SNP positions — Pre allele]
>FBte0000626
[100bp of TE sequence from database]
```

## Design Decisions

- **Anchor at `te_qstart`** (LEFT) or `ref_qstart` (RIGHT): trusts the TE
  boundary, which is more conserved than the reference flank.
- **Abs from reference genome**: diagnostic k-mers inherently contain TE
  sequence, so they cannot match non-TE reads regardless of SNP state in the
  flank. No IUPAC needed for Abs.
- **IUPAC in Pre consensus**: `extract_junction_kmers.py` expands these to
  all concrete k-mer variants before writing `pre_kmers.fa` for BBDuk.

## Archived Scripts

- `scripts/archive/run_te_assembly_spades.sh` — original SPAdes-based pipeline
- `scripts/archive/reconstruct_junctions.sh` — prototype that motivated this change

## Modules Removed

- `SPAdes/3.15.4` — no longer needed
- `minimap2/2.28` — no longer needed

## Verification

```bash
cd /dfs7/adl/tdlong/Sarah/local_TE
git pull

srun --pty bash
bash scripts/run_te_assembly.sh \
    "chr3L:8710861-8744900" \
    "temp_work/chr3L_test" \
    "/path/to/sample.bam"

# Expect:
#   junction_*.fasta files with correct 4-record format
#   WT_REF headers contain insertion= and te= fields
#   Pre sequences contain IUPAC codes at SNP positions
#   Junction count ≥ what SPAdes found on the same region

python scripts/build_junctions_ref.py temp_work junctions.fa
python scripts/extract_junction_kmers.py junctions.fa kmer_work
```

# Genomic Region Annotator — Design Report (BioGeMT)

This report documents the major design decisions and the resulting pipeline.

## Background / motivation

The original motivation was CLASH-like datasets where an input “feature” column (e.g. exon/intron/UTR) can be misleading because it is not transcript-aware.
A genomic interval can overlap multiple transcripts with different region interpretations (CDS vs 3′UTR, etc.).

Therefore, a single-pass “final region” annotator would hide ambiguity and could produce systematically biased labels.

## Core decision: two-step workflow

### Step 1 — evidence only (no assumptions)

**Goal:** Preserve evidence and ambiguity.

Given an interval (chr, start, end, strand):

1. Find all transcripts that overlap the interval sufficiently.
2. For each passing transcript, compute per-nucleotide membership in:
   - CDS, UTR5, UTR3, EXON, INTRON (derived), and INTERGENIC.
3. Union (OR) these per-nt tracks across all passing transcripts.
4. Output:
   - `_input_with_ids.tsv` with stable ids (001…)
   - `_matrix.tsv`: one row per region (CDS/UTR5/UTR3/EXON/INTRON/INTERGENIC) and columns nt_1..nt_L
   - `_transcripts.tsv`: one row per passing transcript per read with overlap bp summaries
   - `_step1_stats.tsv`: dataset-level evidence summary

**Overlap filter:**
- Default is **100% containment**: read must lie inside transcript span.
- Optional `--min-overlap-nt N` keeps transcripts with overlap >= N nt.

This filter is intentionally exposed because different applications need different sensitivity.

### Why keep INTERGENIC in Step 1?

Step 1 models “evidence” at the per-nt level. Even if reads are fully contained within transcripts, INTERGENIC can still be non-zero when:
- no transcript passes the overlap filter for that read, or
- strand / chromosome constraints remove candidates.

Therefore, INTERGENIC is useful to distinguish “no transcript evidence” cases.

### Step 2 — explicit assumptions (transcript selection)

**Goal:** Pick a transcript with an explicit policy, then summarize region composition.

Inputs:
- Step1 `_transcripts.tsv` and `_matrix.tsv`

Outputs:
- `_site_summary.tsv` with:
  - selected transcript id/gene
  - dominant region (selected transcript)
  - multi-region composition (selected transcript)
  - dominant region (UNION evidence)
  - multi-region composition (UNION evidence)
  - an “ambiguous” flag if UNION dominant differs from selected dominant
- `_step2_stats.tsv` with:
  - dominant region distributions
  - top multi-region patterns (e.g. CDS|UTR3)
  - ambiguity rate
  - top genes/transcripts selected

## Dominant region vs multi-region evidence

We explicitly output both because they answer different questions:

- **dominant_region**: “If forced to pick one label, what is the best summary?”
- **regions_present**: “Which regions are involved at all?” (lossless)

This is important for boundary-spanning sites (e.g. CDS|UTR3).

## Dominance modes

### coverage-dominant (recommended)
Pick the region with the most bp overlap.

This prevents a tiny UTR sliver (1–2 nt) from overriding a mostly-CDS site.

### priority-dominant
Pick the first region present by a priority list.

This is useful for some UTR-centric reports, but is intentionally not the default.

## Current transcript selection policy for CLASH

Implemented policy: `clash_utr3_first`

Lexicographic scoring:

1) max UTR3 bp  
2) max CDS bp  
3) max UTR5 bp  
4) max EXON_OTHER bp  
5) max INTRON bp  
6) max EXON bp  
7) max TX overlap bp  
Tie-breakers: contained_100pct, transcript_id

This policy matches typical miRNA assumptions without hard-coding “UTR3 always wins”.

## Output folder conventions

To avoid clutter and keep provenance clear:

- Step 1 writes to: `data/processed/step1/`
- Step 2 writes to: `data/processed/step2/`

If inputs are not in `step1/`, Step 2 writes next to the transcripts file.

## Known limitations / next steps

- Add additional transcript policies (canonical/MANE/protein_coding first, etc.)
- Add optional gene-level aggregation (pick one transcript per gene, etc.)
- Consider parallelization for very large datasets (currently optimized via filtered GTF + caching + binning)

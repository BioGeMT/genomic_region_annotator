# Genomic Region Annotator 

A genomic interval annotation tool built on Ensembl GTF files, designed for **read/peak/interval** datasets (e.g. CLASH, eCLIP, ChIP peaks) where transcript ambiguity is common.

This repository implements a **2-step workflow**:

- **Step 1 (evidence-only)**: compute *all* transcript overlaps and build a per-nt **UNION region matrix** across transcripts (no transcript-selection assumptions).
- **Step 2 (assumption-bearing)**: choose a transcript using an explicit policy (e.g. CLASH-friendly) and summarize region composition (dominant + multi-region).

---

## Why 2 steps

- Step 1 preserves evidence: *what regions are supported by any transcript?*
- Step 2 makes assumptions explicit: *which transcript policy are you applying and why?*

---

## Installation

```bash
conda env create -f environment.yml
conda activate genomic-region-annotator
pip install -e .
```

---

## Input format

Input must be a TSV with:

| column | description |
|---|---|
| chr | Chromosome (e.g. `1`, `chr1`) |
| start | Start coordinate |
| end | End coordinate |
| strand | `+` / `-` (or `.` if unknown) |

---

## Coordinate conventions

Use `--coords`:

- `1-based` (default): 1-based inclusive (GTF/SAM-like)
- `bed`: 0-based half-open (BED)

Internally, normalization to **1-based inclusive** for overlap math.

---

## Step 1 — annotate (evidence-only)

Annotate an input TSV and produce step1 outputs:

```bash
genomic-region-annotator annotate \
  --input data/raw/file_with_intervals.tsv \
  --release <version_of_ensembl> \
  --output data/processed/output.tsv
```

### Outputs (written into `data/processed/step1/`)

The output stem is taken from `--output`. If your output path is `data/processed/output.tsv`, the stem is `output` and you get:

- `data/processed/step1/output_input_with_ids.tsv`
- `data/processed/step1/output_matrix.tsv`
- `data/processed/step1/output_transcripts.tsv`
- `data/processed/step1/output_step1_stats.tsv`

### IDs and row tracking

Step 1 adds a stable zero-padded `id` column to every input row and writes the full input table to `*_input_with_ids.tsv`. This `id` is the join key used by all downstream files:

- `*_input_with_ids.tsv` has one row per original input interval.
- `*_matrix.tsv` has per-region, per-nucleotide union evidence keyed by `id`.
- `*_transcripts.tsv` has one row per passing transcript keyed by `id`; an input interval with no passing transcript will not have transcript rows.

If the input already contains an `id` column, it is preserved as `original_id`, and the pipeline-generated `id` is still used internally for stable row tracking.

### Transcript overlap filter

By default, Step 1 keeps only transcripts where the read is **100% contained** inside the transcript span.

You can relax this with:

```bash
--min-overlap-nt 30
```

meaning: keep transcripts with **≥ 30 nt overlap** with transcript span.

---

## Step 2 — summarize-sites (explicit transcript choice + region summary)

Run transcript selection + site summary:

```bash
genomic-region-annotator summarize-sites \
  --transcripts data/processed/step1/output_transcripts.tsv \
  --matrix data/processed/step1/output_matrix.tsv \
  --input-with-ids data/processed/step1/output_input_with_ids.tsv \
  --policy clash_utr3_first \
  --dominance coverage \
  --report
```

### Outputs (written into `data/processed/step2/`)

With the default output paths, Step 2 writes:

- `data/processed/step2/output_site_summary.tsv`
- `data/processed/step2/output_final.tsv`
- `data/processed/step2/output_step2_stats.tsv`

The compact `*_final.tsv` file contains the most important columns for downstream review; `*_site_summary.tsv` keeps the fuller site-level output.

### Selection status and retained rows

When `--input-with-ids` is provided, Step 2 writes one output row for every input interval from `*_input_with_ids.tsv`. The `selection_status` column explains whether transcript selection was possible:

- `SELECTED_TRANSCRIPT`: at least one transcript passed Step 1, and Step 2 selected one transcript according to `--policy`.
- `NO_TRANSCRIPT`: no transcript passed the Step 1 transcript-overlap filter for that input interval. These rows are retained so the final output keeps the same row universe as the input.

`NO_TRANSCRIPT` is different from `INTERGENIC`:

- `NO_TRANSCRIPT` is a transcript-selection status: there was no passing transcript to select.
- `INTERGENIC` is a region label from the union/evidence matrix.

For `NO_TRANSCRIPT` rows, selected-transcript fields are blank, `dominant_region_selected` and `regions_present_selected` are set to `NO_TRANSCRIPT`, and the union fields (`dominant_region_union`, `regions_present_union`, and `bp_*_union`) still summarize the Step 1 matrix evidence for that input interval.

---

## Region concepts (Step 2)

Step 2 outputs both:

### 1) Dominant region (single label)
- `dominant_region_selected`
- `dominant_region_union`

**Dominance modes:**
- `--dominance coverage` (recommended): label = region with most bp overlap (ties broken by priority)
- `--dominance priority`: label = first region present by priority order (UTR3 > CDS > ...)

### 2) Multi-region evidence (lossless)
- `regions_present_selected`
- `regions_present_union`

Example:
- `CDS|UTR3` means the site spans both CDS and 3′UTR.

---

## CLASH policy currently supported

`--policy clash_utr3_first` selects a transcript per read by:

1. maximize UTR3 overlap bp  
2. then CDS overlap bp  
3. then UTR5 overlap bp  
4. then EXON_OTHER bp (exon excluding CDS/UTRs)  
5. then INTRON bp  
6. then EXON bp  
7. then TX overlap bp  
8. tie-breakers: contained_100pct, then transcript_id

---


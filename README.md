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
  --output data/processed/output_name.tsv
```

### Outputs (written into `data/processed/step1/`)

If your output stem is `output`, you get:

- `data/processed/step1/output_annotated_input_with_ids.tsv`
- `data/processed/step1/output_annotated_matrix.tsv`
- `data/processed/step1/output_annotated_transcripts.tsv`
- `data/processed/step1/output_annotated_step1_stats.tsv`

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
  --matrix data/processed/step1/output_annotated_matrix.tsv \
  --policy clash_utr3_first \
  --dominance coverage \
  --report
```

### Outputs (written into `data/processed/step2/`)

- `data/processed/step2/output_site_summary.tsv`
- `data/processed/step2/output_step2_stats.tsv`

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


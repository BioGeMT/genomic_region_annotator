# Genomic Region Annotator

A **transcript-first genomic interval annotation tool** based on Ensembl GTF files.

> **Key principle**  
> This is a transcript-first region annotator: each interval is assigned to a single transcript, and all region labels and coordinates are interpreted in that transcript’s context.

---

## Overview

`genomic-region-annotator` annotates genomic intervals (e.g. CLASH reads, peaks, binding sites) with **transcript-aware functional regions**, such as:

- five_prime_utr
- three_prime_utr
- CDS
- start_codon / stop_codon
- exon
- intron (derived)
- transcript / gene
- intergenic

The tool is **not CLASH-specific** — it works for any genomic intervals.

### What makes it different?

- **Transcript-first annotation**
  - Each input interval is assigned to **one transcript**, not multiple genes/transcripts.
  - Transcript selection prioritizes **protein-coding transcripts** and **maximum overlap**.
- **Multi-label annotation**
  - You get:
    - `region_all` – all overlapping regions
    - `region_dominant` – dominant region by overlap
    - `region_primary` – single region by biological priority
- **Relative coordinates**
  - Outputs **strand-aware transcript-relative coordinates** (`rel_start`, `rel_end`)
- **Intron inference**
  - Introns are derived when not explicitly present in the GTF
- **Fast**
  - GTFs are cached as Parquet on first use

---

## Installation

### 1. Clone the repository

```bash
git clone git@github.com:BioGeMT/genomic_region_annotator.git
cd genomic_region_annotator
2. Create the conda environment
bash
Copy code
conda env create -f environment.yml
conda activate genomic-region-annotator
3. Install the package (editable mode)
bash
Copy code
pip install -e .
Input format
Input must be a TSV file with the following columns:

column	description
chr	Chromosome (e.g. 1, chr1, MT)
start	Start coordinate
end	End coordinate
strand	+ or -

Coordinate conventions
Use --coords to specify input coordinates:

1-based (default): 1-based inclusive (GTF / GFF / SAM-like)

bed: 0-based half-open (BED format)

Internally, all coordinates are converted to 0-based half-open for correctness and consistency.

Basic usage
Annotate using an Ensembl release (auto-download)
bash
Copy code
genomic-region-annotator annotate \
  --input data/raw/intervals.tsv \
  --release 115 \
  --output annotated.tsv
Annotate using a local GTF
bash
Copy code
genomic-region-annotator annotate \
  --input intervals.tsv \
  --gtf Homo_sapiens.GRCh38.115.gtf \
  --output annotated.tsv
Output columns
column	description
region_all	All overlapping region labels (comma-separated)
region_dominant	Region with maximum bp overlap
region_primary	Single region chosen by priority
transcript_id	Selected transcript
transcript_name	Transcript name
transcript_biotype	Transcript biotype
gene_id	Gene ID
gene_name	Gene name
gene_biotype	Gene biotype
rel_start	Transcript-relative start (1-based, strand-aware)
rel_end	Transcript-relative end (1-based, strand-aware)
annotation_notes	Warnings or fallback behavior notes

Region logic
region_all
All regions overlapped by the interval within the selected transcript.

Example:

Copy code
three_prime_utr,CDS,stop_codon,exon
region_dominant
The region with the largest base-pair overlap.

Tie-breaker: biological priority.

region_primary
A single biologically prioritized label, using this order:

five_prime_utr

three_prime_utr

CDS

start_codon

stop_codon

exon

intron

transcript

gene

intergenic

Relative coordinates (rel_start, rel_end)
Computed relative to the selected transcript

Strand-aware

Reported as 1-based inclusive

Always refer to the transcript orientation

If transcript spans are missing from the GTF:

Relative coordinates are left empty

annotation_notes explains why

Intron handling
GTF files do not explicitly define introns

Introns are derived as:

ini
Copy code
intron_bp = transcript_overlap - exon_overlap
If transcript features are missing:

Introns are estimated from read length

A warning is added to annotation_notes

Statistics and reporting
Print summary to stdout
bash
Copy code
--report
Example:

yaml
Copy code
Summary: 8,873 annotated / 285 intergenic
Top region_primary counts:
  CDS: 3,421
  three_prime_utr: 2,114
Write statistics to file
bash
Copy code
--stats-out stats.tsv
--top-n 20
The stats file includes:

Annotated vs intergenic counts

Transcript biotypes

Region distributions

Top genes and transcripts

Filtering intergenic regions
Drop intergenic rows
bash
Copy code
--drop-intergenic
(Default behavior keeps all rows.)

GTF caching
To speed up repeated runs, parsed GTFs are cached as Parquet files.

Default cache directory:

bash
Copy code
.cache/genomic-region-annotator/
Custom cache location:

bash
Copy code
--cache-dir /path/to/cache
Debugging a specific row
Print a detailed overlap breakdown for one input row:

bash
Copy code
--debug-row-id 7
--debug-n 50
This prints:

Selected transcript

Overlap per feature

Intron derivation

Final region labels

Species and assemblies
You can annotate any Ensembl-supported species and assembly, as long as a compatible GTF is provided.

When using --release, the tool downloads from Ensembl FTP.

Examples:

Human: homo_sapiens, GRCh38

Mouse: mus_musculus, GRCm39

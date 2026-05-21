# Tiny end-to-end example

This directory contains a small fixture-based workflow that can run without downloading a full Ensembl GTF.

It is useful for checking the command flow and output filenames:

```bash
genomic-region-annotator annotate \
  --input examples/mini_intervals.tsv \
  --gtf examples/mini.gtf \
  --output examples/output/mini.tsv \
  --report

genomic-region-annotator summarize-sites \
  --transcripts examples/output/step1/mini_transcripts.tsv \
  --matrix examples/output/step1/mini_matrix.tsv \
  --input-with-ids examples/output/step1/mini_input_with_ids.tsv \
  --report
```

Expected output files:

```text
examples/output/step1/mini_input_with_ids.tsv
examples/output/step1/mini_matrix.tsv
examples/output/step1/mini_transcripts.tsv
examples/output/step1/mini_step1_stats.tsv
examples/output/step2/mini_site_summary.tsv
examples/output/step2/mini_final.tsv
examples/output/step2/mini_step2_stats.tsv
```

The fixture includes plus- and minus-strand transcripts with CDS and 3' UTR annotations. For real analyses, use `download-gtf` or `--release` with the target Ensembl release instead of this miniature GTF.

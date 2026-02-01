import click

from genomic_region_annotator.annotate import run as run_annotate
from genomic_region_annotator.download import download_gtf
from genomic_region_annotator.summarize_sites import run as run_summarize_sites


@click.group()
def cli():
    """Genomic region annotation utilities."""
    pass


@cli.command("download-gtf")
@click.option("--release", required=True, type=int, help="Ensembl release number (e.g. 111, 112, 115)")
@click.option("--outdir", default="data/annotations", show_default=True, type=click.Path())
@click.option(
    "--species",
    default="homo_sapiens",
    show_default=True,
    help="Ensembl species directory name (e.g. homo_sapiens, mus_musculus).",
)
@click.option(
    "--assembly",
    default=None,
    help="Optional assembly substring to select the correct GTF when multiple exist (e.g. GRCh38, GRCm39).",
)
def download_gtf_cmd(release: int, outdir: str, species: str, assembly: str | None):
    path = download_gtf(release=release, out_dir=outdir, species=species, assembly=assembly)
    click.echo(str(path))


@cli.command("annotate")
@click.option(
    "--input",
    "input_path",
    required=True,
    type=click.Path(exists=True),
    help="Input TSV with genomic coordinates (must include chr,start,end,strand).",
)
@click.option(
    "--output",
    "output_path",
    required=True,
    type=click.Path(),
    help="Output base path (with or without .tsv). Produces <base>_matrix.tsv, <base>_transcripts.tsv, <base>_stats.tsv.",
)

# GTF source (exactly one)
@click.option("--gtf", "gtf_path", type=click.Path(exists=True), help="Path to an Ensembl GTF (alternative to --release).")
@click.option("--release", type=int, help="Ensembl release to download/use (alternative to --gtf).")
@click.option("--outdir", default="data/annotations", show_default=True, type=click.Path(), help="Where to store downloaded GTF.")

# Download selection
@click.option(
    "--species",
    default="homo_sapiens",
    show_default=True,
    help="Ensembl species directory name (used only with --release).",
)
@click.option(
    "--assembly",
    default=None,
    help="Optional assembly substring filter for selecting the correct GTF (used only with --release).",
)

# Coords
@click.option(
    "--coords",
    default="1-based",
    show_default=True,
    type=click.Choice(["1-based", "bed"], case_sensitive=False),
    help="Coordinate convention: '1-based' (inclusive) or 'bed' (0-based half-open).",
)

# Transcript overlap filter
@click.option(
    "--min-overlap-nt",
    default=None,
    type=int,
    help="If set, keep transcripts with >= this many overlapping nucleotides with the transcript span. Default is 100% containment.",
)

# Intergenic handling
@click.option(
    "--drop-intergenic",
    is_flag=True,
    default=False,
    help="Do not output the INTERGENIC row in the matrix output.",
)

# Cache / report
@click.option(
    "--cache-dir",
    default=".cache/genomic-region-annotator",
    show_default=True,
    type=click.Path(),
    help="Directory to store/read cached transcript models.",
)
@click.option(
    "--report",
    is_flag=True,
    default=False,
    help="Print a small report to stdout (also writes stats TSV).",
)

# Debug
@click.option(
    "--debug-row-id",
    type=int,
    default=None,
    help="Debug a specific row (0-based row index in the input TSV).",
)
@click.option("--debug-n", type=int, default=20, show_default=True, help="Reserved (compatibility).")

# Stats
@click.option("--stats-out", default=None, type=click.Path(), help="Write run statistics to this TSV file.")
@click.option("--top-n", default=20, show_default=True, type=int, help="Top N genes to include in stats output.")
def annotate_cmd(
    input_path: str,
    output_path: str,
    gtf_path: str | None,
    release: int | None,
    outdir: str,
    species: str,
    assembly: str | None,
    coords: str,
    min_overlap_nt: int | None,
    drop_intergenic: bool,
    cache_dir: str,
    report: bool,
    debug_row_id: int | None,
    debug_n: int,
    stats_out: str | None,
    top_n: int,
):
    # validate gtf source
    if (gtf_path is None) == (release is None):
        raise click.UsageError("Provide exactly one of --gtf or --release.")

    if release is not None:
        gtf_path = str(download_gtf(release=release, out_dir=outdir, species=species, assembly=assembly))
        click.echo(
            f"[INFO] Using Ensembl GTF (release={release}, species={species}, assembly={assembly or 'auto'}): {gtf_path}"
        )
    else:
        click.echo(f"[INFO] Using provided GTF: {gtf_path}")

    run_annotate(
        input_path=input_path,
        gtf_path=gtf_path,
        output_tsv=output_path,
        coords=coords,
        transcript_first=True,
        min_overlap_nt=min_overlap_nt,
        debug_row_id=debug_row_id,
        debug_n=debug_n,
        drop_intergenic=drop_intergenic,
        cache_dir=cache_dir,
        report=report,
        stats_out=stats_out,
        top_n=top_n,
    )


@cli.command("summarize-sites")
@click.option(
    "--transcripts",
    "transcripts_tsv",
    required=True,
    type=click.Path(exists=True),
    help="Path to <base>_transcripts.tsv produced by the annotate step.",
)
@click.option(
    "--matrix",
    "matrix_tsv",
    required=True,
    type=click.Path(exists=True),
    help="Path to <base>_matrix.tsv produced by the annotate step (used for UNION composition).",
)
@click.option(
    "--output",
    "output_tsv",
    required=False,
    default=None,
    type=click.Path(),
    help="Output path. Default: derive from transcripts path as <base>_site_summary.tsv",
)
@click.option(
    "--policy",
    default="clash_utr3_first",
    show_default=True,
    type=click.Choice(["clash_utr3_first"], case_sensitive=False),
    help="Transcript selection policy (explicit assumptions).",
)
@click.option(
    "--dominance",
    default="coverage",
    show_default=True,
    type=click.Choice(["coverage", "priority"], case_sensitive=False),
    help="How to define dominant region from region composition.",
)
def summarize_sites_cmd(
    transcripts_tsv: str,
    matrix_tsv: str,
    output_tsv: str | None,
    policy: str,
    dominance: str,
):
    run_summarize_sites(
        transcripts_tsv=transcripts_tsv,
        matrix_tsv=matrix_tsv,
        output_tsv=output_tsv,
        policy=policy,
        dominance=dominance,
    )


if __name__ == "__main__":
    cli()

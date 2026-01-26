import click

from genomic_region_annotator.annotate import run as run_annotate
from genomic_region_annotator.download import download_gtf


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
    help="Output TSV with region annotations.",
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

# Intergenic handling
@click.option(
    "--drop-intergenic",
    is_flag=True,
    default=False,
    help="Drop rows with region_primary == 'intergenic' from output.",
)
@click.option(
    "--keep-intergenic/--no-keep-intergenic",
    default=True,
    show_default=True,
    help="Keep intergenic rows in output (default true). Use --drop-intergenic to remove them.",
)

# Cache / report
@click.option(
    "--cache-dir",
    default=".cache/genomic-region-annotator",
    show_default=True,
    type=click.Path(),
    help="Directory to store/read cached GTF parquet.",
)
@click.option(
    "--report",
    is_flag=True,
    default=False,
    help="Print a small report (top region_primary counts).",
)

# Debug
@click.option(
    "--debug-row-id",
    type=int,
    default=None,
    help="Debug a specific row (0-based row index in the input TSV). Prints overlap breakdown.",
)
@click.option("--debug-n", type=int, default=20, show_default=True, help="Max number of overlap rows to print in debug output.")

# Stats
@click.option("--stats-out", default=None, type=click.Path(), help="Write run statistics to this TSV file.")
@click.option("--top-n", default=20, show_default=True, type=int, help="Top N genes/transcripts to include in stats output.")
def annotate_cmd(
    input_path: str,
    output_path: str,
    gtf_path: str | None,
    release: int | None,
    outdir: str,
    species: str,
    assembly: str | None,
    coords: str,
    drop_intergenic: bool,
    keep_intergenic: bool,
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

    # intergenic flags: drop beats keep
    if drop_intergenic:
        keep_intergenic = False

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
        debug_row_id=debug_row_id,
        debug_n=debug_n,
        drop_intergenic=drop_intergenic,
        cache_dir=cache_dir,
        report=report,
        stats_out=stats_out,
        top_n=top_n,
    )


if __name__ == "__main__":
    cli()

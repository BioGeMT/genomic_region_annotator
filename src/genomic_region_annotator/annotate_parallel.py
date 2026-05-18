# src/genomic_region_annotator/annotate_parallel.py
from __future__ import annotations

import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path
from typing import Optional

import pandas as pd

from genomic_region_annotator.annotate import (
    MATRIX_BASE_COLUMNS,
    REGIONS,
    TRANSCRIPT_COLUMNS,
    _append_rows,
    _build_bin_index,
    _empty_output,
    _load_or_build_models_filtered,
    _make_ids,
    _normalize_query_coords,
    _output_base,
    _process_chunk,
    _safe_mean,
    _safe_median,
    log,
)


def run_parallel(
    *,
    input_path: str,
    gtf_path: str,
    output_tsv: str,
    coords: str = "1-based",
    transcript_first: bool = True,
    min_overlap_nt: Optional[int] = None,
    debug_row_id: Optional[int] = None,
    debug_n: int = 20,
    drop_intergenic: bool = False,
    cache_dir: str = ".cache/genomic-region-annotator",
    report: bool = False,
    stats_out: Optional[str] = None,
    top_n: int = 20,
    write_chunk_size: int = 10_000,
    jobs: int = 2,
) -> None:
    """Run Step 1 with bounded chunk-level parallelism.

    This deliberately reuses the same per-read processing function as the serial
    runner and writes chunk results in input order. The matrix and transcript
    schemas are unchanged.
    """
    del transcript_first, debug_n, top_n  # kept for API compatibility

    jobs = max(1, int(jobs))
    if jobs <= 1:
        from genomic_region_annotator.annotate import run

        run(
            input_path=input_path,
            gtf_path=gtf_path,
            output_tsv=output_tsv,
            coords=coords,
            min_overlap_nt=min_overlap_nt,
            debug_row_id=debug_row_id,
            drop_intergenic=drop_intergenic,
            cache_dir=cache_dir,
            report=report,
            stats_out=stats_out,
            write_chunk_size=write_chunk_size,
        )
        return

    base = _output_base(output_tsv)
    base_dir = base.parent if str(base.parent) not in {"", "."} else Path(".")
    stem = base.name

    step1_dir = base_dir / "step1"
    step1_dir.mkdir(parents=True, exist_ok=True)

    input_with_ids_path = step1_dir / f"{stem}_input_with_ids.tsv"
    matrix_path = step1_dir / f"{stem}_matrix.tsv"
    transcripts_path = step1_dir / f"{stem}_transcripts.tsv"
    stats_path = Path(stats_out) if stats_out else (step1_dir / f"{stem}_step1_stats.tsv")

    log("Stage 1/7: read input intervals")
    df = pd.read_csv(Path(input_path), sep="\t", dtype={"chr": "string", "strand": "string"})
    required = {"chr", "start", "end", "strand"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Input TSV missing required columns: {sorted(missing)}")

    if "id" in df.columns:
        df = df.rename(columns={"id": "original_id"})

    df.insert(0, "id", _make_ids(len(df)))
    df["strand"] = df["strand"].fillna(".").astype("string")
    df.loc[~df["strand"].isin(["+", "-", "."]), "strand"] = "."

    starts_1b: list[int] = []
    ends_1b: list[int] = []
    lens: list[int] = []
    for _, r in df.iterrows():
        s, e = _normalize_query_coords(int(r["start"]), int(r["end"]), coords=coords)
        starts_1b.append(s)
        ends_1b.append(e)
        lens.append(e - s + 1)

    df["start_1based"] = starts_1b
    df["end_1based"] = ends_1b
    df["read_len"] = lens

    log(f"Writing: {input_with_ids_path}")
    df.to_csv(input_with_ids_path, sep="\t", index=False)

    chroms_needed = set(df["chr"].dropna().astype(str).unique().tolist())
    log(f"Input rows: {len(df):,}")
    log(f"Chromosomes in input: {sorted(chroms_needed)}")

    log("Stage 2/7: load/build transcript models (filtered)")
    models = _load_or_build_models_filtered(gtf_path=gtf_path, cache_dir=cache_dir, chroms_needed=chroms_needed)

    log("Stage 3/7: build bin index")
    bin_size = 100_000
    index = _build_bin_index(models=models, bin_size=bin_size)

    if min_overlap_nt is None:
        log("Transcript overlap filter: 100% containment (default)")
    else:
        log(f"Transcript overlap filter: overlap >= {min_overlap_nt} nt")

    log(f"Stage 4/7: per-read transcript filtering + matrix build (parallel jobs={jobs})")
    t0 = time.time()

    n = len(df)
    progress_every = 2_000 if n < 50_000 else 10_000
    write_chunk_size = max(1, int(write_chunk_size))
    max_read_len = int(df["read_len"].max()) if n else 0
    matrix_columns = MATRIX_BASE_COLUMNS + [f"nt_{i}" for i in range(1, max_read_len + 1)]

    matrix_header = True
    tx_header = True
    matrix_wrote_rows = False
    tx_wrote_rows = False

    for p in [matrix_path, transcripts_path]:
        if p.exists():
            p.unlink()

    tx_count_per_read: list[int] = []
    reads_with_any_tx = 0
    reads_with_any_cds = 0
    reads_with_any_utr3 = 0
    reads_with_any_utr5 = 0
    reads_with_any_exon = 0
    reads_with_any_intron = 0
    region_total_ones = {r: 0 for r in REGIONS}

    def make_chunk(start: int) -> pd.DataFrame:
        return df.iloc[start : start + write_chunk_size].copy()

    def submit(executor: ThreadPoolExecutor, start: int):
        chunk = make_chunk(start)
        return executor.submit(
            _process_chunk,
            chunk,
            models=models,
            index=index,
            bin_size=bin_size,
            min_overlap_nt=min_overlap_nt,
            drop_intergenic=drop_intergenic,
            debug_row_id=debug_row_id,
        )

    starts = list(range(0, n, write_chunk_size))
    next_submit = 0
    next_write = 0
    processed = 0
    buffered = {}
    futures = {}

    with ThreadPoolExecutor(max_workers=jobs) as executor:
        while next_submit < len(starts) and len(futures) < jobs:
            fut = submit(executor, starts[next_submit])
            futures[fut] = next_submit
            next_submit += 1

        while futures:
            done, _ = wait(futures.keys(), return_when=FIRST_COMPLETED)
            for fut in done:
                idx = futures.pop(fut)
                buffered[idx] = fut.result()
                if next_submit < len(starts):
                    new_fut = submit(executor, starts[next_submit])
                    futures[new_fut] = next_submit
                    next_submit += 1

            while next_write in buffered:
                result = buffered.pop(next_write)
                matrix_header = _append_rows(matrix_path, result.matrix_rows, matrix_columns, write_header=matrix_header)
                if result.matrix_rows:
                    matrix_wrote_rows = True
                tx_header = _append_rows(transcripts_path, result.tx_rows, TRANSCRIPT_COLUMNS, write_header=tx_header)
                if result.tx_rows:
                    tx_wrote_rows = True

                tx_count_per_read.extend(result.tx_count_per_read)
                reads_with_any_tx += result.reads_with_any_tx
                reads_with_any_cds += result.reads_with_any_cds
                reads_with_any_utr3 += result.reads_with_any_utr3
                reads_with_any_utr5 += result.reads_with_any_utr5
                reads_with_any_exon += result.reads_with_any_exon
                reads_with_any_intron += result.reads_with_any_intron
                for region, ones in result.region_total_ones.items():
                    region_total_ones[region] += ones
                for msg in result.debug_logs:
                    log(msg)

                processed += len(result.tx_count_per_read)
                if processed % progress_every == 0 or processed == n:
                    elapsed = time.time() - t0
                    rate = processed / elapsed if elapsed > 0 else 0.0
                    log(f"Processed {processed:,}/{n:,} reads ({rate:,.1f} reads/s)")
                next_write += 1

    log("Stage 5/7: finalize streamed outputs")
    if not matrix_wrote_rows:
        _empty_output(matrix_path, matrix_columns)
    if not tx_wrote_rows:
        _empty_output(transcripts_path, TRANSCRIPT_COLUMNS)

    log("Stage 6/7: compute stats + write stats TSV")
    total_reads = len(df)
    total_nt = int(sum(df["read_len"].astype(int).tolist()))
    region_frac = {str(region): float(ones / total_nt) if total_nt else 0.0 for region, ones in region_total_ones.items()}

    stats_rows: list[dict[str, object]] = [
        {"metric": "step", "value": "step1_evidence_only"},
        {"metric": "output_stem", "value": stem},
        {"metric": "total_reads", "value": total_reads},
        {"metric": "total_nt", "value": total_nt},
        {"metric": "min_overlap_nt", "value": "None(100% containment)" if min_overlap_nt is None else int(min_overlap_nt)},
        {"metric": "reads_with_any_transcript", "value": reads_with_any_tx},
        {"metric": "reads_with_no_transcript", "value": total_reads - reads_with_any_tx},
        {"metric": "mean_transcripts_per_read", "value": round(_safe_mean(tx_count_per_read), 4)},
        {"metric": "median_transcripts_per_read", "value": round(_safe_median(tx_count_per_read), 4)},
        {"metric": "max_transcripts_per_read", "value": max(tx_count_per_read) if tx_count_per_read else 0},
        {"metric": "reads_with_any_CDS_union", "value": reads_with_any_cds},
        {"metric": "reads_with_any_UTR3_union", "value": reads_with_any_utr3},
        {"metric": "reads_with_any_UTR5_union", "value": reads_with_any_utr5},
        {"metric": "reads_with_any_EXON_union", "value": reads_with_any_exon},
        {"metric": "reads_with_any_INTRON_union", "value": reads_with_any_intron},
    ]

    for region in REGIONS:
        if drop_intergenic and region == "INTERGENIC":
            continue
        stats_rows.append({"metric": f"total_ones_{region}", "value": region_total_ones[region]})
        stats_rows.append({"metric": f"fraction_ones_{region}", "value": round(region_frac[region], 6)})

    pd.DataFrame(stats_rows).to_csv(stats_path, sep="\t", index=False)

    if report:
        log("Stage 7/7: report")
        log(f"Reads with >=1 passing transcript: {reads_with_any_tx:,}/{total_reads:,}")
        log(f"Mean transcripts/read: {_safe_mean(tx_count_per_read):.3f} (median={_safe_median(tx_count_per_read):.1f}, max={max(tx_count_per_read) if tx_count_per_read else 0})")
        if region_frac:
            log("Region coverage fractions (UNION; across all nts):")
            for region in REGIONS:
                if drop_intergenic and region == "INTERGENIC":
                    continue
                if region in region_frac:
                    log(f"  {region}: {region_frac[region]:.4f}")

# src/genomic_region_annotator/annotate.py
from __future__ import annotations

import hashlib
import statistics
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Optional

import pandas as pd


# ----------------------------
# Logging
# ----------------------------
def _ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str) -> None:
    print(f"[{_ts()}] [INFO] {msg}", flush=True)


# ----------------------------
# Transcript model
# ----------------------------
@dataclass(frozen=True)
class TranscriptModel:
    chrom: str
    strand: str  # '+' or '-'
    transcript_id: str
    gene_id: str
    gene_name: str

    tx_start: int  # 1-based inclusive
    tx_end: int    # 1-based inclusive

    exons: tuple[tuple[int, int], ...]  # 1-based inclusive merged
    cds: tuple[tuple[int, int], ...]    # 1-based inclusive merged
    utr5: tuple[tuple[int, int], ...]   # 1-based inclusive merged (may be inferred)
    utr3: tuple[tuple[int, int], ...]   # 1-based inclusive merged (may be inferred)


# ----------------------------
# Interval helpers (1-based inclusive)
# ----------------------------
def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ps, pe = merged[-1]
        if s <= pe + 1:
            merged[-1] = (ps, max(pe, e))
        else:
            merged.append((s, e))
    return merged


def _within(q: tuple[int, int], iv: tuple[int, int]) -> bool:
    return q[0] >= iv[0] and q[1] <= iv[1]


def _ovl_bp(a: tuple[int, int], b: tuple[int, int]) -> int:
    s = max(a[0], b[0])
    e = min(a[1], b[1])
    return max(0, e - s + 1)


def _sum_ovl_bp(q: tuple[int, int], ivs: tuple[tuple[int, int], ...]) -> int:
    return sum(_ovl_bp(q, iv) for iv in ivs)


def _iter_bins(start: int, end: int, bin_size: int) -> Iterable[int]:
    return range(start // bin_size, end // bin_size + 1)


def _pos_in_intervals(pos: int, intervals: tuple[tuple[int, int], ...]) -> bool:
    # intervals are merged & sorted; linear scan is fine (exons/CDS/UTRs per transcript are small)
    for s, e in intervals:
        if pos < s:
            return False
        if s <= pos <= e:
            return True
    return False


# ----------------------------
# UTR inference (if missing)
# ----------------------------
def _infer_utrs_from_exons_and_cds(
    exons: list[tuple[int, int]],
    cds: list[tuple[int, int]],
    strand: str,
) -> tuple[list[tuple[int, int]], list[tuple[int, int]]]:
    """
    Infer UTRs from exon/CDS when explicit UTR features are absent.
      '+' strand: utr5 = exonic bases < CDS_start; utr3 = exonic bases > CDS_end
      '-' strand: utr5 = exonic bases > CDS_end; utr3 = exonic bases < CDS_start
    """
    ex_m = _merge_intervals(exons)
    cds_m = _merge_intervals(cds)
    if not ex_m or not cds_m:
        return [], []

    cds_start = min(s for s, _ in cds_m)
    cds_end = max(e for _, e in cds_m)

    utr5: list[tuple[int, int]] = []
    utr3: list[tuple[int, int]] = []

    def clip_left(iv: tuple[int, int], bound: int) -> Optional[tuple[int, int]]:
        s, e = iv
        if e < bound:
            return (s, e)
        if s < bound <= e:
            return (s, bound - 1)
        return None

    def clip_right(iv: tuple[int, int], bound: int) -> Optional[tuple[int, int]]:
        s, e = iv
        if s > bound:
            return (s, e)
        if s <= bound < e:
            return (bound + 1, e)
        return None

    for ex in ex_m:
        left = clip_left(ex, cds_start)
        right = clip_right(ex, cds_end)
        if strand == "+":
            if left:
                utr5.append(left)
            if right:
                utr3.append(right)
        else:
            if right:
                utr5.append(right)
            if left:
                utr3.append(left)

    return _merge_intervals(utr5), _merge_intervals(utr3)


# ----------------------------
# GTF loading (filtered, chunked, fast attribute extraction)
# ----------------------------
def _cache_key_for_gtf(gtf_path: str) -> str:
    p = Path(gtf_path)
    st = p.stat()
    s = f"{p.resolve()}|{st.st_size}|{int(st.st_mtime)}"
    return hashlib.sha256(s.encode("utf-8")).hexdigest()[:16]


def _load_gtf_dataframe_filtered_fast(gtf_path: str, chroms_needed: set[str]) -> pd.DataFrame:
    """
    Load only GTF rows on chromosomes present in input.
    Keeps features: exon, CDS, five_prime_utr, three_prime_utr.
    Extracts transcript_id / gene_id / gene_name via vectorized regex.
    """
    cols = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
    wanted_features = {"exon", "CDS", "five_prime_utr", "three_prime_utr"}

    log("Reading GTF in chunks (filtered to input chromosomes)...")
    t0 = time.time()

    chunks: list[pd.DataFrame] = []
    chunksize = 2_000_000
    kept_total = 0

    reader = pd.read_csv(
        gtf_path,
        sep="\t",
        comment="#",
        header=None,
        names=cols,
        usecols=[0, 2, 3, 4, 6, 8],  # seqname, feature, start, end, strand, attribute
        dtype={
            "seqname": "string",
            "feature": "string",
            "start": "int64",
            "end": "int64",
            "strand": "string",
            "attribute": "string",
        },
        chunksize=chunksize,
    )

    for k, chunk in enumerate(reader, start=1):
        chunk = chunk.loc[
            chunk["seqname"].astype(str).isin(chroms_needed)
            & chunk["feature"].isin(wanted_features)
        ].copy()

        if chunk.empty:
            if k % 5 == 0:
                log(f"GTF chunk {k}: kept_rows_total={kept_total:,}")
            continue

        attr = chunk["attribute"].astype(str)
        chunk["transcript_id"] = attr.str.extract(r'transcript_id "([^"]+)"', expand=False)
        chunk["gene_id"] = attr.str.extract(r'gene_id "([^"]+)"', expand=False)
        chunk["gene_name"] = attr.str.extract(r'gene_name "([^"]+)"', expand=False)
        chunk = chunk.loc[chunk["transcript_id"].notna()].copy()

        keep = chunk[["seqname", "strand", "feature", "start", "end", "transcript_id", "gene_id", "gene_name"]]
        chunks.append(keep)
        kept_total += len(keep)

        if k % 2 == 0:
            log(f"GTF chunk {k}: kept_rows_total={kept_total:,}")

    if not chunks:
        log(f"No matching GTF rows found for chroms={sorted(chroms_needed)}")
        return pd.DataFrame(columns=["seqname", "strand", "feature", "start", "end", "transcript_id", "gene_id", "gene_name"])

    df = pd.concat(chunks, ignore_index=True)
    log(f"Filtered GTF rows kept: {len(df):,} (took {time.time()-t0:.1f}s)")
    return df


def _build_transcript_models_single_pass(gtf_df: pd.DataFrame) -> list[TranscriptModel]:
    """
    Build TranscriptModels by sorting by transcript_id and scanning once.
    Transcript span is derived from exons.
    """
    log("Building transcript models (single-pass over transcript_id)...")
    t0 = time.time()

    if gtf_df.empty:
        return []

    gtf_df = gtf_df.sort_values(["transcript_id", "start", "end"], kind="mergesort")

    models: list[TranscriptModel] = []

    cur_tid: Optional[str] = None
    cur_chrom = ""
    cur_strand = ""
    cur_gene_id = ""
    cur_gene_name = ""
    exons: list[tuple[int, int]] = []
    cds: list[tuple[int, int]] = []
    utr5: list[tuple[int, int]] = []
    utr3: list[tuple[int, int]] = []

    def flush() -> None:
        nonlocal cur_tid, cur_chrom, cur_strand, cur_gene_id, cur_gene_name, exons, cds, utr5, utr3
        if not cur_tid or not exons:
            exons, cds, utr5, utr3 = [], [], [], []
            return

        ex_m = _merge_intervals(exons)
        tx_start = min(s for s, _ in ex_m)
        tx_end = max(e for _, e in ex_m)

        cds_m = _merge_intervals(cds)
        utr5_m = _merge_intervals(utr5)
        utr3_m = _merge_intervals(utr3)

        if (not utr5_m and not utr3_m) and ex_m and cds_m:
            utr5_m, utr3_m = _infer_utrs_from_exons_and_cds(ex_m, cds_m, cur_strand)

        models.append(
            TranscriptModel(
                chrom=cur_chrom,
                strand=cur_strand,
                transcript_id=cur_tid,
                gene_id=cur_gene_id or "",
                gene_name=cur_gene_name or "",
                tx_start=tx_start,
                tx_end=tx_end,
                exons=tuple(ex_m),
                cds=tuple(cds_m),
                utr5=tuple(utr5_m),
                utr3=tuple(utr3_m),
            )
        )

        exons, cds, utr5, utr3 = [], [], [], []

    for i, r in enumerate(gtf_df.itertuples(index=False), start=1):
        tid = str(r.transcript_id)
        if cur_tid is None:
            cur_tid = tid
            cur_chrom = str(r.seqname)
            cur_strand = str(r.strand)
            cur_gene_id = "" if pd.isna(r.gene_id) else str(r.gene_id)
            cur_gene_name = "" if pd.isna(r.gene_name) else str(r.gene_name)

        if tid != cur_tid:
            flush()
            cur_tid = tid
            cur_chrom = str(r.seqname)
            cur_strand = str(r.strand)
            cur_gene_id = "" if pd.isna(r.gene_id) else str(r.gene_id)
            cur_gene_name = "" if pd.isna(r.gene_name) else str(r.gene_name)

        s, e = int(r.start), int(r.end)
        feat = str(r.feature)
        if feat == "exon":
            exons.append((s, e))
        elif feat == "CDS":
            cds.append((s, e))
        elif feat == "five_prime_utr":
            utr5.append((s, e))
        elif feat == "three_prime_utr":
            utr3.append((s, e))

        if i % 2_000_000 == 0:
            log(f"Model build scan: processed {i:,} filtered GTF rows...")

    flush()
    log(f"Built transcript models: {len(models):,} (took {time.time()-t0:.1f}s)")
    return models


def _load_or_build_models_filtered(gtf_path: str, cache_dir: str, chroms_needed: set[str]) -> list[TranscriptModel]:
    cache_base = Path(cache_dir)
    cache_base.mkdir(parents=True, exist_ok=True)

    chrom_key = hashlib.sha256(",".join(sorted(chroms_needed)).encode("utf-8")).hexdigest()[:10]
    key = f"{_cache_key_for_gtf(gtf_path)}_{chrom_key}"
    pkl_path = cache_base / f"models_filtered_{key}.pkl"

    if pkl_path.exists():
        log(f"Loading cached transcript models: {pkl_path}")
        return pd.read_pickle(pkl_path)

    gtf_df = _load_gtf_dataframe_filtered_fast(gtf_path, chroms_needed=chroms_needed)
    models = _build_transcript_models_single_pass(gtf_df)

    log(f"Saving cache: {pkl_path}")
    pd.to_pickle(models, pkl_path)
    return models


# ----------------------------
# Candidate retrieval via bin index
# ----------------------------
def _build_bin_index(models: list[TranscriptModel], bin_size: int = 100_000) -> dict[tuple[str, str, int], list[int]]:
    """
    Map (chrom, strand, bin) -> transcript indices.
    Index transcript span only (we're requiring 100% containment by default).
    """
    log(f"Building bin index (bin_size={bin_size})...")
    t0 = time.time()
    index: dict[tuple[str, str, int], list[int]] = {}
    for i, m in enumerate(models):
        for b in _iter_bins(m.tx_start, m.tx_end, bin_size):
            index.setdefault((m.chrom, m.strand, b), []).append(i)
    log(f"Bin index keys: {len(index):,} (took {time.time()-t0:.1f}s)")
    return index


def _candidate_indices(chrom: str, qstrand: str, qiv: tuple[int, int], index: dict[tuple[str, str, int], list[int]], bin_size: int) -> set[int]:
    s, e = qiv
    bins = list(_iter_bins(s, e, bin_size))
    cand: set[int] = set()

    if qstrand in {"+", "-"}:
        for b in bins:
            cand.update(index.get((chrom, qstrand, b), []))
    else:
        for b in bins:
            cand.update(index.get((chrom, "+", b), []))
            cand.update(index.get((chrom, "-", b), []))
    return cand


# ----------------------------
# Matrix construction
# ----------------------------
REGIONS = ["CDS", "UTR5", "UTR3", "EXON", "INTRON", "INTERGENIC"]


def _normalize_query_coords(start: int, end: int, coords: str) -> tuple[int, int]:
    if coords.lower() == "1-based":
        s, e = int(start), int(end)
        if s > e:
            s, e = e, s
        return s, e
    if coords.lower() == "bed":
        s0, e0 = int(start), int(end)
        if s0 > e0:
            s0, e0 = e0, s0
        s1 = s0 + 1
        e1 = e0
        if e1 < s1:
            e1 = s1
        return s1, e1
    raise ValueError(f"Unknown coords convention: {coords}")


def _per_nt_region_flags_for_transcript(m: TranscriptModel, qiv: tuple[int, int]) -> dict[str, list[int]]:
    """
    For a transcript, compute per-nt flags on the query interval.
    Raw sets:
      - CDS, UTR5, UTR3, EXON (any exon), INTRON (within tx span but not exon)
    """
    s, e = qiv
    L = e - s + 1
    cds_flags = [0] * L
    utr5_flags = [0] * L
    utr3_flags = [0] * L
    exon_flags = [0] * L
    intron_flags = [0] * L

    for i in range(L):
        pos = s + i
        in_cds = _pos_in_intervals(pos, m.cds)
        in_utr5 = _pos_in_intervals(pos, m.utr5)
        in_utr3 = _pos_in_intervals(pos, m.utr3)
        in_exon = _pos_in_intervals(pos, m.exons)
        in_intron = (m.tx_start <= pos <= m.tx_end) and (not in_exon)

        cds_flags[i] = 1 if in_cds else 0
        utr5_flags[i] = 1 if in_utr5 else 0
        utr3_flags[i] = 1 if in_utr3 else 0
        exon_flags[i] = 1 if in_exon else 0
        intron_flags[i] = 1 if in_intron else 0

    return {"CDS": cds_flags, "UTR5": utr5_flags, "UTR3": utr3_flags, "EXON": exon_flags, "INTRON": intron_flags}


def _or_region_matrices(mats: list[dict[str, list[int]]], L: int) -> dict[str, list[int]]:
    """
    Union across transcripts: OR per cell.
    INTERGENIC is 1 where no other region is 1.
    """
    out = {r: [0] * L for r in REGIONS}
    if not mats:
        out["INTERGENIC"] = [1] * L
        return out

    for r in ["CDS", "UTR5", "UTR3", "EXON", "INTRON"]:
        for i in range(L):
            out[r][i] = 1 if any(m[r][i] == 1 for m in mats) else 0

    for i in range(L):
        out["INTERGENIC"][i] = 1 if (
            out["CDS"][i] == 0
            and out["UTR5"][i] == 0
            and out["UTR3"][i] == 0
            and out["EXON"][i] == 0
            and out["INTRON"][i] == 0
        ) else 0

    return out


def _make_ids(n: int) -> list[str]:
    width = max(3, len(str(n)))
    return [f"{i+1:0{width}d}" for i in range(n)]


def _output_base(path: str) -> Path:
    p = Path(path)
    if p.suffix.lower() == ".tsv":
        return p.with_suffix("")
    return p


def _safe_mean(xs: list[int]) -> float:
    return float(sum(xs) / len(xs)) if xs else 0.0


def _safe_median(xs: list[int]) -> float:
    return float(statistics.median(xs)) if xs else 0.0


# ----------------------------
# Public API (CLI entrypoint)
# ----------------------------
def run(
    *,
    input_path: str,
    gtf_path: str,
    output_tsv: str,
    coords: str = "1-based",
    transcript_first: bool = True,  # kept for CLI compatibility; transcript choice comes next
    min_overlap_nt: Optional[int] = None,  # None => require 100% containment; else keep transcripts with >= this many bp overlap with transcript span
    debug_row_id: Optional[int] = None,
    debug_n: int = 20,  # unused; kept for CLI compatibility
    drop_intergenic: bool = False,  # if True, do not output INTERGENIC rows in matrix
    cache_dir: str = ".cache/genomic-region-annotator",
    report: bool = False,
    stats_out: Optional[str] = None,  # if None, writes <output_base>_stats.tsv
    top_n: int = 20,  # top genes in stats output
) -> None:
    """
    Outputs:

    1) input_with_ids.tsv
        - input rows with generated ids: 001, 002, ...
        - if input already had an 'id' column, it is preserved as 'original_id'

    2) <output_base>_matrix.tsv
        - one row per (id, region), columns nt_1..nt_L
        - cell = 1 if at least one passing transcript assigns that base to that region

    3) <output_base>_transcripts.tsv
        - one row per (id, transcript) where transcript passes the overlap filter

    4) <output_base>_stats.tsv (or stats_out if provided)
        - run summary + some dataset-level metrics

    Transcript overlap filter (strand + chromosome must match):
      - Default (min_overlap_nt is None): require 100% containment (read fully within transcript span [tx_start, tx_end]).
      - If min_overlap_nt is set: keep transcripts with overlap(read, transcript span) >= min_overlap_nt.
    """
    base = _output_base(output_tsv)
    out_dir = base.parent if str(base.parent) != "" else Path(".")
    out_dir.mkdir(parents=True, exist_ok=True)

    input_with_ids_path = out_dir / "input_with_ids.tsv"
    matrix_path = Path(str(base) + "_matrix.tsv")
    transcripts_path = Path(str(base) + "_transcripts.tsv")
    stats_path = Path(stats_out) if stats_out else Path(str(base) + "_stats.tsv")

    log("Stage 1/7: read input intervals")
    df = pd.read_csv(Path(input_path), sep="\t", dtype={"chr": "string", "strand": "string"})
    required = {"chr", "start", "end", "strand"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Input TSV missing required columns: {sorted(missing)}")

    if "id" in df.columns:
        df = df.rename(columns={"id": "original_id"})

    df.insert(0, "id", _make_ids(len(df)))
    log(f"Assigned ids: {df['id'].iloc[0]} .. {df['id'].iloc[-1]}")

    df["strand"] = df["strand"].fillna(".").astype("string")
    df.loc[~df["strand"].isin(["+", "-", "."]), "strand"] = "."

    # Normalize coords and compute read length for the id-anchored input file
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
    log(f"Transcript models loaded: {len(models):,}")

    log("Stage 3/7: build bin index")
    bin_size = 100_000
    index = _build_bin_index(models=models, bin_size=bin_size)

    if min_overlap_nt is None:
        log("Transcript overlap filter: 100% containment (default)")
    else:
        log(f"Transcript overlap filter: overlap >= {min_overlap_nt} nt")

    log("Stage 4/7: per-read transcript filtering + matrix build")
    t0 = time.time()

    matrix_rows: list[dict[str, Any]] = []
    tx_rows: list[dict[str, Any]] = []

    n = len(df)
    progress_every = 2_000 if n < 50_000 else 10_000

    # Stats accumulators
    tx_count_per_read: list[int] = []
    reads_with_any_tx = 0
    reads_with_any_cds = 0
    reads_with_any_utr3 = 0
    reads_with_any_utr5 = 0
    reads_with_any_exon = 0
    reads_with_any_intron = 0

    for row_i, row in df.iterrows():
        rid = str(row["id"])
        chrom = str(row["chr"])
        qstrand = str(row["strand"])
        s = int(row["start_1based"])
        e = int(row["end_1based"])
        qiv = (s, e)
        L = int(row["read_len"])

        cand = _candidate_indices(chrom, qstrand, qiv, index=index, bin_size=bin_size)

        passing: list[TranscriptModel] = []
        for idx in cand:
            m = models[idx]
            if m.chrom != chrom:
                continue
            if qstrand in {"+", "-"} and m.strand != qstrand:
                continue
            tx_span = (m.tx_start, m.tx_end)
            ov_tx = _ovl_bp(qiv, tx_span)
            if min_overlap_nt is None:
                if _within(qiv, tx_span):
                    passing.append(m)
            else:
                if ov_tx >= int(min_overlap_nt):
                    passing.append(m)

        tx_count_per_read.append(len(passing))
        if passing:
            reads_with_any_tx += 1

        # transcript report
        for m in passing:
            tx_rows.append(
                {
                    "id": rid,
                    "chr": chrom,
                    "start": s,
                    "end": e,
                    "strand": qstrand,
                    "read_len": L,
                    "transcript_id": m.transcript_id,
                    "gene_id": m.gene_id,
                    "gene_name": m.gene_name,
                    "transcript_strand": m.strand,
                    "tx_start": m.tx_start,
                    "tx_end": m.tx_end,
                    "overlap_tx_bp": _ovl_bp(qiv, (m.tx_start, m.tx_end)),
                    "contained_100pct": 1 if _within(qiv, (m.tx_start, m.tx_end)) else 0,
                    "overlap_exon_bp": _sum_ovl_bp(qiv, m.exons),
                    "overlap_cds_bp": _sum_ovl_bp(qiv, m.cds),
                    "overlap_utr5_bp": _sum_ovl_bp(qiv, m.utr5),
                    "overlap_utr3_bp": _sum_ovl_bp(qiv, m.utr3),
                }
            )

        per_tx_mats = [_per_nt_region_flags_for_transcript(m, qiv) for m in passing]
        union = _or_region_matrices(per_tx_mats, L=L)

        # read-level region presence stats (union)
        if sum(union["CDS"]) > 0:
            reads_with_any_cds += 1
        if sum(union["UTR3"]) > 0:
            reads_with_any_utr3 += 1
        if sum(union["UTR5"]) > 0:
            reads_with_any_utr5 += 1
        if sum(union["EXON"]) > 0:
            reads_with_any_exon += 1
        if sum(union["INTRON"]) > 0:
            reads_with_any_intron += 1

        for region in REGIONS:
            if drop_intergenic and region == "INTERGENIC":
                continue
            r = {"id": rid, "chr": chrom, "start": s, "end": e, "strand": qstrand, "region": region, "read_len": L}
            for j in range(L):
                r[f"nt_{j+1}"] = union[region][j]
            matrix_rows.append(r)

        if debug_row_id is not None and row_i == debug_row_id:
            log(f"[DEBUG] id={rid} {chrom}:{s}-{e} strand={qstrand} L={L}")
            log(f"[DEBUG] candidate_transcripts_in_bins={len(cand)} passing_filter={len(passing)}")
            log("[DEBUG] union region sums:")
            for region in REGIONS:
                log(f"  {region}: {sum(union[region])}/{L}")

        if (row_i + 1) % progress_every == 0:
            elapsed = time.time() - t0
            rate = (row_i + 1) / elapsed if elapsed > 0 else 0.0
            log(f"Processed {row_i+1:,}/{n:,} reads ({rate:,.1f} reads/s)")

    log("Stage 5/7: write outputs")
    out_df = pd.DataFrame(matrix_rows)

    # Fill missing nt columns (because reads can have different lengths)
    nt_cols = [c for c in out_df.columns if c.startswith("nt_")]
    if nt_cols:
        out_df[nt_cols] = out_df[nt_cols].fillna(0).astype(int)

    out_df.to_csv(matrix_path, sep="\t", index=False)
    log(f"Wrote matrix: {matrix_path} (rows={len(out_df):,})")

    tx_df = pd.DataFrame(tx_rows)
    tx_df.to_csv(transcripts_path, sep="\t", index=False)
    log(f"Wrote transcripts: {transcripts_path} (rows={len(tx_df):,})")

    log("Stage 6/7: compute stats + write stats TSV")
    read_lens = df["read_len"].astype(int).tolist()
    total_reads = len(df)
    total_nt = int(sum(read_lens))

    # Region coverage totals from matrix
    region_total_ones: dict[str, int] = {}
    region_frac: dict[str, float] = {}
    if not out_df.empty and nt_cols:
        for region, grp in out_df.groupby("region"):
            ones = int(grp[nt_cols].sum().sum())
            region_total_ones[str(region)] = ones
            region_frac[str(region)] = float(ones / total_nt) if total_nt else 0.0

    # Top genes by transcript hits
    top_genes: list[tuple[str, int]] = []
    if not tx_df.empty and "gene_name" in tx_df.columns:
        gcounts = tx_df["gene_name"].fillna("").astype(str)
        gcounts = gcounts.where(gcounts != "", tx_df["gene_id"].fillna("").astype(str))
        top = gcounts.value_counts().head(max(1, int(top_n)))
        top_genes = list(zip(top.index.tolist(), top.values.tolist()))

    stats_rows: list[dict[str, Any]] = []
    stats_rows.append({"metric": "total_reads", "value": total_reads})
    stats_rows.append({"metric": "total_nt", "value": total_nt})
    stats_rows.append({"metric": "min_overlap_nt", "value": "None(100% containment)" if min_overlap_nt is None else int(min_overlap_nt)})
    stats_rows.append({"metric": "reads_with_any_transcript", "value": reads_with_any_tx})
    stats_rows.append({"metric": "reads_with_no_transcript", "value": total_reads - reads_with_any_tx})
    stats_rows.append({"metric": "mean_transcripts_per_read", "value": round(_safe_mean(tx_count_per_read), 4)})
    stats_rows.append({"metric": "median_transcripts_per_read", "value": round(_safe_median(tx_count_per_read), 4)})
    stats_rows.append({"metric": "max_transcripts_per_read", "value": max(tx_count_per_read) if tx_count_per_read else 0})

    stats_rows.append({"metric": "reads_with_any_CDS_union", "value": reads_with_any_cds})
    stats_rows.append({"metric": "reads_with_any_UTR3_union", "value": reads_with_any_utr3})
    stats_rows.append({"metric": "reads_with_any_UTR5_union", "value": reads_with_any_utr5})
    stats_rows.append({"metric": "reads_with_any_EXON_union", "value": reads_with_any_exon})
    stats_rows.append({"metric": "reads_with_any_INTRON_union", "value": reads_with_any_intron})

    for region in REGIONS:
        if drop_intergenic and region == "INTERGENIC":
            continue
        if region in region_total_ones:
            stats_rows.append({"metric": f"total_ones_{region}", "value": region_total_ones[region]})
            stats_rows.append({"metric": f"fraction_ones_{region}", "value": round(region_frac[region], 6)})

    # Add top genes as separate rows (wide in a simple TSV)
    for i, (g, c) in enumerate(top_genes, start=1):
        stats_rows.append({"metric": f"top_gene_{i}", "value": g})
        stats_rows.append({"metric": f"top_gene_{i}_transcript_rows", "value": c})

    stats_df = pd.DataFrame(stats_rows)
    stats_df.to_csv(stats_path, sep="\t", index=False)
    log(f"Wrote stats: {stats_path}")

    if report:
        log("Stage 7/7: report")
        log(f"Reads with >=1 passing transcript: {reads_with_any_tx:,}/{total_reads:,}")
        log(f"Mean transcripts/read: {_safe_mean(tx_count_per_read):.3f}  (median={_safe_median(tx_count_per_read):.1f}, max={max(tx_count_per_read) if tx_count_per_read else 0})")
        if region_frac:
            log("Region coverage fractions (union; across all nts):")
            for region in REGIONS:
                if drop_intergenic and region == "INTERGENIC":
                    continue
                if region in region_frac:
                    log(f"  {region}: {region_frac[region]:.4f}")

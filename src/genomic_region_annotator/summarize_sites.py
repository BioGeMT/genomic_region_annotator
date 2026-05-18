# src/genomic_region_annotator/summarize_sites.py
from __future__ import annotations

import time
from collections import Counter
from pathlib import Path
from typing import Any, Optional

import pandas as pd


def _ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str) -> None:
    print(f"[{_ts()}] [INFO] {msg}", flush=True)


PRIORITY = ["UTR3", "CDS", "UTR5", "EXON_OTHER", "INTRON", "INTERGENIC"]

SUMMARY_COLUMNS = [
    "id",
    "read_len",
    "policy",
    "dominance",
    "selected_transcript_id",
    "selected_gene_id",
    "selected_gene_name",
    "selected_tx_start",
    "selected_tx_end",
    "selected_read_start_in_tx_1based",
    "selected_read_end_in_tx_1based",
    "selected_overlap_start_genome_1based",
    "selected_overlap_end_genome_1based",
    "selected_overlap_start_in_tx_1based",
    "selected_overlap_end_in_tx_1based",
    "dominant_region_selected",
    "regions_present_selected",
    "bp_utr3_selected",
    "bp_cds_selected",
    "bp_utr5_selected",
    "bp_exon_other_selected",
    "bp_intron_selected",
    "bp_intergenic_selected",
    "dominant_region_union",
    "regions_present_union",
    "bp_utr3_union",
    "bp_cds_union",
    "bp_utr5_union",
    "bp_exon_other_union",
    "bp_intron_union",
    "bp_intergenic_union",
    "ambiguous_union_vs_selected",
    "n_passing_transcripts",
]

FINAL_COLUMNS = [
    "id",
    "chr",
    "start",
    "end",
    "strand",
    "selected_gene_name",
    "selected_transcript_id",
    "selected_read_start_in_tx_1based",
    "selected_read_end_in_tx_1based",
    "dominant_region_selected",
    "regions_present_selected",
    "ambiguous_union_vs_selected",
]



def _derive_step2_paths(
    transcripts_tsv: str,
    output_path: Optional[str],
    stats_out: Optional[str],
    input_with_ids_tsv: Optional[str],
) -> tuple[str, str, Optional[str]]:
    txp = Path(transcripts_tsv)
    name = txp.name
    stem = name.replace("_transcripts.tsv", "") if name.endswith("_transcripts.tsv") else txp.stem

    if txp.parent.name == "step1":
        step2_dir = txp.parent.parent / "step2"
        inferred_input = str(txp.parent / f"{stem}_input_with_ids.tsv")
    else:
        step2_dir = txp.parent
        inferred_input = str(txp.parent / f"{stem}_input_with_ids.tsv")

    step2_dir.mkdir(parents=True, exist_ok=True)

    out_default = str(step2_dir / f"{stem}_site_summary.tsv")
    stats_default = str(step2_dir / f"{stem}_step2_stats.tsv")

    return output_path or out_default, stats_out or stats_default, (input_with_ids_tsv or inferred_input)


def _dominant_region(bp: dict[str, int], dominance: str) -> str:
    if dominance == "priority":
        for r in PRIORITY:
            if int(bp.get(r, 0)) > 0:
                return r
        return "INTERGENIC"

    best_r = "INTERGENIC"
    best_v = -1
    for r in PRIORITY:
        v = int(bp.get(r, 0))
        if v > best_v:
            best_v = v
            best_r = r
    return best_r


def _regions_present(bp: dict[str, int]) -> str:
    present = [r for r in PRIORITY if int(bp.get(r, 0)) > 0]
    return "|".join(present) if present else "INTERGENIC"


def _ensure_int(x: Any) -> int:
    try:
        if pd.isna(x):
            return 0
    except Exception:
        pass
    return int(x)


def _maybe_int(row: pd.Series, col: str) -> Optional[int]:
    """
    Return int(row[col]) if the column exists and is not NA; else None.
    Using None keeps blanks in the TSV if the Step1 transcripts.tsv is missing the column.
    """
    if col not in row.index:
        return None
    v = row.get(col)
    try:
        if pd.isna(v):
            return None
    except Exception:
        pass
    try:
        return int(v)
    except Exception:
        return None


def _select_transcript_clash_utr3_first(g: pd.DataFrame) -> pd.Series:
    df = g.copy()
    for c in ["overlap_utr3_bp", "overlap_cds_bp", "overlap_utr5_bp", "overlap_exon_bp", "overlap_tx_bp"]:
        df[c] = df[c].fillna(0).astype(int)

    df["exon_other_bp"] = (
        df["overlap_exon_bp"] - df["overlap_cds_bp"] - df["overlap_utr3_bp"] - df["overlap_utr5_bp"]
    ).clip(lower=0)
    df["intron_bp"] = (df["overlap_tx_bp"] - df["overlap_exon_bp"]).clip(lower=0)

    if "contained_100pct" in df.columns:
        df["contained_100pct"] = df["contained_100pct"].fillna(0).astype(int)
    else:
        df["contained_100pct"] = 0

    df = df.sort_values(
        by=[
            "overlap_utr3_bp",
            "overlap_cds_bp",
            "overlap_utr5_bp",
            "exon_other_bp",
            "intron_bp",
            "overlap_exon_bp",
            "overlap_tx_bp",
            "contained_100pct",
            "transcript_id",
        ],
        ascending=[False, False, False, False, False, False, False, False, True],
        kind="mergesort",
    )
    return df.iloc[0]


def _bp_from_selected_row(row: pd.Series, read_len: int) -> dict[str, int]:
    utr3 = _ensure_int(row.get("overlap_utr3_bp", 0))
    cds = _ensure_int(row.get("overlap_cds_bp", 0))
    utr5 = _ensure_int(row.get("overlap_utr5_bp", 0))
    exon = _ensure_int(row.get("overlap_exon_bp", 0))
    tx = _ensure_int(row.get("overlap_tx_bp", 0))

    exon_other = max(0, exon - cds - utr3 - utr5)
    intron = max(0, tx - exon)

    covered = min(read_len, exon + intron)
    intergenic = max(0, read_len - covered)

    return {"UTR3": utr3, "CDS": cds, "UTR5": utr5, "EXON_OTHER": exon_other, "INTRON": intron, "INTERGENIC": intergenic}


def _bp_from_matrix_rows(g: pd.DataFrame, read_len: int, nt_cols: list[str]) -> dict[str, int]:
    if not nt_cols:
        raise ValueError("Matrix TSV is missing nt_* columns.")
    gg = g.copy()
    gg[nt_cols] = gg[nt_cols].fillna(0).astype(int)

    def row_vec(region: str) -> Optional[pd.Series]:
        r = gg.loc[gg["region"] == region]
        if r.empty:
            return None
        return r[nt_cols].sum(axis=0).clip(upper=1)

    def zeros_vec() -> pd.Series:
        return pd.Series([0] * len(nt_cols), index=nt_cols)

    def vec_or_zero(v: Optional[pd.Series]) -> pd.Series:
        return v if v is not None else zeros_vec()

    utr3_v = vec_or_zero(row_vec("UTR3"))
    cds_v = vec_or_zero(row_vec("CDS"))
    utr5_v = vec_or_zero(row_vec("UTR5"))
    exon_v = row_vec("EXON")
    intron_v = row_vec("INTRON")
    intergenic_v = row_vec("INTERGENIC")

    utr3 = int(utr3_v.sum())
    cds = int(cds_v.sum())
    utr5 = int(utr5_v.sum())

    exon = int(exon_v.sum()) if exon_v is not None else 0
    exon_other = int(((exon_v == 1) & (cds_v == 0) & (utr3_v == 0) & (utr5_v == 0)).sum()) if exon_v is not None else 0
    intron = int(intron_v.sum()) if intron_v is not None else 0
    intergenic = int(intergenic_v.sum()) if intergenic_v is not None else 0

    covered = min(read_len, exon + intron)
    intergenic = max(intergenic, read_len - covered)

    return {"UTR3": utr3, "CDS": cds, "UTR5": utr5, "EXON_OTHER": exon_other, "INTRON": intron, "INTERGENIC": intergenic}


def _precompute_union_bp_from_matrix(mx: pd.DataFrame, read_len_by_id: dict[str, int]) -> dict[str, dict[str, int]]:
    """Precompute UNION bp summaries from the Step 1 matrix once.

    This intentionally delegates to _bp_from_matrix_rows so the calculations stay
    identical to the previous per-read implementation; it only moves the matrix
    scan outside the transcript-selection loop.
    """
    nt_cols = [c for c in mx.columns if c.startswith("nt_")]
    if not nt_cols:
        raise ValueError("Matrix TSV is missing nt_* columns.")

    union_by_id: dict[str, dict[str, int]] = {}
    for rid, g in mx.groupby("id", sort=False):
        rid_str = str(rid)
        read_len = read_len_by_id.get(rid_str)
        if read_len is None:
            if "read_len" not in g.columns or g["read_len"].empty:
                raise ValueError(f"Cannot determine read_len for matrix id={rid_str}.")
            read_len = int(g["read_len"].iloc[0])
        union_by_id[rid_str] = _bp_from_matrix_rows(g, read_len=read_len, nt_cols=nt_cols)
    return union_by_id


def _stats_from_summary(df: pd.DataFrame, top_n: int = 20) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    rows.append({"metric": "step", "value": "step2_transcript_selection"})
    rows.append({"metric": "n_reads", "value": int(len(df))})

    if "ambiguous_union_vs_selected" in df.columns:
        amb = int(df["ambiguous_union_vs_selected"].fillna(0).astype(int).sum())
        rows.append({"metric": "n_ambiguous_union_vs_selected", "value": amb})
        rows.append({"metric": "fraction_ambiguous_union_vs_selected", "value": round(amb / len(df), 6) if len(df) else 0.0})

    for col in ["dominant_region_selected", "dominant_region_union"]:
        if col in df.columns:
            vc = df[col].fillna("NA").astype(str).value_counts()
            for k, v in vc.items():
                rows.append({"metric": f"{col}_count__{k}", "value": int(v)})
                rows.append({"metric": f"{col}_fraction__{k}", "value": round(int(v) / len(df), 6) if len(df) else 0.0})

    for col in ["regions_present_selected", "regions_present_union"]:
        if col in df.columns:
            vc = df[col].fillna("NA").astype(str).value_counts().head(max(1, int(top_n)))
            for i, (k, v) in enumerate(vc.items(), start=1):
                rows.append({"metric": f"top_{col}_{i}", "value": k})
                rows.append({"metric": f"top_{col}_{i}_count", "value": int(v)})

    for col in ["selected_gene_name", "selected_transcript_id"]:
        if col in df.columns:
            vc = df[col].fillna("NA").astype(str).value_counts().head(max(1, int(top_n)))
            for i, (k, v) in enumerate(vc.items(), start=1):
                rows.append({"metric": f"top_{col}_{i}", "value": k})
                rows.append({"metric": f"top_{col}_{i}_count", "value": int(v)})

    return pd.DataFrame(rows)



def _append_rows(path: Path, rows: list[dict[str, Any]], columns: list[str], *, header: bool) -> bool:
    """Append rows to a TSV in bounded chunks, preserving a stable column order."""
    if not rows:
        return header
    mode = "w" if header else "a"
    pd.DataFrame(rows).reindex(columns=columns).to_csv(path, sep="\t", index=False, mode=mode, header=header)
    return False


def _write_empty(path: Path, columns: list[str]) -> None:
    pd.DataFrame(columns=columns).to_csv(path, sep="\t", index=False)


def _make_summary_row(
    *,
    rid: str,
    g: pd.DataFrame,
    union_bp: dict[str, int],
    policy: str,
    dominance: str,
) -> dict[str, Any]:
    """Build one Step 2 summary row using the existing transcript-selection and bp logic."""
    read_len = int(g["read_len"].iloc[0])

    sel = _select_transcript_clash_utr3_first(g)
    sel_bp = _bp_from_selected_row(sel, read_len=read_len)

    dom_sel = _dominant_region(sel_bp, dominance=dominance)
    dom_union = _dominant_region(union_bp, dominance=dominance)

    return {
        "id": rid,
        "read_len": read_len,
        "policy": policy,
        "dominance": dominance,
        "selected_transcript_id": str(sel.get("transcript_id", "")),
        "selected_gene_id": str(sel.get("gene_id", "")),
        "selected_gene_name": str(sel.get("gene_name", "")),
        "selected_tx_start": _maybe_int(sel, "tx_start"),
        "selected_tx_end": _maybe_int(sel, "tx_end"),
        "selected_read_start_in_tx_1based": _maybe_int(sel, "read_start_in_tx_1based"),
        "selected_read_end_in_tx_1based": _maybe_int(sel, "read_end_in_tx_1based"),
        "selected_overlap_start_genome_1based": _maybe_int(sel, "overlap_start_genome_1based"),
        "selected_overlap_end_genome_1based": _maybe_int(sel, "overlap_end_genome_1based"),
        "selected_overlap_start_in_tx_1based": _maybe_int(sel, "overlap_start_in_tx_1based"),
        "selected_overlap_end_in_tx_1based": _maybe_int(sel, "overlap_end_in_tx_1based"),
        "dominant_region_selected": dom_sel,
        "regions_present_selected": _regions_present(sel_bp),
        "bp_utr3_selected": sel_bp["UTR3"],
        "bp_cds_selected": sel_bp["CDS"],
        "bp_utr5_selected": sel_bp["UTR5"],
        "bp_exon_other_selected": sel_bp["EXON_OTHER"],
        "bp_intron_selected": sel_bp["INTRON"],
        "bp_intergenic_selected": sel_bp["INTERGENIC"],
        "dominant_region_union": dom_union,
        "regions_present_union": _regions_present(union_bp),
        "bp_utr3_union": union_bp["UTR3"],
        "bp_cds_union": union_bp["CDS"],
        "bp_utr5_union": union_bp["UTR5"],
        "bp_exon_other_union": union_bp["EXON_OTHER"],
        "bp_intron_union": union_bp["INTRON"],
        "bp_intergenic_union": union_bp["INTERGENIC"],
        "ambiguous_union_vs_selected": 1 if dom_union != dom_sel else 0,
        "n_passing_transcripts": int(len(g)),
    }


def _init_step2_counters() -> tuple[dict[str, Counter], dict[str, int]]:
    counters = {
        "dominant_region_selected": Counter(),
        "dominant_region_union": Counter(),
        "regions_present_selected": Counter(),
        "regions_present_union": Counter(),
        "selected_gene_name": Counter(),
        "selected_transcript_id": Counter(),
    }
    totals = {"n_reads": 0, "ambiguous_union_vs_selected": 0}
    return counters, totals


def _update_step2_counters(row: dict[str, Any], counters: dict[str, Counter], totals: dict[str, int]) -> None:
    totals["n_reads"] += 1
    totals["ambiguous_union_vs_selected"] += int(row.get("ambiguous_union_vs_selected", 0) or 0)

    for col in [
        "dominant_region_selected",
        "dominant_region_union",
        "regions_present_selected",
        "regions_present_union",
        "selected_gene_name",
        "selected_transcript_id",
    ]:
        value = row.get(col)
        counters[col]["NA" if value is None or pd.isna(value) else str(value)] += 1


def _stats_from_counters(counters: dict[str, Counter], totals: dict[str, int], top_n: int = 20) -> pd.DataFrame:
    n_reads = int(totals.get("n_reads", 0))
    amb = int(totals.get("ambiguous_union_vs_selected", 0))

    rows: list[dict[str, Any]] = []
    rows.append({"metric": "step", "value": "step2_transcript_selection"})
    rows.append({"metric": "n_reads", "value": n_reads})
    rows.append({"metric": "n_ambiguous_union_vs_selected", "value": amb})
    rows.append({"metric": "fraction_ambiguous_union_vs_selected", "value": round(amb / n_reads, 6) if n_reads else 0.0})

    for col in ["dominant_region_selected", "dominant_region_union"]:
        for k, v in counters[col].most_common():
            rows.append({"metric": f"{col}_count__{k}", "value": int(v)})
            rows.append({"metric": f"{col}_fraction__{k}", "value": round(int(v) / n_reads, 6) if n_reads else 0.0})

    for col in ["regions_present_selected", "regions_present_union"]:
        for i, (k, v) in enumerate(counters[col].most_common(max(1, int(top_n))), start=1):
            rows.append({"metric": f"top_{col}_{i}", "value": k})
            rows.append({"metric": f"top_{col}_{i}_count", "value": int(v)})

    for col in ["selected_gene_name", "selected_transcript_id"]:
        for i, (k, v) in enumerate(counters[col].most_common(max(1, int(top_n))), start=1):
            rows.append({"metric": f"top_{col}_{i}", "value": k})
            rows.append({"metric": f"top_{col}_{i}_count", "value": int(v)})

    return pd.DataFrame(rows)


def run(
    *,
    transcripts_tsv: str,
    matrix_tsv: str,
    output_tsv: Optional[str] = None,
    stats_out: Optional[str] = None,
    input_with_ids_tsv: Optional[str] = None,
    policy: str = "clash_utr3_first",
    dominance: str = "coverage",
    report: bool = False,
    top_n: int = 20,
) -> None:
    if dominance not in {"coverage", "priority"}:
        raise ValueError("dominance must be one of: coverage, priority")
    if policy != "clash_utr3_first":
        raise ValueError("Currently supported policy: clash_utr3_first")

    out_path, stats_path, input_ids_path = _derive_step2_paths(transcripts_tsv, output_tsv, stats_out, input_with_ids_tsv)
    out_path_obj = Path(out_path)
    final_path_obj = out_path_obj.with_name(out_path_obj.name.replace("_site_summary.tsv", "_final.tsv"))

    log(f"Reading transcripts: {transcripts_tsv}")
    tx = pd.read_csv(transcripts_tsv, sep="\t", dtype={"id": "string"})
    if tx.empty:
        raise ValueError("transcripts_tsv is empty (no passing transcripts).")

    log(f"Reading matrix (UNION composition): {matrix_tsv}")
    mx = pd.read_csv(matrix_tsv, sep="\t", dtype={"id": "string", "region": "string"})
    if mx.empty:
        raise ValueError("matrix_tsv is empty.")

    input_by_id: dict[str, dict[str, Any]] = {}
    original_cols: list[str] = []
    if input_ids_path and Path(input_ids_path).exists():
        log(f"Reading input_with_ids (for original columns): {input_ids_path}")
        input_df = pd.read_csv(input_ids_path, sep="\t", dtype={"id": "string"})
        original_cols = [c for c in input_df.columns if c not in SUMMARY_COLUMNS]
        input_by_id = input_df.set_index("id", drop=False)[original_cols].to_dict(orient="index")
    else:
        if input_ids_path:
            log(f"input_with_ids not found (skipping merge): {input_ids_path}")

    site_columns = ["id"] + [c for c in original_cols if c != "id"] + [c for c in SUMMARY_COLUMNS if c != "id"]
    final_columns = [c for c in FINAL_COLUMNS if c in site_columns]

    log("Precomputing union bp from matrix...")
    read_len_by_id = tx.groupby("id", sort=False)["read_len"].first().astype(int).to_dict()
    union_by_id = _precompute_union_bp_from_matrix(mx, read_len_by_id=read_len_by_id)

    log("Selecting transcript per read and computing site summaries...")
    t0 = time.time()

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)

    out_buffer: list[dict[str, Any]] = []
    final_buffer: list[dict[str, Any]] = []
    out_header = True
    final_header = True
    wrote_out = False
    wrote_final = False
    write_chunk_size = 10_000

    counters, totals = _init_step2_counters()

    for rid, g in tx.groupby("id", sort=False):
        rid_str = str(rid)
        union_bp = union_by_id.get(rid_str)
        if union_bp is None:
            raise ValueError(f"Matrix missing id={rid_str} present in transcripts.")

        row = _make_summary_row(
            rid=rid_str,
            g=g,
            union_bp=union_bp,
            policy=policy,
            dominance=dominance,
        )

        if input_by_id:
            merged = {"id": rid_str}
            merged.update(input_by_id.get(rid_str, {}))
            merged.update(row)
            row = merged

        _update_step2_counters(row, counters, totals)

        out_buffer.append(row)
        final_buffer.append({c: row.get(c) for c in final_columns})

        if len(out_buffer) >= write_chunk_size:
            out_header = _append_rows(out_path_obj, out_buffer, site_columns, header=out_header)
            final_header = _append_rows(final_path_obj, final_buffer, final_columns, header=final_header)
            wrote_out = True
            wrote_final = True
            out_buffer = []
            final_buffer = []

    if out_buffer:
        out_header = _append_rows(out_path_obj, out_buffer, site_columns, header=out_header)
        final_header = _append_rows(final_path_obj, final_buffer, final_columns, header=final_header)
        wrote_out = True
        wrote_final = True

    if not wrote_out:
        _write_empty(out_path_obj, site_columns)
    if not wrote_final:
        _write_empty(final_path_obj, final_columns)

    log(f"Writing: {out_path}")
    log(f"Writing compact final table: {final_path_obj}")

    log(f"Computing stats: {stats_path}")
    _stats_from_counters(counters, totals, top_n=top_n).to_csv(stats_path, sep="\t", index=False)

    n_written = int(totals["n_reads"])
    log(f"Done. Wrote {n_written:,} rows in {time.time()-t0:.1f}s")

    if report:
        log("Report (Step 2)")
        vc = counters["dominant_region_selected"]
        if len(vc):
            log("Top dominant_region_selected:")
            for k, v in vc.most_common(10):
                log(f"  {k}: {v}")
        amb = int(totals["ambiguous_union_vs_selected"])
        log(f"Ambiguous union vs selected: {amb}/{n_written} ({(amb/n_written):.3f})" if n_written else "Ambiguous union vs selected: 0/0 (0.000)")

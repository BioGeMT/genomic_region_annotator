# src/genomic_region_annotator/summarize_sites.py
from __future__ import annotations

import time
from pathlib import Path
from typing import Any, Optional

import pandas as pd


def _ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str) -> None:
    print(f"[{_ts()}] [INFO] {msg}", flush=True)


PRIORITY = ["UTR3", "CDS", "UTR5", "EXON_OTHER", "INTRON", "INTERGENIC"]


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


def _bp_from_matrix_rows(g: pd.DataFrame, read_len: int) -> dict[str, int]:
    nt_cols = [c for c in g.columns if c.startswith("nt_")]
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

    log(f"Reading transcripts: {transcripts_tsv}")
    tx = pd.read_csv(transcripts_tsv, sep="\t", dtype={"id": "string"})
    if tx.empty:
        raise ValueError("transcripts_tsv is empty (no passing transcripts).")

    log(f"Reading matrix (UNION composition): {matrix_tsv}")
    mx = pd.read_csv(matrix_tsv, sep="\t", dtype={"id": "string", "region": "string"})
    if mx.empty:
        raise ValueError("matrix_tsv is empty.")

    input_df: Optional[pd.DataFrame] = None
    if input_ids_path and Path(input_ids_path).exists():
        log(f"Reading input_with_ids (for original columns): {input_ids_path}")
        input_df = pd.read_csv(input_ids_path, sep="\t", dtype={"id": "string"})
    else:
        if input_ids_path:
            log(f"input_with_ids not found (skipping merge): {input_ids_path}")

    log("Selecting transcript per read and computing site summaries...")
    t0 = time.time()

    mx_groups = {k: g for k, g in mx.groupby("id", sort=False)}

    rows: list[dict[str, Any]] = []
    for rid, g in tx.groupby("id", sort=False):
        read_len = int(g["read_len"].iloc[0])

        sel = _select_transcript_clash_utr3_first(g)
        sel_bp = _bp_from_selected_row(sel, read_len=read_len)

        mxg = mx_groups.get(rid)
        if mxg is None:
            raise ValueError(f"Matrix missing id={rid} present in transcripts.")
        union_bp = _bp_from_matrix_rows(mxg, read_len=read_len)

        dom_sel = _dominant_region(sel_bp, dominance=dominance)
        dom_union = _dominant_region(union_bp, dominance=dominance)

        rows.append(
            {
                "id": str(rid),
                "read_len": read_len,
                "policy": policy,
                "dominance": dominance,
                "selected_transcript_id": str(sel.get("transcript_id", "")),
                "selected_gene_id": str(sel.get("gene_id", "")),
                "selected_gene_name": str(sel.get("gene_name", "")),
                # Transcript-relative coordinates for the SELECTED transcript (from Step 1 transcripts.tsv)
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
        )

    out = pd.DataFrame(rows).sort_values("id", kind="mergesort")

    if input_df is not None and not input_df.empty:
        keep = [c for c in input_df.columns if c not in out.columns]
        out = input_df[["id"] + keep].merge(out, on="id", how="right")

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    log(f"Writing: {out_path}")
    out.to_csv(out_path, sep="\t", index=False)
    # Write compact FINAL table (reduced columns)
    final_path = str(Path(out_path).with_name(
        Path(out_path).name.replace("_site_summary.tsv", "_final.tsv")
    ))

    important_cols = [
        "id",
        "chr",
        "start",
        "end",
        "strand",
        "selected_gene_name",
        "selected_transcript_id",
        # relative transcript coordinates
        "selected_read_start_in_tx_1based",
        "selected_read_end_in_tx_1based",
        # region summary
        "dominant_region_selected",
        "regions_present_selected",
        # ambiguity flag
        "ambiguous_union_vs_selected",
    ]

    keep = [c for c in important_cols if c in out.columns]
    final_df = out[keep].copy()

    log(f"Writing compact final table: {final_path}")
    final_df.to_csv(final_path, sep="\t", index=False)

    log(f"Computing stats: {stats_path}")
    _stats_from_summary(out, top_n=top_n).to_csv(stats_path, sep="\t", index=False)

    log(f"Done. Wrote {len(out):,} rows in {time.time()-t0:.1f}s")

    if report:
        log("Report (Step 2)")
        vc = out["dominant_region_selected"].value_counts() if "dominant_region_selected" in out.columns else {}
        if len(vc):
            log("Top dominant_region_selected:")
            for k, v in vc.head(10).items():
                log(f"  {k}: {v}")
        if "ambiguous_union_vs_selected" in out.columns:
            amb = int(out["ambiguous_union_vs_selected"].sum())
            log(f"Ambiguous union vs selected: {amb}/{len(out)} ({(amb/len(out)):.3f})")

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


# Region priority for CLASH-like interpretation
PRIORITY = ["UTR3", "CDS", "UTR5", "EXON_OTHER", "INTRON", "INTERGENIC"]


def _derive_output_path(transcripts_tsv: str, output_path: Optional[str]) -> str:
    if output_path:
        return output_path
    p = Path(transcripts_tsv)
    name = p.name
    if name.endswith("_transcripts.tsv"):
        return str(p.with_name(name.replace("_transcripts.tsv", "_site_summary.tsv")))
    return str(p.with_suffix("")) + "_site_summary.tsv"


def _dominant_region(bp: dict[str, int], dominance: str) -> str:
    if dominance == "priority":
        for r in PRIORITY:
            if int(bp.get(r, 0)) > 0:
                return r
        return "INTERGENIC"
    # coverage-dominant
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


def _select_transcript_clash_utr3_first(g: pd.DataFrame) -> pd.Series:
    """
    CLASH-friendly lexicographic selection:
      max UTR3 bp
      then max CDS bp
      then max UTR5 bp
      then max EXON_OTHER bp (exon - cds - utr3 - utr5)
      then max INTRON bp (tx - exon)
      then max exon bp
      then max tx overlap
    Tie-breakers:
      - prefer contained_100pct == 1 (if present)
      - stable by transcript_id (lexicographic)
    """
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

    return {
        "UTR3": int(utr3),
        "CDS": int(cds),
        "UTR5": int(utr5),
        "EXON_OTHER": int(exon_other),
        "INTRON": int(intron),
        "INTERGENIC": int(intergenic),
    }


def _bp_from_matrix_rows(g: pd.DataFrame, read_len: int) -> dict[str, int]:
    """
    Compute UNION region composition from matrix rows (exact per-nt OR already done in annotate):
      UTR3/CDS/UTR5/INTRON/INTERGENIC: sum of nt_* for that region row
      EXON_OTHER: exon nt_* where cds/utr3/utr5 are 0 (per-nt exact)
    """
    nt_cols = [c for c in g.columns if c.startswith("nt_")]
    if not nt_cols:
        raise ValueError("Matrix TSV is missing nt_* columns.")

    gg = g.copy()
    for c in nt_cols:
        gg[c] = gg[c].fillna(0).astype(int)

    def row_vec(region: str) -> Optional[pd.Series]:
        r = gg.loc[gg["region"] == region]
        if r.empty:
            return None
        # if duplicates exist, OR them (sum then clip)
        return r[nt_cols].sum(axis=0).clip(upper=1)

    def zeros_vec() -> pd.Series:
        return pd.Series([0] * len(nt_cols), index=nt_cols)

    def vec_or_zero(v: Optional[pd.Series]) -> pd.Series:
        return v if v is not None else zeros_vec()

    def row_sum(region: str) -> int:
        v = row_vec(region)
        return int(v.sum()) if v is not None else 0

    utr3_v = vec_or_zero(row_vec("UTR3"))
    cds_v = vec_or_zero(row_vec("CDS"))
    utr5_v = vec_or_zero(row_vec("UTR5"))
    exon_v = row_vec("EXON")
    intron = row_sum("INTRON")
    intergenic = row_sum("INTERGENIC")

    utr3 = int(utr3_v.sum())
    cds = int(cds_v.sum())
    utr5 = int(utr5_v.sum())

    if exon_v is None:
        exon_other = 0
        exon = 0
    else:
        exon = int(exon_v.sum())
        exon_other = int(((exon_v == 1) & (cds_v == 0) & (utr3_v == 0) & (utr5_v == 0)).sum())

    covered = min(read_len, exon + intron)
    intergenic = max(intergenic, read_len - covered)

    return {
        "UTR3": int(utr3),
        "CDS": int(cds),
        "UTR5": int(utr5),
        "EXON_OTHER": int(exon_other),
        "INTRON": int(intron),
        "INTERGENIC": int(intergenic),
    }


def run(
    *,
    transcripts_tsv: str,
    matrix_tsv: str,
    output_tsv: Optional[str] = None,
    policy: str = "clash_utr3_first",
    dominance: str = "coverage",
) -> None:
    """
    Explicit Step 2 (assumption-bearing):

    Inputs:
      - transcripts_tsv: <base>_transcripts.tsv from annotate
      - matrix_tsv: <base>_matrix.tsv from annotate (used for UNION composition)

    Output:
      - <base>_site_summary.tsv (unless output_tsv is provided)

    dominance:
      - 'coverage' (recommended): dominant region = max bp (ties broken by PRIORITY)
      - 'priority': dominant region = first present in PRIORITY
    """
    if dominance not in {"coverage", "priority"}:
        raise ValueError("dominance must be one of: coverage, priority")
    if policy != "clash_utr3_first":
        raise ValueError("Currently supported policy: clash_utr3_first")

    out_path = _derive_output_path(transcripts_tsv, output_tsv)

    log(f"Reading transcripts: {transcripts_tsv}")
    tx = pd.read_csv(
        transcripts_tsv,
        sep="\t",
        dtype={"id": "string", "transcript_id": "string", "gene_id": "string", "gene_name": "string"},
    )
    if tx.empty:
        raise ValueError("transcripts_tsv is empty (no passing transcripts).")

    required_tx = {
        "id",
        "read_len",
        "transcript_id",
        "overlap_tx_bp",
        "overlap_exon_bp",
        "overlap_cds_bp",
        "overlap_utr5_bp",
        "overlap_utr3_bp",
    }
    missing_tx = required_tx - set(tx.columns)
    if missing_tx:
        raise ValueError(f"transcripts_tsv missing required columns: {sorted(missing_tx)}")

    log(f"Reading matrix (UNION composition): {matrix_tsv}")
    mx = pd.read_csv(matrix_tsv, sep="\t", dtype={"id": "string", "region": "string"})
    if mx.empty:
        raise ValueError("matrix_tsv is empty.")
    required_mx = {"id", "region", "read_len"}
    missing_mx = required_mx - set(mx.columns)
    if missing_mx:
        raise ValueError(f"matrix_tsv missing required columns: {sorted(missing_mx)}")

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
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)

    log(f"Writing: {out_path}")
    out.to_csv(out_path, sep="\t", index=False)

    log(f"Done. Wrote {len(out):,} rows in {time.time()-t0:.1f}s")

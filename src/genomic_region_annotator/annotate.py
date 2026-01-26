from __future__ import annotations

from pathlib import Path
import hashlib

import click
import pandas as pd
import pyranges as pr


def _first_existing(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def _norm_chr_series(s: pd.Series) -> pd.Series:
    # strip chr prefix
    x = s.astype(str).str.replace("^chr", "", regex=True)
    # common mito naming
    x = x.replace({"M": "MT", "chrM": "MT", "chrMT": "MT"})
    return x


def _compute_overlap_bp(a_start: pd.Series, a_end: pd.Series, b_start: pd.Series, b_end: pd.Series) -> pd.Series:
    left = pd.concat([a_start, b_start], axis=1).max(axis=1)
    right = pd.concat([a_end, b_end], axis=1).min(axis=1)
    return (right - left).clip(lower=0)


def _is_protein_coding_row(df: pd.DataFrame) -> pd.Series:
    tb = _first_existing(df, ["transcript_biotype", "transcript_type"])
    gb = _first_existing(df, ["gene_biotype", "gene_type"])
    if tb and tb in df.columns:
        return df[tb].astype(str).eq("protein_coding")
    if gb and gb in df.columns:
        return df[gb].astype(str).eq("protein_coding")
    return pd.Series(False, index=df.index)


REGION_PRIORITY: list[str] = [
    "five_prime_utr",
    "three_prime_utr",
    "CDS",
    "start_codon",
    "stop_codon",
    "exon",
    "intron",
    "transcript",
    "gene",
]


def _priority_rank(label: str) -> int:
    if label in REGION_PRIORITY:
        return REGION_PRIORITY.index(label)
    if label == "intergenic":
        return 10_000
    return 9_000


def _sorted_unique_labels(labels: list[str]) -> list[str]:
    return sorted(set(labels), key=_priority_rank)


def _debug_print(
    read_row: pd.Series,
    best_tx_row: pd.Series | None,
    feat_breakdown: pd.DataFrame,
    intron_bp_val: int,
    region_all: str,
    region_dom: str,
    region_primary: str,
    debug_n: int,
) -> None:
    click.echo("\n" + "=" * 80)
    click.echo(f"[DEBUG] Row ID: {int(read_row['_row_id'])}")
    click.echo(
        f"[DEBUG] Input: chr={read_row['chr']} start={read_row['start']} end={read_row['end']} strand={read_row['strand']}"
    )
    click.echo(
        f"[DEBUG] Internal (0-based half-open): Start={int(read_row['Start'])} End={int(read_row['End'])} len={int(read_row['End']-read_row['Start'])}"
    )

    if best_tx_row is None or best_tx_row.empty:
        click.echo("[DEBUG] No transcript selected (intergenic or no overlap).")
    else:
        click.echo("[DEBUG] Selected transcript (coding-prioritized):")
        for k in [
            "transcript_id",
            "transcript_name",
            "transcript_biotype",
            "gene_id",
            "gene_name",
            "gene_biotype",
        ]:
            if k in best_tx_row.index:
                click.echo(f"  - {k}: {best_tx_row[k]}")

    click.echo("\n[DEBUG] Overlap breakdown (Feature -> bp overlap), sorted:")
    if feat_breakdown.empty:
        click.echo("  (no overlaps)")
    else:
        show = feat_breakdown.sort_values(["overlap_bp", "Feature"], ascending=[False, True]).head(debug_n)
        for _, r in show.iterrows():
            click.echo(f"  - {r['Feature']}: {int(r['overlap_bp'])} bp")
        if len(feat_breakdown) > debug_n:
            click.echo(f"  ... ({len(feat_breakdown)-debug_n} more rows not shown)")

    click.echo(f"\n[DEBUG] Derived intron_bp: {int(intron_bp_val)} bp")
    click.echo(f"[DEBUG] region_all: {region_all}")
    click.echo(f"[DEBUG] region_dominant: {region_dom}")
    click.echo(f"[DEBUG] region_primary: {region_primary}")
    click.echo("=" * 80 + "\n")


def _write_stats_tsv(out_df: pd.DataFrame, path: str, top_n: int = 20) -> None:

    rows: list[dict[str, str]] = []

    def add(section: str, key: str, value) -> None:
        rows.append({"section": section, "key": str(key), "value": "" if value is None else str(value)})

    n_total = len(out_df)

    # intergenic vs annotated
    if "region_primary" in out_df.columns:
        n_intergenic = int((out_df["region_primary"] == "intergenic").sum())
    else:
        n_intergenic = 0
    n_annot = n_total - n_intergenic
    pct_intergenic = (n_intergenic / n_total * 100.0) if n_total else 0.0

    add("summary", "n_total", n_total)
    add("summary", "n_annotated", n_annot)
    add("summary", "n_intergenic", n_intergenic)
    add("summary", "pct_intergenic", f"{pct_intergenic:.3f}")

    # transcript assignment
    if "transcript_id" in out_df.columns:
        n_with_tx = int(out_df["transcript_id"].notna().sum())
    else:
        n_with_tx = 0
    add("transcripts", "n_with_transcript", n_with_tx)

    # transcript biotypes
    if "transcript_biotype" in out_df.columns:
        vc_tb = out_df["transcript_biotype"].fillna("NA").value_counts()
        for k, v in vc_tb.items():
            add("transcripts.by_biotype", k, int(v))

        n_pc = int((out_df["transcript_biotype"] == "protein_coding").sum())
        add("transcripts", "n_protein_coding", n_pc)
        add(
            "transcripts",
            "pct_protein_coding_of_assigned",
            f"{(n_pc / n_with_tx * 100.0):.3f}" if n_with_tx else "0.000",
        )

    # region distributions
    for col in ["region_primary", "region_dominant"]:
        if col in out_df.columns:
            vc = out_df[col].fillna("NA").value_counts()
            for k, v in vc.items():
                add(f"{col}.counts", k, int(v))

    # top genes/transcripts (simple frequency)
    if "gene_name" in out_df.columns:
        vc_g = out_df["gene_name"].dropna().value_counts().head(top_n)
        for k, v in vc_g.items():
            add("top_genes.by_gene_name", k, int(v))

    if "gene_id" in out_df.columns:
        vc_gid = out_df["gene_id"].dropna().value_counts().head(top_n)
        for k, v in vc_gid.items():
            add("top_genes.by_gene_id", k, int(v))

    if "transcript_id" in out_df.columns:
        vc_tid = out_df["transcript_id"].dropna().value_counts().head(top_n)
        for k, v in vc_tid.items():
            add("top_transcripts.by_transcript_id", k, int(v))

    if "transcript_name" in out_df.columns:
        vc_tn = out_df["transcript_name"].dropna().value_counts().head(top_n)
        for k, v in vc_tn.items():
            add("top_transcripts.by_transcript_name", k, int(v))

    stats_df = pd.DataFrame(rows)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    stats_df.to_csv(path, sep="\t", index=False)


def _cache_key_for_gtf(gtf_path: str) -> str:
    p = Path(gtf_path)
    st = p.stat()
    h = hashlib.sha256()
    h.update(str(p.resolve()).encode())
    h.update(str(st.st_size).encode())
    h.update(str(int(st.st_mtime)).encode())
    return h.hexdigest()[:16]


def _load_gtf_cached(gtf_path: str, cache_dir: str) -> pr.PyRanges:
    cache_root = Path(cache_dir)
    cache_root.mkdir(parents=True, exist_ok=True)

    key = _cache_key_for_gtf(gtf_path)
    parquet_path = cache_root / f"gtf_cache_{key}.parquet"

    if parquet_path.exists():
        click.echo(f"[INFO] Loading cached GTF parquet: {parquet_path}")
        df = pd.read_parquet(parquet_path)
        return pr.PyRanges(df)

    click.echo("[INFO] Parsing GTF (first run for this file; will be cached)...")
    gtf = pr.read_gtf(gtf_path)
    gtf.Chromosome = _norm_chr_series(gtf.Chromosome)
    gtf_df = gtf.df

    # Keep only columns we actually use (future-proof with fallbacks)
    keep = [
        "Chromosome", "Start", "End", "Strand", "Feature",
        "gene_id", "gene_name", "gene_biotype",
        "transcript_id", "transcript_name", "transcript_biotype",
        "gene_type", "transcript_type",
    ]
    keep = [c for c in keep if c in gtf_df.columns]
    gtf_df = gtf_df[keep].copy()

    # Save parquet
    gtf_df.to_parquet(parquet_path, index=False)
    click.echo(f"[INFO] Cached GTF parquet written: {parquet_path}")

    return pr.PyRanges(gtf_df)


def run(
    input_path: str,
    gtf_path: str,
    output_tsv: str,
    coords: str = "1-based",
    transcript_first: bool = True,
    debug_row_id: int | None = None,
    debug_n: int = 20,
    drop_intergenic: bool = False,
    cache_dir: str = ".cache/genomic-region-annotator",
    report: bool = False,
    stats_out: str | None = None,
    top_n: int = 20,
) -> None:
    click.echo(f"[INFO] Reading input TSV: {input_path}")
    df = pd.read_csv(input_path, sep="\t")
    click.echo(f"[INFO] Rows loaded: {len(df):,}")

    required = {"chr", "start", "end", "strand"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in TSV: {sorted(missing)}")

    # basic validations
    df = df.copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    if (df["end"] < df["start"]).any():
        bad = df[df["end"] < df["start"]].head(5)
        raise ValueError(f"Found rows where end < start (showing up to 5):\n{bad[['chr','start','end','strand']]}")

    df["_row_id"] = range(len(df))
    df["Chromosome"] = _norm_chr_series(df["chr"])
    df["Strand"] = df["strand"].astype(str)
    # accept + / -;
    if (~df["Strand"].isin(["+", "-"])).any():
        click.echo("[WARN] Some strand values are not '+' or '-' (continuing).")

    if coords.lower() == "1-based":
        # 1-based inclusive -> 0-based half-open
        df["Start"] = df["start"] - 1
        df["End"] = df["end"]
    elif coords.lower() == "bed":
        df["Start"] = df["start"]
        df["End"] = df["end"]
    else:
        raise ValueError("coords must be '1-based' or 'bed'")

    reads = pr.PyRanges(df)

    click.echo(f"[INFO] Loading Ensembl GTF: {gtf_path}")
    feats = _load_gtf_cached(gtf_path, cache_dir=cache_dir)
    gtf_df = feats.df
    click.echo(f"[INFO] GTF rows (cached view): {len(gtf_df):,}")

    # Robust column detection
    col_gene_id = _first_existing(gtf_df, ["gene_id"])
    col_gene_name = _first_existing(gtf_df, ["gene_name", "Name", "gene"])
    col_gene_bio = _first_existing(gtf_df, ["gene_biotype", "gene_type"])
    col_tx_id = _first_existing(gtf_df, ["transcript_id"])
    col_tx_name = _first_existing(gtf_df, ["transcript_name", "transcript"])
    col_tx_bio = _first_existing(gtf_df, ["transcript_biotype", "transcript_type"])

    if not col_tx_id:
        raise ValueError("GTF missing transcript_id column; cannot do transcript-first annotation.")

    # Keep transcript-associated features
    feat_df = gtf_df[gtf_df[col_tx_id].notna()].copy()
    feats_tx = pr.PyRanges(feat_df)

    click.echo("[INFO] Joining reads to GTF features...")
    j = reads.join(feats_tx).df

    notes_global: list[str] = []

    if j.empty:
        click.echo("[WARN] No overlaps with any transcript-associated GTF features. All rows will be intergenic.")
        out0 = df.drop(columns=["Chromosome", "Start", "End", "Strand"], errors="ignore").copy()
        out0["region_primary"] = "intergenic"
        out0["region_dominant"] = "intergenic"
        out0["region_all"] = "intergenic"
        out0["gene_id"] = pd.NA
        out0["gene_name"] = pd.NA
        out0["gene_biotype"] = pd.NA
        out0["transcript_id"] = pd.NA
        out0["transcript_name"] = pd.NA
        out0["transcript_biotype"] = pd.NA
        out0["annotation_notes"] = "no_transcript_overlap"
        out0.to_csv(output_tsv, sep="\t", index=False)
        click.echo("[DONE] Finished annotation.")
        return

    if "Start_b" not in j.columns or "End_b" not in j.columns:
        raise RuntimeError("Unexpected PyRanges join output: missing Start_b/End_b columns.")

    j["overlap_bp"] = _compute_overlap_bp(j["Start"], j["End"], j["Start_b"], j["End_b"])
    j = j[j["overlap_bp"] > 0].copy()

    # transcript choice
    j_tx = j.groupby(["_row_id", col_tx_id], as_index=False)["overlap_bp"].sum()

    meta_cols = [col_tx_id]
    for c in [col_tx_name, col_tx_bio, col_gene_id, col_gene_name, col_gene_bio]:
        if c and c in j.columns:
            meta_cols.append(c)

    tx_meta = j[meta_cols].drop_duplicates(subset=[col_tx_id]).copy()
    j_tx = j_tx.merge(tx_meta, on=col_tx_id, how="left")
    j_tx["is_coding"] = _is_protein_coding_row(j_tx).astype(int)

    j_tx = j_tx.sort_values(by=["_row_id", "is_coding", "overlap_bp", col_tx_id],
                            ascending=[True, False, False, True])
    best_tx = j_tx.drop_duplicates(subset=["_row_id"], keep="first").copy()

    out = df.copy()

    out["transcript_id"] = out["_row_id"].map(best_tx.set_index("_row_id")[col_tx_id])
    out["transcript_name"] = out["_row_id"].map(best_tx.set_index("_row_id")[col_tx_name]) if col_tx_name and col_tx_name in best_tx.columns else pd.NA
    out["transcript_biotype"] = out["_row_id"].map(best_tx.set_index("_row_id")[col_tx_bio]) if col_tx_bio and col_tx_bio in best_tx.columns else pd.NA

    out["gene_id"] = out["_row_id"].map(best_tx.set_index("_row_id")[col_gene_id]) if col_gene_id and col_gene_id in best_tx.columns else pd.NA
    out["gene_name"] = out["_row_id"].map(best_tx.set_index("_row_id")[col_gene_name]) if col_gene_name and col_gene_name in best_tx.columns else pd.NA
    out["gene_biotype"] = out["_row_id"].map(best_tx.set_index("_row_id")[col_gene_bio]) if col_gene_bio and col_gene_bio in best_tx.columns else pd.NA

    n_tx = out["transcript_id"].notna().sum()
    click.echo(f"[INFO] Selected a primary transcript for {n_tx:,} / {len(out):,} reads")

    # keep only selected transcript hits
    sel = out[["_row_id", "transcript_id"]].dropna()
    j_sel = j.merge(sel, left_on=["_row_id", col_tx_id], right_on=["_row_id", "transcript_id"], how="inner")

    # region_all labels (all features overlapped)
    region_labels = j_sel.groupby("_row_id")["Feature"].apply(lambda s: _sorted_unique_labels([str(x) for x in s.tolist()]))

    # intron derivation
    exon_df = feat_df[feat_df["Feature"] == "exon"].copy()
    exon_pr = pr.PyRanges(exon_df) if not exon_df.empty else None

    txspan_df = feat_df[feat_df["Feature"] == "transcript"].copy()
    txspan_pr = pr.PyRanges(txspan_df) if not txspan_df.empty else None
    if txspan_pr is None:
        click.echo("[WARN] No 'transcript' features in GTF; intron is estimated from read length - exon overlap.")
        notes_global.append("gtf_missing_transcript_feature")

    out["rel_start"] = pd.NA
    out["rel_end"] = pd.NA

    if txspan_pr is None:
        notes_global.append("no_transcript_span_for_relative_coords")
    else:
        # Build transcript spans per transcript_id (some GTFs may have multiple transcript rows; be safe)
        tx_coords = (
            txspan_df.groupby(col_tx_id, as_index=False)
            .agg(tx_start=("Start", "min"), tx_end=("End", "max"), tx_strand=("Strand", "first"))
        )

        # Join transcript span info onto each read (only reads with transcript_id)
        tx_coords = tx_coords.rename(columns={col_tx_id: "transcript_id"})
        out = out.merge(tx_coords, how="left", on="transcript_id")
        missing_span = out["transcript_id"].notna() & out["tx_start"].isna()
        if missing_span.any():
            click.echo(f"[WARN] {int(missing_span.sum()):,} reads have transcript_id but no transcript span found in GTF.")
            notes_global.append("some_missing_transcript_spans")

        # Compute relative coordinates using internal 0-based half-open [Start, End)
        # Convert to 1-based inclusive on transcript orientation.
        # prefer transcript strand; fall back to read strand
        strand = out["tx_strand"].fillna(out["Strand"]).astype(str)

        plus = strand.eq("+")
        minus = strand.eq("-")

        rel_start_0 = pd.Series(pd.NA, index=out.index, dtype="object")
        rel_end_0 = pd.Series(pd.NA, index=out.index, dtype="object")

        # plus
        rel_start_0.loc[plus] = (out.loc[plus, "Start"] - out.loc[plus, "tx_start"]).astype("Int64")
        rel_end_0.loc[plus] = ((out.loc[plus, "End"] - 1) - out.loc[plus, "tx_start"]).astype("Int64")
        # minus
        rel_start_0.loc[minus] = (out.loc[minus, "tx_end"] - out.loc[minus, "End"]).astype("Int64")
        rel_end_0.loc[minus] = (out.loc[minus, "tx_end"] - out.loc[minus, "Start"] - 1).astype("Int64")

        # to 1-based inclusive
        out["rel_start"] = (rel_start_0.astype("Int64") + 1).astype("Int64")
        out["rel_end"] = (rel_end_0.astype("Int64") + 1).astype("Int64")

        out = out.drop(columns=["tx_start", "tx_end", "tx_strand"], errors="ignore")

    exon_bp = pd.Series(0, index=out["_row_id"], dtype="int64")
    if exon_pr is not None:
        reads_sel_df = out.dropna(subset=["transcript_id"]).copy()
        reads_sel_df[col_tx_id] = reads_sel_df["transcript_id"]
        reads_sel = pr.PyRanges(reads_sel_df)

        j_ex = reads_sel.join(exon_pr).df
        if not j_ex.empty:
            j_ex["overlap_bp"] = _compute_overlap_bp(j_ex["Start"], j_ex["End"], j_ex["Start_b"], j_ex["End_b"])
            j_ex = j_ex[j_ex["overlap_bp"] > 0]
            exon_bp.update(j_ex.groupby("_row_id")["overlap_bp"].sum().astype("int64"))

    read_len = (out["End"] - out["Start"]).astype("int64")
    tx_bp = read_len.copy()

    if txspan_pr is not None:
        reads_sel_df = out.dropna(subset=["transcript_id"]).copy()
        reads_sel_df[col_tx_id] = reads_sel_df["transcript_id"]
        reads_sel = pr.PyRanges(reads_sel_df)

        j_txspan = reads_sel.join(txspan_pr).df
        if not j_txspan.empty:
            j_txspan["overlap_bp"] = _compute_overlap_bp(j_txspan["Start"], j_txspan["End"], j_txspan["Start_b"], j_txspan["End_b"])
            j_txspan = j_txspan[j_txspan["overlap_bp"] > 0]
            tx_span_bp = j_txspan.groupby("_row_id")["overlap_bp"].max().astype("int64")
            tx_bp.update(tx_span_bp)

    intron_bp = (tx_bp - exon_bp).clip(lower=0).astype("int64")
    intron_rows = intron_bp[intron_bp > 0].index.tolist()

    def _labels_for_row(rid: int) -> list[str]:
        base = region_labels.get(rid, [])
        labs = list(base)
        if rid in intron_rows and "intron" not in labs:
            labs.append("intron")
        return _sorted_unique_labels(labs) if labs else ["intergenic"]

    out["region_all"] = out["_row_id"].apply(lambda rid: ",".join(_labels_for_row(rid)))

    # dominant
    INFORMATIVE = {"five_prime_utr", "three_prime_utr", "CDS", "start_codon", "stop_codon", "intron"}

    feat_bp_full = j_sel.groupby(["_row_id", "Feature"], as_index=False)["overlap_bp"].sum()
    if intron_rows:
        feat_bp_full = pd.concat(
            [feat_bp_full,
             pd.DataFrame({"_row_id": intron_rows, "Feature": "intron", "overlap_bp": intron_bp.loc[intron_rows].values})],
            ignore_index=True,
        )

    feat_bp_inf = feat_bp_full[feat_bp_full["Feature"].isin(INFORMATIVE)].copy()
    dom_inf = pd.Series(dtype=object)
    if not feat_bp_inf.empty:
        feat_bp_inf["prio"] = feat_bp_inf["Feature"].map(_priority_rank)
        feat_bp_inf = feat_bp_inf.sort_values(by=["_row_id", "overlap_bp", "prio"], ascending=[True, False, True])
        dom_inf = feat_bp_inf.drop_duplicates(subset=["_row_id"], keep="first").set_index("_row_id")["Feature"]

    feat_bp_fallback = feat_bp_full[~feat_bp_full["Feature"].isin({"transcript", "gene"})].copy()
    feat_bp_fallback["prio"] = feat_bp_fallback["Feature"].map(_priority_rank)
    feat_bp_fallback = feat_bp_fallback.sort_values(by=["_row_id", "overlap_bp", "prio"], ascending=[True, False, True])
    dom_fb = feat_bp_fallback.drop_duplicates(subset=["_row_id"], keep="first").set_index("_row_id")["Feature"]

    out["region_dominant"] = out["_row_id"].map(dom_inf)
    out["region_dominant"] = out["region_dominant"].fillna(out["_row_id"].map(dom_fb)).fillna("intergenic")

    def _primary_from_all(s: str) -> str:
        if not isinstance(s, str) or s.strip() == "" or s == "intergenic":
            return "intergenic"
        labs = _sorted_unique_labels(s.split(","))
        return labs[0] if labs else "intergenic"

    out["region_primary"] = out["region_all"].apply(_primary_from_all)

    # annotation notes (global)
    out["annotation_notes"] = ";".join(notes_global) if notes_global else ""

    # DEBUG (before dropping internal cols)
    if debug_row_id is not None:
        if debug_row_id < 0 or debug_row_id >= len(out):
            raise ValueError(f"--debug-row-id must be between 0 and {len(out)-1}")
        rid = int(debug_row_id)

        read_row = df.loc[df["_row_id"] == rid].iloc[0]
        best_tx_row = out.loc[out["_row_id"] == rid, [
            "transcript_id", "transcript_name", "transcript_biotype",
            "gene_id", "gene_name", "gene_biotype"
        ]].iloc[0]

        fb = feat_bp_full[feat_bp_full["_row_id"] == rid][["Feature", "overlap_bp"]].copy()
        intr_val = int(intron_bp.loc[rid]) if rid in intron_bp.index else 0
        r_all = out.loc[out["_row_id"] == rid, "region_all"].iloc[0]
        r_dom = out.loc[out["_row_id"] == rid, "region_dominant"].iloc[0]
        r_pri = out.loc[out["_row_id"] == rid, "region_primary"].iloc[0]

        _debug_print(read_row, best_tx_row, fb, intr_val, r_all, r_dom, r_pri, debug_n)

    # summary: annotated vs intergenic
    n_intergenic = (out["region_primary"] == "intergenic").sum()
    n_annot = len(out) - n_intergenic
    click.echo(f"[INFO] Summary: {n_annot:,} annotated / {n_intergenic:,} intergenic")

    if report:
        click.echo("[INFO] Top region_primary counts:")
        vc = out["region_primary"].value_counts().head(15)
        for k, v in vc.items():
            click.echo(f"  - {k}: {v:,}")

    if drop_intergenic:
        before = len(out)
        out = out[out["region_primary"] != "intergenic"].copy()
        click.echo(f"[INFO] Dropped intergenic rows: {before - len(out):,} (kept {len(out):,})")

    if stats_out:
        click.echo(f"[INFO] Writing stats TSV: {stats_out}")
        _write_stats_tsv(out, stats_out, top_n=top_n)

    out = out.drop(columns=["Chromosome", "Start", "End", "Strand", "_row_id"], errors="ignore")

    click.echo(f"[INFO] Writing output TSV: {output_tsv}")
    out.to_csv(output_tsv, sep="\t", index=False)
    click.echo("[DONE] Finished annotation.")

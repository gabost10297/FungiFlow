import glob
import os
import subprocess
from io import BytesIO

import pandas as pd
import plotly.express as px
import streamlit as st

BLAST_DIR = "/data/blast_results"
BLAST_LOG = f"{BLAST_DIR}/blast_run.log"
BLAST_DONE_FLAG = f"{BLAST_DIR}/blast_done.flag"
BLAST_SCRIPT = "/data/scripts/run_blast.sh"

PLOTLY_LAYOUT = dict(
    template="plotly_white",
    font=dict(color="#0a160a"),
    colorway=["#1b4332", "#2d6a4f", "#40916c", "#52b788", "#74c69d"],
)

OVERVIEW_CHART_HEIGHT = 400
OVERVIEW_CHART_MARGIN = dict(l=20, r=20, t=60, b=40)


def page_hero(title: str, subtitle: str):
    st.markdown(
        f'<div class="ff-page-hero"><h1>{title}</h1><p>{subtitle}</p></div>',
        unsafe_allow_html=True,
    )



def overview_chart_hint(text: str) -> None:
    """Fixed-height hint above Overview charts so both columns align."""
    st.markdown(
        f'<p class="ff-overview-chart-hint">{text}</p>',
        unsafe_allow_html=True,
    )


def sample_id_from_filename(filename: str) -> str:
    return filename.replace("_blast_summary_strict.tsv", "").replace("_blast_summary.tsv", "")


def list_blast_files(strict_only: bool) -> list[str]:
    if strict_only:
        strict = sorted(glob.glob(os.path.join(BLAST_DIR, "*_blast_summary_strict.tsv")))
        if strict:
            return [os.path.basename(p) for p in strict]
    return sorted(
        os.path.basename(p) for p in glob.glob(os.path.join(BLAST_DIR, "*_blast_summary.tsv"))
    )


def filename_for_sample(sample_id: str, strict: bool, available: list[str]) -> str:
    preferred = (
        f"{sample_id}_blast_summary_strict.tsv"
        if strict
        else f"{sample_id}_blast_summary.tsv"
    )
    if preferred in available:
        return preferred
    return available[0] if available else preferred


def on_strict_toggle_change():
    """Keep same barcode; swap strict <-> full TSV; reset filters tied to old file."""
    strict = st.session_state.get("blast_strict", False)
    available = list_blast_files(strict)
    if not available:
        return

    current = st.session_state.get("blast_file", available[0])
    sid = sample_id_from_filename(current)
    st.session_state["blast_file"] = filename_for_sample(sid, strict, available)

    for key in list(st.session_state.keys()):
        if key.startswith(("blast_conf_", "blast_genus_", "blast_min_pid_")):
            del st.session_state[key]


def sync_blast_file_to_toggle():
    """Ensure session selection matches toggle (handles stale _strict filenames)."""
    strict = st.session_state.get("blast_strict", False)
    available = list_blast_files(strict)
    if not available:
        return

    current = st.session_state.get("blast_file")
    sid = sample_id_from_filename(current) if current else sample_id_from_filename(available[0])
    st.session_state["blast_file"] = filename_for_sample(sid, strict, available)


@st.cache_data
def load_and_clean_data(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")

    def parse_taxonomy(tax_string):
        try:
            parts = str(tax_string).split("|")
            tax_raw = parts[-1] if len(parts) > 1 else str(tax_string)
            tax_dict = {}
            for r in tax_raw.split(";"):
                if "__" in r:
                    lvl, val = r.split("__", 1)
                    tax_dict[lvl] = val
            return pd.Series(tax_dict)
        except Exception:
            return pd.Series({})

    if "Species_Name" in df.columns:
        parsed = df["Species_Name"].apply(parse_taxonomy).rename(
            columns={
                "k": "Kingdom",
                "p": "Phylum",
                "c": "Class",
                "o": "Order",
                "f": "Family",
                "g": "Genus",
                "s": "Species",
            }
        )
        df = pd.concat([df, parsed], axis=1)

    if "Cluster_Name" in df.columns:
        df = df.drop_duplicates(subset=["Cluster_Name"], keep="first")

    for col in ("Percent_Identity", "Alignment_Length", "Query_Coverage(%)", "Query_Length"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


def apply_filters(df: pd.DataFrame, genus: str, min_pident: float, conf_filter: list[str]) -> pd.DataFrame:
    out = df.copy()
    if genus != "(all)" and "Genus" in out.columns:
        out = out[out["Genus"].fillna("").astype(str) == genus]
    if "Percent_Identity" in out.columns and min_pident > 0:
        out = out[out["Percent_Identity"] >= min_pident]
    if conf_filter and "Confidence" in out.columns:
        out = out[out["Confidence"].astype(str).isin(conf_filter)]
    return out


def _hover_columns(df: pd.DataFrame) -> list[str]:
    candidates = [
        "Cluster_Name",
        "Genus",
        "Species",
        "Assigned_Level",
        "Confidence",
        "Query_Length",
        "Query_Coverage(%)",
        "Length_Tier",
    ]
    return [c for c in candidates if c in df.columns]


def _hover_template(hover_cols: list[str], x_line: str, y_line: str | None = None) -> str:
    idx = {name: i for i, name in enumerate(hover_cols)}
    lines = []
    if "Cluster_Name" in idx:
        lines.append(f"<b>Cluster</b>: %{{customdata[{idx['Cluster_Name']}]}}")
    if "Genus" in idx:
        lines.append(f"<b>Genus</b>: %{{customdata[{idx['Genus']}]}}")
    if "Species" in idx:
        lines.append(f"<b>Species</b>: %{{customdata[{idx['Species']}]}}")
    if "Assigned_Level" in idx:
        lines.append(f"Assigned: %{{customdata[{idx['Assigned_Level']}]}}")
    if "Confidence" in idx:
        lines.append(f"Confidence: %{{customdata[{idx['Confidence']}]}}")
    if "Query_Length" in idx:
        lines.append(f"Query len: %{{customdata[{idx['Query_Length']}]}} bp")
    if "Query_Coverage(%)" in idx:
        lines.append(f"Query cov: %{{customdata[{idx['Query_Coverage(%)']}]}}%")
    if "Length_Tier" in idx:
        lines.append(f"Length tier: %{{customdata[{idx['Length_Tier']}]}}")
    lines.append(x_line)
    if y_line:
        lines.append(y_line)
    return "<br>".join(lines) + "<extra></extra>"


def build_genus_bar(df: pd.DataFrame):
    counts = df["Genus"].fillna("Unknown").value_counts().head(12).reset_index()
    counts.columns = ["Genus", "count"]
    n_clusters = len(df)
    subtitle = f"n={n_clusters} clusters · top {len(counts)} genera · hover bars for full name"
    fig = px.bar(counts, x="count", y="Genus", orientation="h")
    fig.update_layout(
        **PLOTLY_LAYOUT,
        height=OVERVIEW_CHART_HEIGHT,
        margin=OVERVIEW_CHART_MARGIN,
        title=dict(text=f"Top genera<br><sup>{subtitle}</sup>", x=0, xanchor="left"),
        xaxis_title="Clusters",
        yaxis=dict(categoryorder="total ascending"),
    )
    return fig


def build_identity_histogram(df: pd.DataFrame, min_pid: float, height: int = OVERVIEW_CHART_HEIGHT):
    s = df["Percent_Identity"].dropna()
    n = len(s)
    subtitle = ""
    if n:
        subtitle = (
            f"n={n} · mean {s.mean():.1f}% · ≥97%: {int((s >= 97).sum())} "
            f"({100 * (s >= 97).sum() / n:.0f}%)"
        )

    fig = px.histogram(
        df,
        x="Percent_Identity",
        nbins=min(25, max(8, n // 4)) if n else 20,
        color_discrete_sequence=["#2d6a4f"],
        opacity=0.75,
    )

    hover_cols = _hover_columns(df)
    if n and hover_cols:
        rug = df[hover_cols + ["Percent_Identity"]]
        fig.add_scatter(
            x=rug["Percent_Identity"],
            y=[0] * n,
            mode="markers",
            marker=dict(
                size=9,
                color="rgba(27,67,50,0.5)",
                line=dict(width=1, color="#1b4332"),
            ),
            customdata=rug[hover_cols].values,
            hovertemplate=_hover_template(hover_cols, "Identity: %{x:.2f}%"),
            name="Consensuses (hover)",
        )
        fig.update_yaxes(visible=False)

    fig.add_vline(x=97, line_dash="dash", line_color="#bc6c25", annotation_text="97%")
    if min_pid != 97:
        fig.add_vline(
            x=min_pid,
            line_dash="dot",
            line_color="#9b2226",
            annotation_text=f"{min_pid:.0f}%",
        )

    title_text = (
        f"Identity distribution<br><sup>{subtitle}</sup>" if subtitle else "Identity distribution"
    )
    fig.update_layout(
        **PLOTLY_LAYOUT,
        height=height,
        margin=OVERVIEW_CHART_MARGIN,
        title=dict(text=title_text, x=0, xanchor="left"),
        xaxis_title="BLAST identity %",
        yaxis_title="Clusters",
    )
    return fig


def build_quality_scatter(df: pd.DataFrame):
    color = "Confidence" if "Confidence" in df.columns else ("Genus" if "Genus" in df.columns else None)
    hover_cols = _hover_columns(df)

    fig = px.scatter(
        df,
        x="Percent_Identity",
        y="Alignment_Length",
        color=color,
        opacity=0.82,
        hover_data=hover_cols,
        title=f"Consensus quality · n={len(df)}",
    )
    fig.update_traces(
        hovertemplate=_hover_template(
            hover_cols,
            "Identity: %{x:.2f}%",
            "Alignment: %{y} bp",
        )
    )
    fig.update_xaxes(autorange="reversed", title="BLAST identity %")
    fig.update_yaxes(title="Alignment length (bp)")
    fig.update_layout(**PLOTLY_LAYOUT, height=420, margin=OVERVIEW_CHART_MARGIN)
    return fig


def render_summary_metrics(df: pd.DataFrame, *, file_label: str, strict_mode: bool):
    n = len(df)
    genera = df["Genus"].nunique() if "Genus" in df.columns else 0
    mean_pid = df["Percent_Identity"].mean() if "Percent_Identity" in df.columns and n else 0.0
    if "Confidence" in df.columns:
        confident = int(df["Confidence"].isin(["high", "medium"]).sum())
    elif "Percent_Identity" in df.columns:
        confident = int((df["Percent_Identity"] >= 97).sum())
    else:
        confident = 0

    mode = "strict TSV" if strict_mode else "full TSV"
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Clusters", n, help=f"Filtered rows · {mode}")
    c2.metric("Genera", genera)
    c3.metric("Mean identity", f"{mean_pid:.1f}%" if n else "—")
    c4.metric("Confident", confident)
    st.caption(f"Source: `{file_label}` · {mode}")


def run_blast_panel():
    st.subheader("Run BLAST pipeline")
    if os.path.exists(BLAST_LOG) and not os.path.exists(BLAST_DONE_FLAG):
        st.warning("BLAST is running…")
        if st.button("Refresh status", key="blast_refresh"):
            st.rerun()
        try:
            with open(BLAST_LOG, encoding="utf-8", errors="replace") as f:
                st.code("".join(f.readlines()[-15:]), language="log")
        except OSError:
            pass
        return

    c1, c2 = st.columns([1, 2])
    with c1:
        if st.button("Run BLAST", type="primary", key="blast_run_btn"):
            os.makedirs(BLAST_DIR, exist_ok=True)
            if os.path.exists(BLAST_DONE_FLAG):
                os.remove(BLAST_DONE_FLAG)
            with open(BLAST_LOG, "w", encoding="utf-8") as log_file:
                subprocess.run(
                    ["bash", BLAST_SCRIPT],
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
            with open(BLAST_DONE_FLAG, "w", encoding="utf-8") as f:
                f.write("done")
            st.cache_data.clear()
            st.rerun()
    with c2:
        st.caption(f"Runs `{BLAST_SCRIPT}` on all samples in `consensus_results/`.")


def show_blast_page():
    page_hero(
        "BLAST Results Explorer",
        "Toggle strict/full TSV per barcode · hover charts for cluster, genus & species",
    )

    paths = list_blast_files(strict_only=False)
    if not paths and not list_blast_files(strict_only=True):
        run_blast_panel()
        st.info("No BLAST tables in `blast_results/`.")
        return

    top1, top2 = st.columns([3, 1])

    with top1:
        st.toggle(
            "Prefer strict TSV",
            value=False,
            key="blast_strict",
            on_change=on_strict_toggle_change,
            help="Uses *_blast_summary_strict.tsv when available. Keeps the same barcode when toggling.",
        )
        sync_blast_file_to_toggle()

        strict = st.session_state["blast_strict"]
        names = list_blast_files(strict)

        selected = st.selectbox(
            "Sample",
            names,
            key="blast_file",
        )

    with top2:
        sid = sample_id_from_filename(selected)
        st.metric("Barcode", sid)
        st.caption("strict" if strict and selected.endswith("_strict.tsv") else "full")

    file_path = os.path.join(BLAST_DIR, selected)
    df_raw = load_and_clean_data(file_path)

    with st.container(border=True):
        st.markdown("**Filters**")
        genera_opts = (
            ["(all)"] + sorted(df_raw["Genus"].dropna().astype(str).unique().tolist())
            if "Genus" in df_raw.columns
            else ["(all)"]
        )
        f1, f2, f3 = st.columns(3)
        with f1:
            genus = st.selectbox("Genus", genera_opts, key=f"blast_genus_{selected}")
        with f2:
            min_pid = st.slider(
                "Min. identity %",
                0.0,
                100.0,
                90.0,
                0.5,
                key=f"blast_min_pid_{selected}",
            )
        with f3:
            if "Confidence" in df_raw.columns:
                opts = sorted(df_raw["Confidence"].dropna().astype(str).unique().tolist())
                conf = st.multiselect(
                    "Confidence",
                    opts,
                    default=opts,
                    key=f"blast_conf_{selected}",
                )
            else:
                conf = []

    df = apply_filters(df_raw, genus, min_pid, conf)
    render_summary_metrics(df, file_label=selected, strict_mode=strict)

    tab_overview, tab_charts, tab_table, tab_run = st.tabs(
        ["Overview", "Charts", "Results table", "Run BLAST"]
    )

    with tab_overview:
        if df.empty:
            st.warning("No rows match filters.")
        else:
            col_genus, col_identity = st.columns(2)
            with col_genus:
                if "Genus" in df.columns:
                    overview_chart_hint(
                        "Bar length = cluster count per genus in the current filter.<br>"
                        "Up to 12 genera by count · hover bars for the full taxon name."
                    )
                    st.plotly_chart(
                        build_genus_bar(df),
                        width="stretch",
                        key=f"genus_bar_{selected}",
                    )
            with col_identity:
                if "Percent_Identity" in df.columns:
                    overview_chart_hint(
                        "Hover green dots under the bars for cluster · genus · species.<br>"
                        "Title shows n, mean identity %, and the share of clusters ≥97%."
                    )
                    st.plotly_chart(
                        build_identity_histogram(df, min_pid, height=OVERVIEW_CHART_HEIGHT),
                        width="stretch",
                        key=f"identity_hist_{selected}",
                    )
            buf = BytesIO()
            df.to_csv(buf, sep="\t", index=False)
            st.download_button(
                "Download filtered TSV",
                buf.getvalue(),
                file_name=selected.replace(".tsv", "_filtered.tsv"),
                mime="text/tab-separated-values",
            )

    with tab_charts:
        if df.empty:
            st.warning("No rows match filters.")
        elif "Percent_Identity" in df.columns and "Alignment_Length" in df.columns:
            st.caption("Hover any point: cluster, genus, species, confidence, query stats")
            st.plotly_chart(build_quality_scatter(df), width="stretch")

    with tab_table:
        st.caption(
            f"All {len(df_raw)} clusters (unfiltered) · scroll the table to browse rows and taxonomy"
        )
        if df_raw.empty:
            st.warning("No rows to display.")
        else:
            cols = [
                c
                for c in [
                    "Cluster_Name",
                    "Query_Length",
                    "Length_Tier",
                    "Genus",
                    "Species",
                    "Percent_Identity",
                    "Query_Coverage(%)",
                    "Alignment_Length",
                    "Confidence",
                    "Assigned_Level",
                    "Species_Name",
                ]
                if c in df_raw.columns
            ]
            view = df_raw[cols] if cols else df_raw
            cfg = {}
            if "Percent_Identity" in view.columns:
                cfg["Percent_Identity"] = st.column_config.ProgressColumn(
                    "Identity %", min_value=0, max_value=100, format="%.1f%%"
                )
            if "Species_Name" in view.columns:
                cfg["Species_Name"] = st.column_config.TextColumn(
                    "Species name",
                    width=480,
                    help="Scroll horizontally in the table to read the full hit string.",
                )
            st.dataframe(
                view,
                column_config=cfg,
                hide_index=True,
                width="stretch",
                height=520,
            )

    with tab_run:
        run_blast_panel()

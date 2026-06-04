import glob
import os
import shutil
import subprocess

import streamlit as st

INTERMEDIATE_DIR = "/data/intermediate_data"
MAFFT_DIR = f"{INTERMEDIATE_DIR}/mafft"
MAFFT_SCRIPT = "/data/scripts/run_mafft.sh"
R_SCRIPT = "/data/scripts/plot_tree.R"
IQTREE_THREADS = os.environ.get("IQTREE_THREADS", "AUTO")


import re

import pandas as pd

from mafft_app import blast_tsv_path, load_blast_table, parse_taxonomy_from_species_name





def tree_blast_context(treefile: str) -> tuple[str, str] | None:
    """Parse sample + strict/full from IQ-TREE output name."""
    base = os.path.basename(treefile)
    match = re.match(
        r"sample_(.+?)_(strict|full)_.+_trimmed\.fasta\.treefile$",
        base,
    )
    if match:
        return match.group(1), match.group(2)
    return None


def _tip_label_text(cluster: str, genus: str, species: str) -> str:
    g = (genus or "").strip() or "Unknown"
    sp = (species or "").strip()
    if sp and sp.lower() not in g.lower():
        return f"{cluster}\n{g} — {sp}"
    return f"{cluster}\n{g}"


def build_tip_metadata_csv(treefile: str) -> str | None:
    """Write tip metadata for R (cluster, genus, species, multiline tip_label)."""
    ctx = tree_blast_context(treefile)
    if not ctx:
        return None
    sample, mode = ctx
    tsv_path = blast_tsv_path(sample, mode)
    if not os.path.isfile(tsv_path):
        return None

    df = load_blast_table(tsv_path, os.path.getmtime(tsv_path))
    if df.empty or "Cluster_Name" not in df.columns:
        return None

    rows = []
    for _, row in df.iterrows():
        cluster = str(row["Cluster_Name"]).strip()
        genus = str(row.get("Genus", "") or "").strip()
        species = str(row.get("Species", "") or "").strip()
        if (not genus or genus == "nan") and "Species_Name" in row:
            parsed = parse_taxonomy_from_species_name(row.get("Species_Name"))
            genus = parsed.get("Genus", genus) or genus
            species = parsed.get("Species", species) or species
        rows.append(
            {
                "cluster": cluster,
                "genus": genus or "Unknown",
                "species": species,
                "tip_label": _tip_label_text(cluster, genus, species),
            }
        )

    if not rows:
        return None

    csv_path = f"{treefile}.tip_meta.csv"
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    return csv_path


def page_hero(title: str, subtitle: str) -> None:
    st.markdown(
        f'<div class="ff-page-hero"><h1>{title}</h1><p>{subtitle}</p></div>',
        unsafe_allow_html=True,
    )


def find_tree_files() -> list[str]:
    patterns = [
        os.path.join(INTERMEDIATE_DIR, "*.treefile"),
        os.path.join(INTERMEDIATE_DIR, "**", "*.treefile"),
    ]
    found: list[str] = []
    for pattern in patterns:
        found.extend(glob.glob(pattern, recursive=True))
    return sorted(set(found))


def tree_asset_paths(treefile: str) -> dict[str, str]:
    base = treefile.replace(".treefile", "")
    return {
        "rect_png": f"{base}_tree_rect.png",
        "circ_png": f"{base}_tree_circ.png",
        "rect_pdf": f"{base}_tree_rect.pdf",
        "circ_pdf": f"{base}_tree_circ.pdf",
        "legacy_png": f"{base}_tree.png",
        "legacy_pdf": f"{base}_tree.pdf",
        "report": f"{base}.iqtree",
        "log": f"{base}.log",
    }


def count_fasta_seqs(path: str) -> int:
    with open(path, encoding="utf-8") as f:
        return sum(1 for line in f if line.startswith(">"))


def render_tree_plot(treefile: str, layout: str) -> str:
    """Run R script; return path to PNG."""
    assets = tree_asset_paths(treefile)
    out_png = assets["circ_png"] if layout == "circular" else assets["rect_png"]
    meta_csv = build_tip_metadata_csv(treefile)
    cmd = ["Rscript", R_SCRIPT, treefile, out_png, layout]
    if meta_csv:
        cmd.append(meta_csv)
    subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True,
    )
    return out_png


def ensure_tree_plots(treefile: str, layouts: tuple[str, ...] = ("rectangular", "circular")) -> None:
    assets = tree_asset_paths(treefile)
    for layout in layouts:
        target = assets["circ_png"] if layout == "circular" else assets["rect_png"]
        if not os.path.isfile(target):
            render_tree_plot(treefile, layout)


def display_tree_image(treefile: str, layout_choice: str) -> None:
    assets = tree_asset_paths(treefile)
    if layout_choice == "Circular":
        img = assets["circ_png"]
        pdf = assets["circ_pdf"]
    else:
        img = assets["rect_png"] if os.path.isfile(assets["rect_png"]) else assets["legacy_png"]
        pdf = assets["rect_pdf"] if os.path.isfile(assets["rect_pdf"]) else assets["legacy_pdf"]

    if not os.path.isfile(img):
        with st.spinner(f"Rendering {layout_choice.lower()} tree…"):
            try:
                render_tree_plot(
                    treefile,
                    "circular" if layout_choice == "Circular" else "rectangular",
                )
            except subprocess.CalledProcessError as exc:
                st.error(f"R/ggtree failed: {exc.stderr or exc}")
                return

    if os.path.isfile(img):
        st.image(img, use_container_width=True)
    else:
        st.warning("Tree image could not be generated.")

    if os.path.isfile(pdf):
        with open(pdf, "rb") as f:
            st.download_button(
                f"Download {layout_choice} PDF",
                f.read(),
                file_name=os.path.basename(pdf),
                mime="application/pdf",
                width="stretch",
            )


def build_iqtree_command(
    alignment: str,
    prefix: str,
    *,
    model: str,
    bootstrap: int,
    alrt: int,
) -> list[str]:
    cmd = [
        "iqtree2" if shutil.which("iqtree2") else "iqtree",
        "-s",
        alignment,
        "-T",
        str(IQTREE_THREADS),
        "-pre",
        prefix,
        "-redo",
    ]
    if model == "MFP":
        cmd.extend(["-m", "MFP"])
    else:
        cmd.extend(["-m", "HKY+F+R5"])

    if bootstrap > 0:
        cmd.extend(["-bb", str(bootstrap), "-bnni"])
    if alrt > 0:
        cmd.extend(["-alrt", str(alrt)])
    return cmd


def iqtree_running(prefix: str) -> bool:
    treefile = f"{prefix}.treefile"
    logfile = f"{prefix}.log"
    return os.path.isfile(logfile) and not os.path.isfile(treefile)


def launch_iqtree(cmd: list[str]) -> None:
    log_path = cmd[cmd.index("-pre") + 1] + ".log"
    with open(log_path, "w", encoding="utf-8") as log:
        subprocess.Popen(
            cmd,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )


def render_explorer() -> str | None:
    trees = find_tree_files()
    if not trees:
        return None

    labels = {os.path.basename(p): p for p in trees}
    selected_name = st.selectbox(
        "Select tree run",
        options=sorted(labels.keys()),
        key="iqtree_select_tree",
    )
    treefile = labels[selected_name]
    if os.path.isfile(treefile):
        try:
            ensure_tree_plots(treefile)
        except subprocess.CalledProcessError:
            pass
    assets = tree_asset_paths(treefile)

    c1, c2, c3 = st.columns(3)
    with open(treefile, encoding="utf-8") as f:
        c1.download_button(
            "Download Newick (.treefile)",
            f.read(),
            file_name=os.path.basename(treefile),
            mime="text/plain",
            width="stretch",
        )

    layout_choice = st.radio(
        "Tree layout",
        ["Rectangular", "Circular"],
        horizontal=True,
        key="iqtree_layout_view",
    )

    if c2.button("Regenerate both layouts", width="stretch"):
        for layout in ("rectangular", "circular"):
            try:
                render_tree_plot(treefile, layout)
                st.success(f"Updated {layout} plot.")
            except subprocess.CalledProcessError as exc:
                st.error(exc.stderr or str(exc))

    with st.container(border=True):
        display_tree_image(treefile, layout_choice)

    st.markdown("---")
    st.subheader("Run reports")
    if os.path.isfile(assets["report"]):
        with st.expander("IQ-TREE report (.iqtree)"):
            with open(assets["report"], encoding="utf-8", errors="replace") as f:
                st.text(f.read())
    if os.path.isfile(assets["log"]):
        with st.expander("Execution log (.log)"):
            with open(assets["log"], encoding="utf-8", errors="replace") as f:
                st.text("".join(f.readlines()[-40:]))

    return treefile


def render_generate_panel(*, has_trees: bool) -> None:
    st.markdown("### Generate new tree")

    analysis_mode = st.radio(
        "Scope",
        ["Single sample (BLAST + MAFFT)", "Legacy global alignment"],
        horizontal=True,
        help="Per-sample trees use BLAST-filtered MAFFT alignments from the Consensuses module.",
    )

    selected_sample = ""
    blast_set = "strict"
    if analysis_mode.startswith("Single"):
        consensus_root = "/data/consensus_results"
        if not os.path.isdir(consensus_root):
            st.warning("No `consensus_results/` folder found.")
            return
        samples = sorted(
            d
            for d in os.listdir(consensus_root)
            if os.path.isdir(os.path.join(consensus_root, d))
        )
        if not samples:
            st.warning("No samples found.")
            return
        selected_sample = st.selectbox("Sample (barcode)", samples, key="iqtree_sample")
        blast_set = st.radio(
            "BLAST set",
            ["strict", "full"],
            format_func=lambda m: "Strict (high + medium)" if m == "strict" else "Full",
            horizontal=True,
            key="iqtree_blast_set",
        )

    run_suffix = (
        st.text_input("Run name", value="default", key="iqtree_run_suffix")
        .strip()
        .replace(" ", "_")
    )

    with st.container(border=True):
        st.markdown("**IQ-TREE settings**")
        o1, o2, o3 = st.columns(3)
        model = o1.selectbox(
            "Substitution model",
            ["MFP", "HKY+F+R5"],
            help="MFP = ModelFinder (recommended for ITS).",
        )
        bootstrap = o2.selectbox(
            "Ultrafast bootstrap (UFBoot)",
            [1000, 500, 200, 0],
            index=0,
            help="Lower values run faster; 0 skips bootstrap.",
        )
        alrt = o3.selectbox("SH-aLRT replicates", [1000, 500, 0], index=0)

    if analysis_mode.startswith("Single"):
        new_prefix = f"{INTERMEDIATE_DIR}/sample_{selected_sample}_{blast_set}_{run_suffix}"
        mafft_raw = f"{MAFFT_DIR}/{selected_sample}_{blast_set}_mafft.fasta"
        mafft_trim = f"{MAFFT_DIR}/{selected_sample}_{blast_set}_mafft_trimmed.fasta"
    else:
        new_prefix = f"{INTERMEDIATE_DIR}/global_{run_suffix}"
        mafft_raw = f"{INTERMEDIATE_DIR}/mafft_alignment.fasta"
        mafft_trim = mafft_raw

    trimmed = f"{new_prefix}_trimmed.fasta"
    iq_prefix = f"{new_prefix}_trimmed.fasta"
    treefile = f"{iq_prefix}.treefile"
    logfile = f"{iq_prefix}.log"

    st.caption(f"Outputs: `{os.path.basename(iq_prefix)}.*`")

    if iqtree_running(iq_prefix):
        st.warning("IQ-TREE is running…")
        if st.button("Refresh", key="iqtree_refresh_run"):
            st.rerun()
        if os.path.isfile(logfile):
            with open(logfile, encoding="utf-8", errors="replace") as f:
                st.code("".join(f.readlines()[-20:]), language="log")
        return

    if os.path.isfile(treefile):
        st.success("This run already finished. Pick it in **Tree explorer** above or delete to re-run.")
        if st.button("Delete this run", type="primary", key="iqtree_delete_run"):
            for pattern in (f"{iq_prefix}*", f"{new_prefix}*"):
                for path in glob.glob(pattern):
                    try:
                        os.remove(path)
                    except OSError:
                        pass
            st.rerun()
        return

    if st.button("Run trimAl + IQ-TREE", type="primary", key="iqtree_run_pipeline"):
        align_input = mafft_trim if os.path.isfile(mafft_trim) else mafft_raw

        if analysis_mode.startswith("Single"):
            if not os.path.isfile(align_input):
                with st.spinner("Running MAFFT…"):
                    subprocess.run(
                        ["bash", MAFFT_SCRIPT, selected_sample, blast_set],
                        check=False,
                    )
                align_input = mafft_trim if os.path.isfile(mafft_trim) else mafft_raw
            if not os.path.isfile(align_input):
                st.error("No MAFFT alignment. Run MAFFT on the Consensuses page first.")
                return
            n = count_fasta_seqs(align_input)
            if n < 3:
                st.error(f"Need ≥3 sequences (found {n}).")
                return
        elif not os.path.isfile(align_input):
            st.error(f"Legacy alignment missing: {align_input}")
            return

        with st.status("Building tree…", expanded=True) as status:
            st.write("trimAl (-gappyout)")
            if align_input.endswith("_trimmed.fasta") or align_input.endswith("_mafft_trimmed.fasta"):
                shutil.copy(align_input, trimmed)
            else:
                subprocess.run(
                    ["trimal", "-in", align_input, "-out", trimmed, "-gappyout"],
                    check=True,
                )

            st.write("IQ-TREE")
            cmd = build_iqtree_command(
                trimmed,
                iq_prefix,
                model=model,
                bootstrap=int(bootstrap),
                alrt=int(alrt),
            )
            status.update(label="IQ-TREE started in background", state="complete")

        launch_iqtree(cmd)
        st.rerun()


def show_iqtree_page() -> None:
    page_hero(
        "Phylogenetic trees (IQ-TREE)",
        "BLAST-aware alignments · optimized IQ-TREE · rectangular and circular ggtree views",
    )

    trees = find_tree_files()
    tab_explore, tab_build = st.tabs(["Tree explorer", "Generate"])

    with tab_explore:
        if not trees:
            st.info("No trees yet. Use **Generate** to run IQ-TREE.")
        else:
            render_explorer()

    with tab_build:
        render_generate_panel(has_trees=bool(trees))

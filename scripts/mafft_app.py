import glob
import os
import re
import subprocess
from collections import Counter

import pandas as pd
import streamlit as st

MAFFT_DIR = "/data/intermediate_data/mafft"
BLAST_DIR = "/data/blast_results"
CONSENSUS_DIR = "/data/consensus_results"
MAFFT_SCRIPT = "/data/scripts/run_mafft.sh"
MANIFEST = f"{MAFFT_DIR}/manifest.tsv"


def page_hero(title: str, subtitle: str):
    st.markdown(
        f'<div class="ff-page-hero"><h1>{title}</h1><p>{subtitle}</p></div>',
        unsafe_allow_html=True,
    )


def list_blast_samples() -> list[str]:
    paths = glob.glob(os.path.join(BLAST_DIR, "*_blast_summary.tsv"))
    samples = []
    for p in paths:
        name = os.path.basename(p)
        if name.endswith("_blast_summary.tsv") and "_strict" not in name:
            samples.append(name.replace("_blast_summary.tsv", ""))
    return sorted(set(samples))


def blast_tsv_path(sample: str, mode: str) -> str:
    if mode == "strict":
        return os.path.join(BLAST_DIR, f"{sample}_blast_summary_strict.tsv")
    return os.path.join(BLAST_DIR, f"{sample}_blast_summary.tsv")


def alignment_paths(sample: str, mode: str) -> dict[str, str]:
    return {
        "raw": os.path.join(MAFFT_DIR, f"{sample}_{mode}_mafft.fasta"),
        "trimmed": os.path.join(MAFFT_DIR, f"{sample}_{mode}_mafft_trimmed.fasta"),
    }


def count_tsv_clusters(tsv_path: str) -> int:
    if not os.path.isfile(tsv_path):
        return 0
    with open(tsv_path, encoding="utf-8") as f:
        return max(0, sum(1 for _ in f) - 1)


def show_mafft_page():
    page_hero(
        "Multiple Sequence Alignment (MAFFT)",
        "Align consensuses chosen by BLAST — strict (high/medium) or full tables per barcode",
    )

    samples = list_blast_samples()
    if not samples:
        st.warning("No BLAST summaries in `blast_results/`. Run BLAST before MAFFT.")
        return

    if "confirmed_seqs" not in st.session_state:
        st.session_state.confirmed_seqs = []
    if "basket" not in st.session_state:
        st.session_state.basket = set()

    top1, top2 = st.columns([2, 1])
    with top1:
        sample = st.selectbox("Sample (barcode)", samples, key="mafft_sample")
        blast_mode = st.radio(
            "BLAST sequence set",
            options=["strict", "full"],
            format_func=lambda m: "Strict (high + medium confidence)"
            if m == "strict"
            else "Full (all BLAST rows)",
            horizontal=True,
            key="mafft_blast_mode",
        )
    with top2:
        tsv = blast_tsv_path(sample, blast_mode)
        n_blast = count_tsv_clusters(tsv)
        st.metric("Clusters in TSV", n_blast)
        st.caption("strict" if blast_mode == "strict" else "full")

    paths = alignment_paths(sample, blast_mode)
    alignment_choice = st.radio(
        "Alignment to view",
        ["Raw (MAFFT)", "Trimmed (trimAl)"],
        horizontal=True,
        key="mafft_view_type",
    )
    fasta_path = paths["raw"] if alignment_choice.startswith("Raw") else paths["trimmed"]

    with st.container(border=True):
        st.markdown("**Run MAFFT**")
        st.caption(
            "Builds one alignment per barcode from clusters listed in the selected BLAST table. "
            "Uses faster MAFFT settings for large sets (>200 / >1000 sequences)."
        )
        c1, c2, c3 = st.columns(3)
        if c1.button("This sample only", type="primary", key="mafft_run_one"):
            subprocess.run(
                ["bash", MAFFT_SCRIPT, sample, blast_mode],
                check=False,
            )
            st.rerun()
        if c2.button("All samples (strict + full)", key="mafft_run_all"):
            subprocess.run(["bash", MAFFT_SCRIPT, "--all", "both"], check=False)
            st.rerun()
        if c3.button("All samples (current set only)", key="mafft_run_all_mode"):
            subprocess.run(
                ["bash", MAFFT_SCRIPT, "--all", blast_mode],
                check=False,
            )
            st.rerun()

    if os.path.isfile(MANIFEST):
        try:
            manifest = pd.read_csv(MANIFEST, sep="\t")
            row = manifest[
                (manifest["sample"] == sample) & (manifest["mode"] == blast_mode)
            ]
            if not row.empty:
                st.success(
                    f"Latest run: **{int(row.iloc[-1]['n_seqs'])}** sequences aligned."
                )
        except Exception:
            pass

    st.markdown("---")
    if os.path.exists(fasta_path):
        st.caption(f"Viewing `{os.path.basename(fasta_path)}`")
        visualize_mafft_alignment(fasta_path)
    else:
        st.info(
            f"No alignment file yet at `{fasta_path}`. "
            f"Run MAFFT above (needs ≥2 FASTAs matching the BLAST table)."
        )
        legacy = "/data/intermediate_data/mafft_alignment.fasta"
        if os.path.exists(legacy):
            with st.expander("Legacy global alignment (pre-BLAST MAFFT)"):
                visualize_mafft_alignment(legacy)


def visualize_mafft_alignment(alignment_fasta_path):
    colors = {
        'A': "#9F1B12", 'T': "#13781b", 'C': "#3c82d8", 'G': "#dea023", 'N': "#7A7A7A"
    }
    
    sequences = {}
    current_seq = ""
    
    with open(alignment_fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_seq = line[1:]
                sequences[current_seq] = ""
            else:
                sequences[current_seq] += line.upper()

    if not sequences:
        st.error("The alignment file is empty.")
        return

    seq_ids = list(sequences.keys())
    aln_length = max(len(seq) for seq in sequences.values())

    st.subheader("Alignment Viewer Controls")
    
    col1, col2 = st.columns([1, 1])
    
    with col1:
        all_barcodes = []
        for seq in seq_ids:
            match = re.search(r'(barcode\d+)', seq.lower())
            if match:
                bc = match.group(1)
                if bc not in all_barcodes:
                    all_barcodes.append(bc)
        
        dominant_clusters = []
        seen_barcodes = set()
        for seq in seq_ids:
            match = re.search(r'(barcode\d+)', seq.lower())
            if match:
                bc = match.group(1)
                if bc not in seen_barcodes:
                    seen_barcodes.add(bc)
                    dominant_clusters.append(seq)

        selected_barcodes = st.multiselect(
            "Filter by Barcode (Optional):",
            options=sorted(all_barcodes)
        )

        if selected_barcodes:
            available_seqs = [seq for seq in seq_ids if any(bc in seq.lower() for bc in selected_barcodes)]
            base_default = available_seqs[:20]
        else:
            available_seqs = seq_ids
            base_default = dominant_clusters[:15] if dominant_clusters else seq_ids[:15]

        if not st.session_state.confirmed_seqs and not st.session_state.basket:
            st.session_state.confirmed_seqs = base_default

        for seq in available_seqs:
            if f"chk_{seq}" not in st.session_state:
                st.session_state[f"chk_{seq}"] = (seq in base_default)

        st.markdown("**Select from list:**")
        def set_all_checkboxes(state):
            for s in available_seqs:
                st.session_state[f"chk_{s}"] = state

        btn_col1, btn_col2 = st.columns(2)
        btn_col1.button("Select All", on_click=set_all_checkboxes, args=(True,), use_container_width=True)
        btn_col2.button("Clear All", on_click=set_all_checkboxes, args=(False,), use_container_width=True)

        current_selections = []
        with st.container(height=250):
            for seq in available_seqs:
                if st.checkbox(seq, key=f"chk_{seq}"):
                    current_selections.append(seq)

        st.markdown("---")

        st.markdown("**Quick Add / Basket:**")
        basket_input = st.text_area(
            "Paste sequence IDs here:", 
            height=68, 
            label_visibility="collapsed", 
            placeholder="Paste copied sequence names here..."
        )
        
        if st.button("➕ Add to Basket", use_container_width=True):
            if basket_input:
                raw_ids = re.split(r'[,\n\s]+', basket_input)
                valid_ids = [s.strip() for s in raw_ids if s.strip() in sequences]
                if valid_ids:
                    for vid in valid_ids:
                        st.session_state.basket.add(vid)
                    st.rerun()
                else:
                    st.error("No valid IDs found. Make sure you copied them correctly.")

        if st.session_state.basket:
            with st.container(border=True):
                st.markdown("<span style='color:#888; font-size: 13px;'>Items currently in Basket:</span>", unsafe_allow_html=True)
                for b_seq in sorted(list(st.session_state.basket)):
                    c1, c2 = st.columns([5, 1])
                    c1.code(b_seq)
                    if c2.button("✖", key=f"del_{b_seq}", help="Remove from basket"):
                        st.session_state.basket.remove(b_seq)
                        st.rerun()
                
                if st.button("Clear Basket", use_container_width=True):
                    st.session_state.basket.clear()
                    st.rerun()

        st.markdown("---")
  
        if st.button("Update Alignment", type="primary", use_container_width=True):
            st.session_state.confirmed_seqs = current_selections
            st.rerun()

    combined_seqs = list(set(st.session_state.confirmed_seqs) | st.session_state.basket)

    with col2:
        start_pos, end_pos = st.slider(
            "4. Select Alignment Region (bp):",
            1, aln_length, (1, min(500, aln_length))
        )
        highlight_mutations = st.checkbox("Highlight mutations only", value=True)
        st.info(f"Displaying: **{len(combined_seqs)}** sequences.")

    selected_seqs = [s for s in combined_seqs if s in sequences]
    
    if not selected_seqs:
        st.warning("Please select sequences from the list or add them to the basket.")
        return

    window_seqs = [sequences[seq_id][start_pos-1:end_pos] for seq_id in selected_seqs]
    consensus_slice = ""
    for i in range(end_pos - start_pos + 1):
        col_bases = [s[i] for s in window_seqs if i < len(s)]
        if not col_bases:
            consensus_slice += "-"
        else:
            consensus_slice += Counter(col_bases).most_common(1)[0][0]

    # RENDERING HTML
    st.markdown("---")
    html_out = "<div class='aln-container'>"
    
    html_out += f"<div class='aln-pos'><strong class='sticky-label'>{'Position':<35}</strong> [{start_pos} ... {end_pos}]</div>"

    html_out += f"<div class='aln-consensus'><strong class='sticky-label'>{'CONSENSUS':<35}</strong> "
    for char in consensus_slice:
        cls = 'base gap' if char == '-' else 'base'
        style = f"style='background-color: {colors.get(char, '#ffffff')}; font-weight: bold;'" if char != '-' else ""
        html_out += f"<span class='{cls}' {style}>{char}</span>"
    html_out += "</div>"

    for seq_id in selected_seqs:
        seq_slice = sequences[seq_id][start_pos-1:end_pos]
        

        html_out += f"<div class='aln-seq'><span class='copyable-seq sticky-label' contenteditable='true' spellcheck='false' title='Click once to select, then Ctrl+C'>{seq_id[:35]:<35}</span> "
        for i, char in enumerate(seq_slice):
            if char == '-':
                html_out += f"<span class='base gap'>-</span>"
            elif highlight_mutations and i < len(consensus_slice) and char == consensus_slice[i]:
                html_out += f"<span class='base match-dot'>.</span>"
            else:
                html_out += f"<span class='base' style='background-color: {colors.get(char, '#ffffff')};'>{char}</span>"
        html_out += "</div>"
    
    html_out += "</div>"
    
    with st.container(border=True):
        st.markdown(html_out, unsafe_allow_html=True)
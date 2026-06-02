import streamlit as st
import os
from collections import Counter
import re

def show_mafft_page():

    alignment_choice = st.radio(
    "Select Alignment Version:",
    ["Raw Alignment (MAFFT)", "Trimmed Alignment (trimAl)"],
    horizontal=True,
    help="Raw shows all sequenced bases and gaps. Trimmed shows only the high-quality columns used for building the phylogenetic tree."
    )
    
    if alignment_choice == "Raw Alignment (MAFFT)":
        fasta_path = "/data/intermediate_data/mafft_alignment.fasta"
    else:
        fasta_path = "/data/intermediate_data/mafft_alignment_trimmed.fasta"
    
    st.title("Multiple Sequence Alignment (MAFFT)")
    st.info("Visual comparison of the generated consensus sequences. Use the list or basket below to select sequences.")
    
    if 'confirmed_seqs' not in st.session_state:
        st.session_state.confirmed_seqs = []
    if 'basket' not in st.session_state:
        st.session_state.basket = set()
    
    
    if os.path.exists(fasta_path):
        visualize_mafft_alignment(fasta_path)
    else:
        st.warning("The file mafft_alignment.fasta could not be found. Please ensure you have run the MAFFT script.")

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
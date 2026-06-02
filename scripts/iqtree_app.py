import streamlit as st
import subprocess
import os
import glob

def show_iqtree_page():
    st.title("Phylogenetic Analysis (IQ-Tree & ggtree)")
    st.info("Build and explore Maximum-Likelihood phylogenetic trees.")

    r_script_path = "/data/scripts/plot_tree.R"
    intermediate_dir = "/data/intermediate_data"

    st.markdown("Tree Explorer & Viewer")
    
    all_tree_files = glob.glob(f"{intermediate_dir}/*.treefile")
    
    if all_tree_files:
        tree_options = {os.path.basename(f): f for f in all_tree_files}
        
        selected_tree_name = st.selectbox(
            "Select a previously generated tree to visualize:", 
            options=sorted(tree_options.keys())
        )
        
        chosen_tree_file = tree_options[selected_tree_name]
        chosen_output_image = chosen_tree_file.replace(".treefile", "_tree.png")
        chosen_output_pdf = chosen_tree_file.replace(".treefile", "_tree.pdf")
        chosen_iqtree_report = chosen_tree_file.replace(".treefile", ".iqtree")
        chosen_log_path = chosen_tree_file.replace(".treefile", ".log")

        c1, c2, c3 = st.columns(3)
        with open(chosen_tree_file, "r") as f:
            c1.download_button("Download Newick (.treefile)", data=f.read(), file_name=selected_tree_name, use_container_width=True)
        
        if os.path.exists(chosen_output_pdf):
            with open(chosen_output_pdf, "rb") as f:
                c2.download_button("Download High-Res PDF", data=f.read(), file_name=selected_tree_name.replace(".treefile", ".pdf"), mime="application/pdf", use_container_width=True)
        else:
            c2.button("PDF not generated yet", disabled=True, use_container_width=True)
            
        if os.path.exists(chosen_output_image):
            with open(chosen_output_image, "rb") as f:
                c3.download_button("Download PNG", data=f.read(), file_name=selected_tree_name.replace(".treefile", ".png"), mime="image/png", use_container_width=True)
        else:
            c3.button("PNG not generated yet", disabled=True, use_container_width=True)

        with st.expander("Toggle to View Phylogenetic Tree Image", expanded=False):
            if not os.path.exists(chosen_output_image):
                with st.spinner("Generating visualization for selected tree (this might take a moment)..."):
                    try:
                        subprocess.run(["Rscript", r_script_path, chosen_tree_file, chosen_output_image], check=True)
                    except subprocess.CalledProcessError:
                        st.error("Failed to render tree image using R.")
            
            if os.path.exists(chosen_output_image):
                st.image(chosen_output_image, use_container_width=True)

        st.markdown("---")
        st.subheader("Selected Tree Reports")
        if os.path.exists(chosen_iqtree_report):
            with st.expander("Show IQ-Tree Statistical Report (.iqtree)"):
                with open(chosen_iqtree_report, "r") as f: st.text(f.read())
        if os.path.exists(chosen_log_path):
            with st.expander("Show Console Execution Log (.log)"):
                with open(chosen_log_path, "r") as f: st.text(f.read())
    else:
        st.info("No computed trees found in the system yet. Use the panel below to generate your first tree.")

    st.markdown("---")

    st.markdown("Generate New Tree")

    analysis_mode = st.radio(
        "Select Analysis Scope for the new run:", 
        ["Global Tree (All Samples)", "Single Sample Tree"], 
        horizontal=True
    )
 
    selected_sample = ""
    if analysis_mode == "Single Sample Tree":
        if os.path.exists("/data/consensus_results"):
            samples = [d for d in os.listdir("/data/consensus_results") if os.path.isdir(os.path.join("/data/consensus_results", d))]
            if samples:
                selected_sample = st.selectbox("Select Sample for local tree:", sorted(samples))
            else:
                st.warning("No samples found.")
                return
        else:
            st.warning("No consensus folder found.")
            return

    st.markdown("**Name your analysis:**")
    run_suffix = st.text_input(
        "Enter unique run name (e.g., run1, HKY_model, trimal_v2):", 
        value="default",
        help="This prefix will prevent overwriting your older trees."
    ).strip().replace(" ", "_") 

    if analysis_mode == "Global Tree (All Samples)":
        new_prefix = f"{intermediate_dir}/global_{run_suffix}"
        mafft_input = f"{intermediate_dir}/mafft_alignment.fasta"
    else:
        new_prefix = f"{intermediate_dir}/sample_{selected_sample}_{run_suffix}"
        combined_fasta = f"{intermediate_dir}/{selected_sample}_combined.fasta"
        mafft_input = f"{intermediate_dir}/{selected_sample}_mafft.fasta"

    new_trimmed_fasta = f"{new_prefix}_trimmed.fasta"
    new_tree_file = f"{new_prefix}_trimmed.fasta.treefile"
    new_log_path = f"{new_prefix}_trimmed.fasta.log"
    new_output_image = f"{new_prefix}_trimmed.fasta_tree.png"

    with st.expander("Launch Computation Panel", expanded=not all_tree_files):
        
        if os.path.exists(new_log_path) and not os.path.exists(new_tree_file):
            st.warning(f"Analysis '{run_suffix}' is currently processing in the background.")
            
            c1, c2 = st.columns(2)
            if c1.button("Refresh Run Status", use_container_width=True): st.rerun()
            if c2.button("Force Abort This Run", type="primary", use_container_width=True):
                if os.path.exists(new_log_path): os.remove(new_log_path)
                st.rerun()
                
            st.markdown("**Live Log Preview for this run:**")
            try:
                with open(new_log_path, "r") as f:
                    st.code("".join(f.readlines()[-15:]))
            except:
                st.write("Reading logs...")

        elif os.path.exists(new_tree_file):
            st.success(f" An analysis named '{run_suffix}' already exists! You can view it in the 'Tree Explorer' dropdown at the top of the page.")
            st.info("If you want to re-run it under the same name, change the name or click below to clear it.")
            if st.button("Delete this specific run to restart"):
                os.remove(new_tree_file)
                if os.path.exists(new_log_path): os.remove(new_log_path)
                st.rerun()

        else:
            st.info(f"New output files will be saved as: `{os.path.basename(new_prefix)}...`")
            
            if st.button("Run Complete Pipeline", type="primary"):
                
                if analysis_mode == "Single Sample Tree":
                    with st.spinner("Preparing local alignment..."):
                        fasta_files = glob.glob(f"/data/consensus_results/{selected_sample}/*.fasta")
                        if len(fasta_files) < 3:
                            st.error("Need at least 3 sequences for a tree!")
                            return
                        with open(combined_fasta, 'w') as outfile:
                            for fname in fasta_files:
                                with open(fname) as infile: outfile.write(infile.read())
                        with open(mafft_input, "w") as mafft_out:
                            subprocess.run(["mafft", "--auto", combined_fasta], stdout=mafft_out, stderr=subprocess.STDOUT, check=True)

                with st.spinner("1/3: Trimming with trimAl..."):
                    trimal_cmd = ["trimal", "-in", mafft_input, "-out", new_trimmed_fasta, "-gappyout"]
                    subprocess.run(trimal_cmd, check=True, capture_output=True)

                iqtree_cmd = [
                    "iqtree", 
                    "-s", new_trimmed_fasta, 
                    "-m", "HKY+F+R5", 
                    "-B", "1000", 
                    "-T", "AUTO", 
                    "-pre", f"{new_prefix}_trimmed.fasta",
                    "-redo"
                ]
                
                subprocess.Popen(iqtree_cmd, start_new_session=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                st.rerun()
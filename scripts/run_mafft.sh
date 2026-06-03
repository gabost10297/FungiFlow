#!/bin/bash
# Align consensus sequences selected by BLAST results (strict and/or full TSV per barcode).
#
# Usage:
#   run_mafft.sh                          # all samples, strict + full
#   run_mafft.sh barcode25                # one sample, strict + full
#   run_mafft.sh barcode25 strict         # one sample, strict only
#   run_mafft.sh --all full               # all samples, full TSV only
#   run_mafft.sh --all both               # explicit default
#
set -euo pipefail

CONSENSUS_DIR="/data/consensus_results"
BLAST_DIR="/data/blast_results"
OUT_DIR="/data/intermediate_data/mafft"
MANIFEST="${OUT_DIR}/manifest.tsv"
THREADS="${MAFFT_THREADS:-${BLAST_THREADS:-4}}"
RUN_TRIMAL="${RUN_TRIMAL:-1}"
MIN_SEQS=2

mkdir -p "$OUT_DIR"

usage() {
    sed -n '2,10p' "$0" | sed 's/^# \?//'
    exit "${1:-0}"
}

run_mafft_align() {
    local n_seqs="$1"
    local input="$2"
    local output="$3"

    if (( n_seqs < MIN_SEQS )); then
        echo "  Skip MAFFT: need at least ${MIN_SEQS} sequences (found ${n_seqs})."
        return 1
    fi

    echo "  MAFFT (${n_seqs} seqs, ${THREADS} threads)..."
    if (( n_seqs <= 200 )); then
        mafft --thread "$THREADS" --auto "$input" > "$output"
    elif (( n_seqs <= 1000 )); then
        mafft --thread "$THREADS" --retree 1 --maxiterate 0 "$input" > "$output"
    else
        mafft --thread "$THREADS" --parttree --retree 1 --maxiterate 0 "$input" > "$output"
    fi
}

build_input_fasta() {
    local sample="$1"
    local tsv="$2"
    local out_fasta="$3"
    local consensus_dir="${CONSENSUS_DIR}/${sample}"
    local cluster fasta_path
    local n_found=0
    local n_missing=0

    : > "$out_fasta"

    if [[ ! -d "$consensus_dir" ]]; then
        echo "  Error: consensus folder not found: ${consensus_dir}" >&2
        return 1
    fi
    if [[ ! -s "$tsv" ]]; then
        echo "  Error: BLAST table missing or empty: ${tsv}" >&2
        return 1
    fi

    while IFS=$'\t' read -r cluster _; do
        [[ "$cluster" == "Cluster_Name" ]] && continue
        [[ -z "$cluster" ]] && continue

        fasta_path="${consensus_dir}/${cluster}.fasta"
        if [[ -f "$fasta_path" ]]; then
            awk -v id="$cluster" '/^>/{print ">" id; next} {print}' "$fasta_path" >> "$out_fasta"
            n_found=$((n_found + 1))
        else
            echo "  Warning: no FASTA for cluster ${cluster}" >&2
            n_missing=$((n_missing + 1))
        fi
    done < <(cut -f1 "$tsv")

    echo "  Clusters: ${n_found} FASTAs written, ${n_missing} missing." >&2
    echo "$n_found"
}

process_sample_mode() {
    local sample="$1"
    local mode="$2"
    local tsv out_raw out_trim input_tmp n_seqs

    case "$mode" in
        strict)
            tsv="${BLAST_DIR}/${sample}_blast_summary_strict.tsv"
            ;;
        full)
            tsv="${BLAST_DIR}/${sample}_blast_summary.tsv"
            ;;
        *)
            echo "Unknown mode: ${mode} (use strict or full)"
            return 1
            ;;
    esac

    if [[ ! -f "$tsv" ]]; then
        echo ">>> ${sample} [${mode}]: no TSV at ${tsv} — skip"
        return 0
    fi

    echo ">>> ${sample} [${mode}] <<<"
    echo "  Source: ${tsv}"

    input_tmp="${OUT_DIR}/.${sample}_${mode}.input.fasta"
    out_raw="${OUT_DIR}/${sample}_${mode}_mafft.fasta"
    out_trim="${OUT_DIR}/${sample}_${mode}_mafft_trimmed.fasta"

    n_seqs=$(build_input_fasta "$sample" "$tsv" "$input_tmp")
    if (( n_seqs < MIN_SEQS )); then
        rm -f "$input_tmp"
        return 0
    fi

    if ! run_mafft_align "$n_seqs" "$input_tmp" "$out_raw"; then
        rm -f "$input_tmp"
        return 0
    fi
    rm -f "$input_tmp"

    if [[ "$RUN_TRIMAL" == "1" ]] && command -v trimal >/dev/null 2>&1; then
        echo "  trimAl (-gappyout)..."
        trimal -in "$out_raw" -out "$out_trim" -gappyout
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample" "$mode" "$n_seqs" "$tsv" "$out_raw" "${out_trim}" >> "$MANIFEST"

    echo "  Saved: ${out_raw}"
    [[ -f "$out_trim" ]] && echo "  Saved: ${out_trim}"
}

ALL_SAMPLES=0
SAMPLE_ARG=""
MODES=("strict" "full")

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help) usage 0 ;;
        --all) ALL_SAMPLES=1; shift ;;
        strict|full|both)
            if [[ "$1" == "both" ]]; then
                MODES=("strict" "full")
            else
                MODES=("$1")
            fi
            shift
            ;;
        barcode*)
            SAMPLE_ARG="$1"
            shift
            ;;
        *)
            echo "Unknown argument: $1"
            usage 1
            ;;
    esac
done

echo "MAFFT (BLAST-selected consensuses)"
echo "====================================================="
echo "Threads: ${THREADS} | trimAl: ${RUN_TRIMAL}"

echo -e "sample\tmode\tn_seqs\tsource_tsv\traw_alignment\ttrimmed_alignment" > "$MANIFEST"

if [[ -n "$SAMPLE_ARG" ]]; then
    samples=("$SAMPLE_ARG")
elif [[ "$ALL_SAMPLES" == 1 ]]; then
    mapfile -t samples < <(
        find "$BLAST_DIR" -maxdepth 1 -name '*_blast_summary.tsv' ! -name '*_strict.tsv' -printf '%f\n' \
            | sed 's/_blast_summary\.tsv$//' | sort -u
    )
else
    mapfile -t samples < <(
        find "$BLAST_DIR" -maxdepth 1 -name '*_blast_summary.tsv' ! -name '*_strict.tsv' -printf '%f\n' \
            | sed 's/_blast_summary\.tsv$//' | sort -u
    )
fi

if [[ ${#samples[@]} -eq 0 ]]; then
    echo "No BLAST summaries found in ${BLAST_DIR}. Run run_blast.sh first."
    exit 1
fi

for sample in "${samples[@]}"; do
    for mode in "${MODES[@]}"; do
        process_sample_mode "$sample" "$mode" || true
    done
done

echo "====================================================="
echo "Done. Manifest: ${MANIFEST}"

#!/usr/bin/env bash
# =============================================================================
# Step 10 -- Plot isoform structures and per-sample usage (Figs. 6 and 7)
# -----------------------------------------------------------------------------
# Manuscript: Fig. 6 (dRNA) and Fig. 7 (cDNA), panels A/C (structures) and
#             B/D (read support per sample)
#
# scripts/bin/plot_isoform_usage_v14.py is a modified version of FLAIR's
# plot_isoform_usage (https://flair.readthedocs.io/en/latest/scripts.html).
# Additions: reference isoforms drawn on top (--ref_isoforms/--ref_on_top),
# sorting by read support, fixed colour palettes (scripts/bin/colorPalettes/,
# one hex colour per line; <N>c_pallette.txt holds N colours), compact row packing.
# Output: <prefix>_isoforms.png and <prefix>_usage.png
#
# The figures are drawn from the final result tables in results/, so this
# step runs without re-doing steps 00-09:
#   bash scripts/10_plot_isoforms.sh [outdir]
# Requires: python3 with matplotlib, seaborn, numpy
# =============================================================================
set -euo pipefail
REPO_DIR="$(cd "$(dirname "$0")/.." && pwd)"
BIN="${REPO_DIR}/scripts/bin"
RES="${REPO_DIR}/results"
OUT="${1:-${REPO_DIR}/figures}"
mkdir -p "${OUT}"

FM=final_models_v2_quant_stringent_checkSplice
RR=renamed.renamed_quant_stringent_checkSplice.renamed_predictProdORF
REF_KLK4="${REPO_DIR}/ref/KLK4_reference_transcripts.chr19-mod2.bed"
REF_FUSION="${REPO_DIR}/ref/KLK4_KLK4T2_KLKP1_reference_transcripts.chr19-mod2.bed"

plot() {  # <bed> <counts> <prefix> <reference bed> <min_reads>
    local n; n=$(($(wc -l < "$2") - 1))           # number of isoforms -> palette size
    local pal="${BIN}/colorPalettes/${n}c_pallette.txt"
    [[ -s "${pal}" ]] || pal="${BIN}/colorPalettes/60c_pallette.txt"
    python3 "${BIN}/plot_isoform_usage_v14.py" "$1" "$2" "${GENE_ID:-ENSG00000167749.12}" \
        -o "$3" --min_reads "$5" --img_resolution 600 --palette "${pal}" \
        --sort_by_reads --ref_isoforms "$4" --ref_on_top
}

# ---- Fig. 6: direct RNA -----------------------------------------------------
d="${RES}/KLK4/dRNA/all"
plot "${d}/ALL.${FM}.${RR}.nofusion.v2.bed" "${d}/ALL.${FM}.counts.renamed.nofusion.v2.tsv" \
     "${OUT}/KLK4_dRNA_all_minReads1" "${REF_KLK4}" 1
d="${RES}/KLK4T2/dRNA/all"
GENE_ID=ENSG00000167749.12 plot "${d}/ALL.${FM}.${RR}.fusion.v2.bed" "${d}/ALL.${FM}.counts.renamed.fusion.v2.tsv" \
     "${OUT}/KLK4T2_dRNA_all_minReads1" "${REF_FUSION}" 1

# ---- Fig. 7: cDNA -----------------------------------------------------------
d="${RES}/KLK4/cDNA/all"
plot "${d}/final_models_v2_quant_stringent_checkSplice.renamed_predictProdORF.nofusion.v2.bed" \
     "${d}/final_models_bamQuant_stringent_checkSplice.counts.renamed.nofusion.v2.tsv" \
     "${OUT}/KLK4_cDNA_all_minReads1" "${REF_KLK4}" 1
d="${RES}/KLK4T2/cDNA/all"
plot "${d}/final_models_v2_quant_stringent_checkSplice.renamed_predictProdORF.fusion.v2.bed" \
     "${d}/final_models_bamQuant_stringent_checkSplice.counts.renamed.fusion.v2.tsv" \
     "${OUT}/KLK4T2_cDNA_all_minReads1" "${REF_FUSION}" 1

ls -1 "${OUT}"

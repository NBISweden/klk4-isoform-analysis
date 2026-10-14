#!/usr/bin/env bash
# =============================================================================
# Step 10 -- Plot isoform structures and per-sample usage (Figs. 6, 7 and S1)
# -----------------------------------------------------------------------------
# Manuscript:
#   Fig. 6   dRNA, 10 most abundant models:  A/B KLK4, C/D KLK4T2-like fusions
#   Fig. 7   cDNA, 10 most abundant models:  A/B KLK4, C/D KLK4T2-like fusions
#            (A/C = <prefix>_isoforms.png, B/D = <prefix>_usage.png)
#   Fig. S1  all final models (structures only):
#            A dRNA KLK4, B dRNA fusions, C cDNA KLK4, D cDNA fusions
#
# scripts/bin/plot_isoform_usage.py is a modified version of FLAIR's
# plot_isoform_usage (https://flair.readthedocs.io/en/latest/scripts.html).
# Additions: reference isoforms drawn on top (--ref_isoforms/--ref_on_top),
# sorting by read support, fixed colour palettes (scripts/bin/colorPalettes/,
# one hex colour per line; <N>c_palette.txt holds N colours), compact row packing.
# Models with a premature termination codon (ORF tag PTC) are drawn hatched.
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
BQ=final_models_bamQuant_stringent_checkSplice
REF_KLK4="${REPO_DIR}/ref/KLK4_reference_transcripts.chr19-mod2.bed"
REF_FUSION="${REPO_DIR}/ref/KLK4_KLK4T2_KLKP1_reference_transcripts.chr19-mod2.bed"

plot() {  # <bed> <counts> <prefix> <reference bed>
    local n; n=$(($(wc -l < "$2") - 1))           # number of isoforms -> palette size
    local pal="${BIN}/colorPalettes/${n}c_palette.txt"
    [[ -s "${pal}" ]] || pal="${BIN}/colorPalettes/60c_palette.txt"
    python3 "${BIN}/plot_isoform_usage.py" "$1" "$2" ENSG00000167749.12 \
        -o "$3" --min_reads 1 --img_resolution 600 --palette "${pal}" \
        --sort_by_reads --ref_isoforms "$4" --ref_on_top
}

# ---- Fig. 6: direct RNA, 10 most abundant models ----------------------------
d="${RES}/KLK4/dRNA/all-top10"
plot "${d}/all-top10.${FM}.${RR}.nofusion.v2.bed" "${d}/all-top10.${FM}.counts.renamed.nofusion.v2.tsv" \
     "${OUT}/Fig6AB_KLK4_dRNA_top10" "${REF_KLK4}"
d="${RES}/KLK4T2/dRNA/all-top10"
plot "${d}/all-top10.${FM}.${RR}.fusion.v2.bed" "${d}/all-top10.${FM}.counts.renamed.fusion.v2.tsv" \
     "${OUT}/Fig6CD_KLK4T2_dRNA_top10" "${REF_FUSION}"

# ---- Fig. 7: cDNA, 10 most abundant models ----------------------------------
d="${RES}/KLK4/cDNA/all-top10"
plot "${d}/all-top10.${BQ}.${RR}.nofusion.v2.bed" "${d}/all-top10.${BQ}.counts_sorted.renamed.renamed.nofusion.v2.tsv" \
     "${OUT}/Fig7AB_KLK4_cDNA_top10" "${REF_KLK4}"
d="${RES}/KLK4T2/cDNA/all-top10"
plot "${d}/all-top10.${BQ}.${RR}.fusion.v2.bed" "${d}/all-top10.${BQ}.counts_sorted.renamed.renamed.fusion.v2.tsv" \
     "${OUT}/Fig7CD_KLK4T2_cDNA_top10" "${REF_FUSION}"

# ---- Fig. S1: all final models (the *_isoforms.png panels) ------------------
d="${RES}/KLK4/dRNA/all"
plot "${d}/ALL.${FM}.${RR}.nofusion.v2.bed" "${d}/ALL.${FM}.counts.renamed.nofusion.v2.tsv" \
     "${OUT}/FigS1A_KLK4_dRNA_all" "${REF_KLK4}"
d="${RES}/KLK4T2/dRNA/all"
plot "${d}/ALL.${FM}.${RR}.fusion.v2.bed" "${d}/ALL.${FM}.counts.renamed.fusion.v2.tsv" \
     "${OUT}/FigS1B_KLK4T2_dRNA_all" "${REF_FUSION}"
d="${RES}/KLK4/cDNA/all"
plot "${d}/${FM}.renamed_predictProdORF.nofusion.v2.bed" "${d}/${BQ}.counts.renamed.nofusion.v2.tsv" \
     "${OUT}/FigS1C_KLK4_cDNA_all" "${REF_KLK4}"
d="${RES}/KLK4T2/cDNA/all"
plot "${d}/${FM}.renamed_predictProdORF.fusion.v2.bed" "${d}/${BQ}.counts.renamed.fusion.v2.tsv" \
     "${OUT}/FigS1D_KLK4T2_cDNA_all" "${REF_FUSION}"

ls -1 "${OUT}"

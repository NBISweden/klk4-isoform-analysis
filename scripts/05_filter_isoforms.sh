#!/usr/bin/env bash
# =============================================================================
# Step 05 -- Remove isoform models outside KLK4 or with very short exons
# -----------------------------------------------------------------------------
# Manuscript: section 2.7 ("Isoform models with exons smaller than 25 bases
#             were filtered out")
#
# Run separately for both model sets from step 04 (corrected / inconsistent):
#   1. filter_isoforms_by_region.sh  keep models with >= 1 exon overlapping the
#                                     gene region (reads were selected by region,
#                                     so some models belong to neighbouring genes)
#   2. filter_isoforms_by_exon_size.sh  drop models containing an exon < 25 bp;
#                                     such models are mostly alignment artefacts
#                                     and cannot be quantified by flair quantify
#
# Input : flair_${TYPE}-${GENE}/${TYPE}-${GENE}_noSE_all_<set>.isoforms.bed
# Output: flair_${TYPE}-${GENE}/${TYPE}-${GENE}_noSE_all_<set>/
#           *.isoforms.matchRegion.kept-exonLen25.bed    (models kept)
#           *.isoforms.matchRegion.disc-exonLen25.bed    (models removed, exon < 25 bp)
#           *.isoforms.filtOutRegion.bed                 (models removed, outside region)
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"
cd "${WORKDIR}/${GENE}/flair_ALL/flair_${TYPE}-${GENE}"

# gene region as BED (0-based start)
awk -v g="${GENE}" 'BEGIN{OFS="\t"} $1==g {print $2,$3-1,$4,$1}' "${GENE_REGIONS}" > "${GENE}_region.bed"

MIN_EXON_LEN=25

for splicing in corrected inconsistent; do
    pre="${TYPE}-${GENE}_noSE_all_${splicing}"
    mkdir -p "${pre}"
    # one folder per model set
    mv "${pre}".* "${pre}"/ 2>/dev/null || true

    bash "${BIN}/filter_isoforms_by_region.sh" \
        "${pre}/${pre}.isoforms.bed" "${GENE}_region.bed" \
        "${pre}/${pre}.isoforms.matchRegion.bed" "${pre}/${pre}.isoforms.filtOutRegion.bed"

    bash "${BIN}/filter_isoforms_by_exon_size.sh" \
        "${pre}/${pre}.isoforms.matchRegion.bed" "${MIN_EXON_LEN}" \
        "${pre}/${pre}.isoforms.matchRegion.kept-exonLen25.bed" \
        "${pre}/${pre}.isoforms.matchRegion.disc-exonLen25.bed"

    wc -l "${pre}/${pre}.isoforms.bed" "${pre}/${pre}.isoforms.matchRegion.kept-exonLen25.bed"
done

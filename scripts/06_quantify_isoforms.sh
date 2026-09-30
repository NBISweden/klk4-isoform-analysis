#!/usr/bin/env bash
# =============================================================================
# Step 06 -- Per-sample read quantification of the filtered isoform models
# -----------------------------------------------------------------------------
# Manuscript: section 2.7 ("Flair quantify was used for isoform read
#             quantification")
#
# Run separately for the corrected and the inconsistent model set.
#   1. collect the reads of the set and split them by sample (sample is part of
#      the read name, see step 03)
#   2. write the flair manifest (sample, condition, batch, fastq)
#   3. flair quantify
#        --stringent     reads must cover >= 25 bp of the first and last exon
#        --check_splice  reads must match the model at every splice site
#        --generate_map  write read -> isoform assignments per sample; these are
#                        the starting point of the QC in step 07
#   4. rank the models by read support and rename them im<rank>-<set>-R<reads>-P<%>
#      (scripts/bin/rank_and_rename_models.py, see docs/isoform_naming.md).
#      Set a (corrected) is numbered first, set b (inconsistent) continues;
#      set b models with < 10 reads are dropped.
#
# Input : step 05 output, merged reads ${TYPE}-${GENE}.fastq (step 04)
# Output: <set>/${TYPE}-${GENE}_noSE_all_<set>.isoforms.matchRegion.kept-exonLen25_quant_stringent_checkSplice.*
#           .counts.tsv                          isoform x sample count matrix
#           .<sample>.<TYPE>.isoform.read.map.txt  read assignment per sample
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"
FLAIR_DIR="${WORKDIR}/${GENE}/flair_ALL/flair_${TYPE}-${GENE}"
ALL_READS="${WORKDIR}/${GENE}/flair_ALL/${TYPE}-${GENE}.fastq"

for splicing in corrected inconsistent; do
    pre="${TYPE}-${GENE}_noSE_all_${splicing}"
    cd "${FLAIR_DIR}/${pre}"

    # --- 1. reads of this set, split by sample --------------------------------
    awk '{print $4}' "${pre}.bed" > "${pre}.readIDs.txt"
    seqkit grep -n -f "${pre}.readIDs.txt" "${ALL_READS}" > "${pre}.fastq"
    for s in ${SAMPLES}; do
        seqkit grep -n -r -p "${TYPE}-${s}-" "${pre}.fastq" > "${pre}_${TYPE}-${s}.fastq"
    done

    # --- 2. flair manifest: sample <tab> condition <tab> batch <tab> reads ---
    for s in ${SAMPLES}; do
        printf "%s\t%s\ts\t%s\n" "${s}" "${TYPE}" "${PWD}/${pre}_${TYPE}-${s}.fastq"
    done > "${TYPE}_quantifyManifest.tsv"
    seqkit stats "${pre}"*.fastq

    # --- 3. quantify the filtered models ------------------------------------
    bed="${pre}.isoforms.matchRegion.kept-exonLen25.bed"
    awk '{print $4}' "${bed}" > "${bed%.bed}.IsoIDs.txt"
    seqkit grep -n -f "${bed%.bed}.IsoIDs.txt" "${pre}.isoforms.fa" > "${bed%.bed}.fa"

    flair quantify -r "${TYPE}_quantifyManifest.tsv" -i "${bed%.bed}.fa" --isoform_bed "${bed}" \
        --stringent --check_splice --generate_map --threads "${THREADS}" \
        --output "${bed%.bed}_quant_stringent_checkSplice"
done

# ---- 4. rank and rename (set a first, set b continues the numbering) -------
start=1
for splicing in corrected inconsistent; do
    pre="${TYPE}-${GENE}_noSE_all_${splicing}"
    cd "${FLAIR_DIR}/${pre}"
    q="${pre}.isoforms.matchRegion.kept-exonLen25_quant_stringent_checkSplice"
    if [[ "${splicing}" == "corrected" ]]; then letter=a; min=0; else letter=b; min=10; fi
    python3 "${BIN}/rank_and_rename_models.py" "${q}.counts.tsv" --set "${letter}" --start "${start}" \
        --min_reads "${min}" --gene_id "${GENE_ID_KLK4}" \
        --fasta "${pre}.isoforms.matchRegion.kept-exonLen25.fa" \
        --bed "${pre}.isoforms.matchRegion.kept-exonLen25.bed"
    start=$((start + $(wc -l < "${q}.counts.replaceName.tsv")))
done

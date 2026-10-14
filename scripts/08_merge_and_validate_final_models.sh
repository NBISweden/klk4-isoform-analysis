#!/usr/bin/env bash
# =============================================================================
# Step 08 -- Merge both model sets, re-quantify and validate the final models
# -----------------------------------------------------------------------------
# Manuscript: sections 2.7 / 3.5
# Final model list: config/${TYPE}_${GENE}_step08_final_models.tsv
# Model descriptions: docs/KLK4_dRNA_final_models.tsv
#
# 8.1 MERGE
#     The polished models of the corrected (a) and inconsistent (b) set were
#     compared in IGV. Identical models were merged and near-duplicates with low
#     read support were removed, leaving 30 dRNA KLK4 models.
#     -> final_models.bed      structure of the models on the modified chr19
#                              (BED12 of the polished model alignment)
#     -> final_models_hap1.fa  haplotype-1 sequence of every model
#                              (= the polished sequence for unphased models)
#     -> final_models_hap2.fa  haplotype-2 sequence (idem)
#
# 8.2 flair quantify of the final models on all KLK4 reads of both sets
#
# 8.3 Alignment-based read counting (stricter than flair):
#     all reads of a sample are aligned to the modified chr19; a read is counted
#     for a model only if it has exactly the model's splice junctions, the same
#     number of exons and at most MAX_5P_CLIP soft-clipped bases at the 5' end
#     (dRNA: 1, cDNA: 5 -- residual adapter sequence in cDNA reads).
#     Models with < 10 supporting reads were removed. flair and alignment-based
#     counts correlate strongly.
#
# 8.4 Final names (dim<rank>-R<reads>, ranked by the flair counts of 8.2;
#     cim<rank>-R<reads>, ranked by the alignment-based counts of 8.3), ORF annotation with flair predictProductivity, and
#     separation of the KLK4T2-like fusion models (dim10, dim18, dim30).
#
# Input : step 07 output, config/${TYPE}_${GENE}_step08_final_models.tsv
# Output: ${FINAL}/final_models_quant_stringent_checkSplice.counts.tsv  (flair counts)
#         ${FINAL}/mapping_bamQuant/<model>.output.filt.{tsv,fastq,bam}  (alignment-based)
#         ${FINAL}/final_models.renamed_predictProdORF.{fusion,nofusion}.bed
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"
FLAIR_DIR="${WORKDIR}/${GENE}/flair_ALL/flair_${TYPE}-${GENE}"
FINAL="${FLAIR_DIR}/final_models"
MODELS="${CONFIG}/${TYPE}_${GENE}_step08_final_models.tsv"
mkdir -p "${FINAL}"
cd "${FINAL}"

# ---- 8.1 merge -------------------------------------------------------------
# columns of ${MODELS}: isoform_id read_set expected_exons model_alignment final_name
: > final_models.bed; : > final_models_hap1.fa; : > final_models_hap2.fa
tail -n +2 "${MODELS}" | while IFS=$'\t' read -r iso set _ aln _; do
    d="${FLAIR_DIR}/${TYPE}-${GENE}_noSE_all_${set}/polishing_quant/${iso}"
    bedtools bamtobed -bed12 -i "${d}/${iso}.${aln}" >> final_models.bed
    if [[ -s "${d}/${iso}.m-hap1.fa" ]]; then
        cat "${d}/${iso}.m-hap1.fa" >> final_models_hap1.fa
        cat "${d}/${iso}.m-hap2.fa" >> final_models_hap2.fa
    else
        cat "${d}/${iso}.m.fa" | tee -a final_models_hap2.fa >> final_models_hap1.fa
    fi
done
# harmonise the sequence names with the BED names
sed -i.bak -E 's/^>([^ ]+).*/>\1/' final_models_hap1.fa final_models_hap2.fa && rm -f ./*.bak

# ---- per-sample reads of both sets + manifest --------------------------------
for s in ${SAMPLES}; do
    cat "${FLAIR_DIR}"/${TYPE}-${GENE}_noSE_all_{corrected,inconsistent}/${TYPE}-${GENE}_noSE_all_*_${TYPE}-${s}.fastq > "${s}.fastq"
    printf "%s\t%s\ts\t%s\n" "${s}" "${TYPE}" "${PWD}/${s}.fastq"
done > quantifyManifest.tsv

# ---- 8.2 flair quantify on the curated models -------------------------------
flair quantify -r quantifyManifest.tsv -i final_models_hap1.fa --isoform_bed final_models.bed \
    --stringent --check_splice --generate_map --threads "${THREADS}" \
    --output final_models_quant_stringent_checkSplice

# ---- 8.3 alignment-based read counting --------------------------------------
if [[ "${TYPE}" == "cDNA" ]]; then MAX_5P_CLIP=5; else MAX_5P_CLIP=1; fi
mkdir -p mapping_bamQuant
for s in ${SAMPLES}; do
    minimap2 -t "${THREADS}" -ax splice -uf -R "@RG\tID:${s}\tSM:${s}" "${REF}" "${s}.fastq" \
        | samtools sort -@ "${THREADS}" --write-index -o "mapping_bamQuant/${s}.bam"
done

tail -n +2 "${MODELS}" | while IFS=$'\t' read -r iso set nex aln _; do
    model_bam="${FLAIR_DIR}/${TYPE}-${GENE}_noSE_all_${set}/polishing_quant/${iso}/${iso}.${aln}"
    for s in ${SAMPLES}; do
        o="mapping_bamQuant/${iso}_${s}"
        python3 "${BIN}/match_splice_sites.py" "${model_bam}" "mapping_bamQuant/${s}.bam" "${o}.output.tsv"
        # keep: same number of exons ($9), 5' soft clip <= MAX_5P_CLIP ($4)
        awk -v e="${nex}" -v m="${MAX_5P_CLIP}" '$9==e && $4<=m {print $5+$6"\t"$0}' "${o}.output.tsv" \
            | sort -k1,1n -k9,9rn -k4,4rn > "${o}.output.filt.tsv"
        cut -f2 "${o}.output.filt.tsv" > "${o}.output.filt.readIDs.txt"
        seqkit grep -n -f "${o}.output.filt.readIDs.txt" "${s}.fastq" > "${o}.output.filt.fastq"
        minimap2 -t "${THREADS}" -ax splice -uf -R "@RG\tID:${iso}_${s}\tSM:${iso}" "${REF}" "${o}.output.filt.fastq" \
            | samtools sort -@ "${THREADS}" --write-index -o "${o}.output.filt.bam"
        printf "%s\t%s\t%s\n" "${iso}" "${s}" "$(wc -l < "${o}.output.filt.readIDs.txt")"
    done
    samtools merge -f -@ "${THREADS}" --write-index "mapping_bamQuant/${iso}.output.filt.bam" \
        mapping_bamQuant/"${iso}"_*.output.filt.bam && rm mapping_bamQuant/"${iso}"_*.output.filt.bam*
done > mapping_bamQuant/alignment_based_counts.long.tsv

# ---- 8.4 final names, ORF annotation, fusion / nofusion ----------------------
# columns of ${MODELS}: isoform_id read_set expected_exons model_alignment final_name
tail -n +2 "${MODELS}" | awk -F'\t' -v g="${GENE_ID_KLK4}" '{print $1"\t"$5"_"g}' > final_names.tsv
awk -F'\t' 'NR==FNR {n[$1]=$2; next} ($4 in n) {$4=n[$4]} 1' OFS='\t' final_names.tsv final_models.bed \
    > final_models.renamed.bed
# per-sample table of the alignment-based counts (8.3)
awk -F'\t' -v samples="${SAMPLES}" 'BEGIN {n = split(samples, s, " "); printf "id"; for (i = 1; i <= n; i++) printf "\t%s", s[i]; print ""}
    {c[$1, $2] = $3; if (!($1 in seen)) {seen[$1] = 1; order[++k] = $1}}
    END {for (j = 1; j <= k; j++) {printf "%s", order[j]; for (i = 1; i <= n; i++) printf "\t%d", c[order[j], s[i]]; print ""}}' \
    mapping_bamQuant/alignment_based_counts.long.tsv > final_models_bamQuant.counts.tsv

# renamed count tables: dRNA results use the flair counts, cDNA results the alignment-based counts
for c in final_models_quant_stringent_checkSplice final_models_bamQuant; do
    awk -F'\t' 'NR==FNR {n[$1]=$2; next} FNR==1 || ($1 in n) {if (FNR>1) $1=n[$1]; print}' OFS='\t' \
        final_names.tsv "${c}.counts.tsv" > "${c}.counts.renamed.tsv"
done

# ORF prediction; the name gets the tag PRO / PTC / NGO / NST (docs/isoform_naming.md)
predictProductivity -i final_models.renamed.bed -f "${REF}" -g "${GTF_CHR19}" \
    -o final_models.renamed_predictProdORF --longestORF

# KLK4T2-like fusion models (spanning into KLKP1, start < 50,900,000) vs. KLK4 models
awk -F'\t' '$2 <  50900000' final_models.renamed_predictProdORF.bed > final_models.renamed_predictProdORF.fusion.bed
awk -F'\t' '$2 >= 50900000' final_models.renamed_predictProdORF.bed > final_models.renamed_predictProdORF.nofusion.bed

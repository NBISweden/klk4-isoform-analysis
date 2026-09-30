#!/usr/bin/env bash
# =============================================================================
# Step 09 -- Haplotype-resolved (phase-aware) read counts of the final models
# -----------------------------------------------------------------------------
# Manuscript: section 2.7
#
# Why: 12 of the 30 dRNA KLK4 models carry heterozygous SNVs (two haplotypes,
#      see config/dRNA_KLK4_step09_phased_models.tsv). To see whether both
#      alleles are expressed in each sample, the reads that passed the
#      alignment-based filter (step 08.3) are aligned to the haplotype-1
#      sequence of their model and SNVs are called per sample (one read group
#      per sample). For every sample the number of reads supporting the
#      reference (= hap1) and the alternative (= hap2) allele is taken as the
#      MINIMUM over all SNV sites of the model (conservative: a read counts only
#      if every site is covered).
#
#   1. write one FASTA per model (and per haplotype if the two differ)
#   2. per phased model: align reads (per sample), freebayes, collect RO/AO
#
# Input : step 08 output, renamed final sequences
#         final_models_hap{1,2}_quant_stringent_checkSplice.renamed.fa (names im1-R2252, ...)
# Output: ${FINAL}/final_fa_set/<model>[-hap1|-hap2].fa
#         ${FINAL}/phasing_results.tsv
#           id | hap1 reads per sample (8) | sum | hap2 reads per sample (8) | sum | total
#         (the study's result is in results/KLK4/dRNA/phasing/phasing_results.v2.tsv)
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"
FINAL="${WORKDIR}/${GENE}/flair_ALL/flair_${TYPE}-${GENE}/final_models"
PHASED="${CONFIG}/${TYPE}_${GENE}_step09_phased_models.tsv"
cd "${FINAL}"

# ---- 1. one FASTA per model / haplotype -------------------------------------
for h in 1 2; do
    sed -e "s:_${GENE_ID_KLK4}::" "final_models_hap${h}_quant_stringent_checkSplice.renamed.fa" > "hap${h}.fa"
    seqkit split -i -O "hap${h}.split" "hap${h}.fa"
done
mkdir -p final_fa_set
for x in $(grep '>' hap1.fa | cut -c2- | cut -d' ' -f1); do
    h1="hap1.split/hap1.part_${x}.fa"; h2="hap2.split/hap2.part_${x}.fa"
    if [[ "$(seqkit seq -s "${h1}")" == "$(seqkit seq -s "${h2}")" ]]; then
        ln -sfr "${h1}" "final_fa_set/${x}.fa"          # single haplotype
    else
        ln -sfr "${h1}" "final_fa_set/${x}-hap1.fa"     # two haplotypes
        ln -sfr "${h2}" "final_fa_set/${x}-hap2.fa"
    fi
done

# ---- 2. per-sample allele counts for every phased model ---------------------
printf "id" > phasing_results.tsv
for s in ${SAMPLES}; do printf "\thap1_%s" "${s}"; done >> phasing_results.tsv
printf "\thap1_sum" >> phasing_results.tsv
for s in ${SAMPLES}; do printf "\thap2_%s" "${s}"; done >> phasing_results.tsv
printf "\thap2_sum\ttotal\n" >> phasing_results.tsv

# columns of ${PHASED}: isoform_id phased_model_name final_name freebayes_min_qual
tail -n +2 "${PHASED}" | while IFS=$'\t' read -r old new final minq; do
    ref="final_fa_set/${new}-hap1.fa"
    d="phasing_${new}"; mkdir -p "${d}"
    for s in ${SAMPLES}; do
        minimap2 -x map-ont -a -t "${THREADS}" -R "@RG\tID:${new}_${s}\tSM:${new}_${s}" "${ref}" \
            "mapping_bamQuant/${old}_${s}.output.filt.fastq" \
            | samtools sort -@ "${THREADS}" --write-index -o "${d}/${new}.reads.${s}.bam"
    done
    samtools merge -f --write-index "${d}/${new}.reads.bam" "${d}/${new}".reads.*-*.bam
    rm "${d}/${new}".reads.*-*.bam*

    samtools faidx "${ref}"
    # QUAL > 300 for the well supported models, > 10 for the models with < 40 reads
    freebayes -f "${ref}" "${d}/${new}.reads.bam" | vcffilter -f "QUAL > ${minq}" \
        | bgzip > "${d}/${new}.freebayes.vcf.gz"
    tabix -f -p vcf "${d}/${new}.freebayes.vcf.gz"

    # RO = reads supporting the reference (hap1) allele, AO = alternative (hap2)
    slist=$(for s in ${SAMPLES}; do printf "%s_%s," "${new}" "${s}"; done); slist=${slist%,}
    bcftools query -s "${slist}" -f '[%RO\t%AO\t]\n' "${d}/${new}.freebayes.vcf.gz" \
    | awk -v id="${final}" -v n=8 '
        { for (i = 1; i <= n; i++) {
              ro = $(2*i-1); ao = $(2*i)
              if (NR == 1 || ro < RO[i]) RO[i] = ro
              if (NR == 1 || ao < AO[i]) AO[i] = ao } }
        END { line = id; rs = as = 0
              for (i = 1; i <= n; i++) { line = line "\t" RO[i]; rs += RO[i] }
              line = line "\t" rs
              for (i = 1; i <= n; i++) { line = line "\t" AO[i]; as += AO[i] }
              print line "\t" as "\t" rs + as }' >> phasing_results.tsv
done

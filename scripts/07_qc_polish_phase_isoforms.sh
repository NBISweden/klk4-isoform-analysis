#!/usr/bin/env bash
# =============================================================================
# Step 07 -- Per-isoform QC, error correction (polishing) and phasing
# -----------------------------------------------------------------------------
# Manuscript: section 2.7
#
# Why: flair models are consensus paths on the reference; their sequence is the
#      reference sequence, not the transcript sequence of the cell lines, and
#      some are artefacts. For every model we therefore
#        QC     : re-align its reads and keep only reads with exactly the model's
#                 splice junctions (scripts/bin/match_splice_sites.py)
#        polish : pick the best read as template and polish it with all other
#                 supporting reads (racon)  -> read-based model sequence *.m.fa
#        phase  : call heterozygous SNVs on the polished model (freebayes),
#                 phase them (whatshap) and build one sequence per haplotype
#                 (bcftools consensus)      -> *.m-hap1.fa / *.m-hap2.fa
#      Models with < 10 supporting reads after QC were discarded.
#
# MANUAL DECISIONS
#   The template read and the expected exon number are chosen for each model
#   from the candidate list (mode "candidates" below). Criteria
#   for the template read, in this order:
#     (1) fewest alignment errors (insertions + deletions),
#     (2) longest soft-clipped poly(A) tail (3' end = left side for KLK4, - strand),
#     (3) no soft clipping at the 5' end (right side).
#   The parameters used are recorded in
#     config/dRNA_KLK4_step06_polishing_parameters.tsv
#     config/cDNA_KLK4_step06_polishing_parameters.tsv
#   (read set, model, expected exons, template read, freebayes QUAL threshold,
#    phasing sample, outcome).
#
# USAGE
#   07_qc_polish_phase_isoforms.sh remap      <set>                 # once per set
#   07_qc_polish_phase_isoforms.sh candidates <set> <model> <n_exons>
#   07_qc_polish_phase_isoforms.sh polish     <set> <model> <n_exons> <template_read> <min_qual> [phasing_sample]
#   07_qc_polish_phase_isoforms.sh all        # run 'polish' for every row of the parameter table
#     <set> = corrected | inconsistent
#
# Output (per model): <set>/polishing_quant/<model>/
#   <model>.output.tsv          read QC table (match_splice_sites.py)
#   <model>.m.fa                polished model sequence
#   <model>.m-spliced.bam       polished model aligned to the modified chr19
#   <model>.m.freebayes.vcf.gz  SNVs (QUAL filtered)
#   <model>.m-hap{1,2}.fa       haplotype sequences (= polished sequence if not phased)
#   <model>.m-spliced-hap{1,2}.bam
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"
FLAIR_DIR="${WORKDIR}/${GENE}/flair_ALL/flair_${TYPE}-${GENE}"
PARAMS="${CONFIG}/${TYPE}_${GENE}_step06_polishing_parameters.tsv"

mode="${1:?mode missing}"; shift || true

setdir() {                         # cd into the folder of a model set
    pre="${TYPE}-${GENE}_noSE_all_$1"
    cd "${FLAIR_DIR}/${pre}"
    q="${pre}.isoforms.matchRegion.kept-exonLen25_quant_stringent_checkSplice"
    models_fa="${q}.renamed.fa"             # renamed flair models (see step 06 NEXT)
    replace="${q}.counts.replaceName.tsv"   # old flair name -> im* name
    reads="${pre}.fastq"
}

# -----------------------------------------------------------------------------
# remap: align the reads assigned to each model by flair quantify (per sample)
#        to the modified chr19 and merge them into one BAM per model
# -----------------------------------------------------------------------------
remap() {
    setdir "$1"
    mkdir -p mapping_quant
    for s in ${SAMPLES}; do
        while IFS=$'\t' read -r iso iso_reads; do
            grep -P "\t${iso}\t" "${pre}.isoforms.matchRegion.kept-exonLen25.bed" > "mapping_quant/${iso}.bed" || true
            tr ',' '\n' <<< "${iso_reads}" > "mapping_quant/${iso}_${s}.readIDs.txt"
            seqkit grep -n -f "mapping_quant/${iso}_${s}.readIDs.txt" "${reads}" > "mapping_quant/${iso}_${s}.fastq"
            minimap2 -t "${THREADS}" -ax splice -uf -R "@RG\tID:${iso}_${s}\tSM:${iso}" "${REF}" \
                "mapping_quant/${iso}_${s}.fastq" | samtools sort -@ "${THREADS}" --write-index -o "mapping_quant/${iso}_${s}.bam"
        done < "${q}.${s}.${TYPE}.isoform.read.map.txt"
    done
    for bed in mapping_quant/*.bed; do
        iso=$(basename "${bed}" .bed)
        samtools merge -f -@ "${THREADS}" --write-index "mapping_quant/${iso}.bam" mapping_quant/"${iso}"_*.bam
        cat mapping_quant/"${iso}"_*.fastq > "mapping_quant/${iso}.fastq"
    done
}

# -----------------------------------------------------------------------------
# qc: compare the splice junctions of every supporting read with the model
# -----------------------------------------------------------------------------
qc() {
    local iso="$1" d="polishing_quant/$1"
    local prev; prev=$(grep -P "\t${iso}(\t|$)" "${replace}" | cut -f1)   # flair name of the model
    mkdir -p "${d}"
    ln -sfr "mapping_quant/${prev}.bam"     "${d}/${iso}.bam"
    ln -sfr "mapping_quant/${prev}.bam.csi" "${d}/${iso}.bam.csi"
    ln -sfr "mapping_quant/${prev}.fastq"   "${d}/${iso}.fastq"
    seqkit grep -n -p "${iso}" "${models_fa}" > "${d}/${iso}.model.fa"
    minimap2 -t "${THREADS}" -ax splice -uf -R "@RG\tID:${iso}\tSM:${iso}" "${REF}" "${d}/${iso}.model.fa" \
        | samtools sort --write-index -o "${d}/${iso}.model.bam"
    # columns: read_id isoform_id softclip_left softclip_right n_insertions
    #          n_deletions delta_left delta_right n_exons
    python3 "${BIN}/match_splice_sites.py" "${d}/${iso}.model.bam" "${d}/${iso}.bam" "${d}/${iso}.output.tsv"
}

# -----------------------------------------------------------------------------
# candidates: list possible template reads, best first
#   keep: n_exons == expected, no 5' soft clip ($4 == 0), poly(A) present ($3 > 0)
#   sort: fewest indels, then most exons, then longest poly(A) soft clip
# -----------------------------------------------------------------------------
candidates() {
    setdir "$1"; local iso="$2" nex="$3"
    [[ -s "polishing_quant/${iso}/${iso}.output.tsv" ]] || qc "${iso}"
    awk -v e="${nex}" '$9==e && $4==0 && $3>0 {print $5+$6"\t"$0}' "polishing_quant/${iso}/${iso}.output.tsv" \
        | sort -k1,1n -k9,9rn -k4,4rn | head -n 30
}

# -----------------------------------------------------------------------------
# polish (+ phase)
# -----------------------------------------------------------------------------
polish() {
    setdir "$1"; local iso="$2" nex="$3" tmpl="$4" minq="$5" phase_s="${6:-}"
    local d="polishing_quant/${iso}"
    [[ -s "${d}/${iso}.output.tsv" ]] || qc "${iso}"
    echo "== ${iso}  (exons=${nex}, template=${tmpl}, QUAL>${minq}, phasing sample=${phase_s:-none})"

    # reads with the model's splice structure, template excluded
    echo "${tmpl}" > "${d}/selected_read_template.txt"
    grep -v -e read_id -e "${tmpl}" "${d}/${iso}.output.tsv" | cut -f1 > "${d}/selected_reads_for_correction.txt"
    seqkit grep -n -f "${d}/selected_reads_for_correction.txt" "${reads}" > "${d}/selected_reads_for_correction.fastq"
    seqkit grep -n -p "${tmpl}" "${reads}" | seqkit fq2fa > "${d}/template.fasta"

    # one round of racon polishing of the template read
    minimap2 -x map-ont -a "${d}/template.fasta" "${d}/selected_reads_for_correction.fastq" \
        | samtools sort -o "${d}/racon_1.sam"
    racon -t "${THREADS}" "${d}/selected_reads_for_correction.fastq" "${d}/racon_1.sam" "${d}/template.fasta" \
        > "${d}/${iso}.c1.fa"
    { echo ">${iso}"; tail -n 1 "${d}/${iso}.c1.fa"; } > "${d}/${iso}.m.fa"
    minimap2 -x splice -a -t "${THREADS}" -R "@RG\tID:${iso}\tSM:${iso}" "${REF}" "${d}/${iso}.m.fa" \
        | samtools sort --write-index -o "${d}/${iso}.m-spliced.bam"

    # align the supporting reads to the polished model, one read group per sample
    for s in ${SAMPLES}; do
        seqkit grep -n -r -p "${TYPE}-${s}-" "${d}/selected_reads_for_correction.fastq" > "${d}/reads_${s}.fastq"
        minimap2 -x map-ont -a -t "${THREADS}" -R "@RG\tID:${iso}_${TYPE}-${s}\tSM:${iso}_${TYPE}-${s}" \
            "${d}/${iso}.m.fa" "${d}/reads_${s}.fastq" | samtools sort --write-index -o "${d}/${iso}.reads.${s}.bam"
    done
    samtools merge -f --write-index "${d}/${iso}.reads.bam" "${d}/${iso}".reads.*-*.bam
    rm "${d}/${iso}".reads.*-*.bam* "${d}"/reads_*.fastq

    # variant calling; QUAL threshold 300 for well supported models, lowered
    # (100 / 1) for models with few reads -- see parameter table
    samtools faidx "${d}/${iso}.m.fa"
    freebayes -f "${d}/${iso}.m.fa" "${d}/${iso}.reads.bam" | vcffilter -f "QUAL > ${minq}" \
        | bgzip > "${d}/${iso}.m.freebayes.vcf.gz"
    tabix -f -p vcf "${d}/${iso}.m.freebayes.vcf.gz"

    if [[ -z "${phase_s}" ]]; then
        # single haplotype: hap1 = hap2 = polished sequence
        echo "   not phased (no heterozygous variants)"
        for h in 1 2; do
            ln -sfr "${d}/${iso}.m.fa" "${d}/${iso}.m-hap${h}.fa"
            ln -sfr "${d}/${iso}.m-spliced.bam" "${d}/${iso}.m-spliced-hap${h}.bam"
            ln -sfr "${d}/${iso}.m-spliced.bam.csi" "${d}/${iso}.m-spliced-hap${h}.bam.csi"
        done
        return
    fi
    # read-based phasing per sample, haplotagging of the reads
    whatshap phase --reference "${d}/${iso}.m.fa" -o "${d}/${iso}.m.phased.vcf.gz" \
        "${d}/${iso}.m.freebayes.vcf.gz" "${d}/${iso}.reads.bam"
    bcftools index -f "${d}/${iso}.m.phased.vcf.gz"
    whatshap haplotag -r "${d}/${iso}.m.fa" --output-haplotag-list "${d}/${iso}.reads.haplotagged.txt" \
        -o "${d}/${iso}.reads.haplotagged.bam" "${d}/${iso}.m.phased.vcf.gz" "${d}/${iso}.reads.bam"
    samtools index "${d}/${iso}.reads.haplotagged.bam"
    # haplotype sequences from the phased genotype of ONE sample (column
    # phasing_sample, mostly C4_10 = the sample with most KLK4 reads)
    for h in 1 2; do
        bcftools consensus -f "${d}/${iso}.m.fa" -s "${iso}_${phase_s}" -H "${h}" \
            "${d}/${iso}.m.phased.vcf.gz" > "${d}/${iso}.m-hap${h}.fa"
        minimap2 -x splice -a -t "${THREADS}" -R "@RG\tID:${iso}\tSM:${iso}" "${REF}" "${d}/${iso}.m-hap${h}.fa" \
            | samtools sort --write-index -o "${d}/${iso}.m-spliced-hap${h}.bam"
    done
}

case "${mode}" in
    remap)      remap "$@" ;;
    candidates) candidates "$@" ;;
    polish)     polish "$@" ;;
    all)
        # columns: read_set isoform_id expected_exons template_read freebayes_min_qual
        #          phasing_sample phased in_final_set final_name note
        tail -n +2 "${PARAMS}" | while IFS=$'\t' read -r set iso nex tmpl minq phase_s _ _ _ _; do
            [[ -n "${tmpl}" && -n "${minq}" ]] || { echo "skip ${iso} (discarded during QC)"; continue; }
            ( polish "${set}" "${iso}" "${nex}" "${tmpl}" "${minq}" "${phase_s}" )
        done ;;
    *) echo "unknown mode ${mode}"; exit 1 ;;
esac

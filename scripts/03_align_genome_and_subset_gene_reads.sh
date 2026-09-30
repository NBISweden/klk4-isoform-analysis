#!/usr/bin/env bash
# =============================================================================
# Step 03 -- Genome-wide alignment and extraction of the reads of one gene
# -----------------------------------------------------------------------------
# Manuscript: section 2.7
#
# Isoform reconstruction is done per gene. All reads of a sample are aligned to
# the complete GRCh38.p14 genome, so that reads of paralogous KLK genes are
# placed at their best locus. Only reads with a primary alignment overlapping
# the target gene region (config/gene_regions.tsv) are kept.
# Read names are set to "<TYPE>-<sample>-<n>" (e.g. dRNA-22_10-103); the sample
# name in the read name is used in all later steps to split reads by sample.
#
# minimap2 options:
#   -ax splice             spliced long-read alignment
#   -ub                    try both strands for splice motifs
#   -k14 -w4               smaller k-mer / window -> more sensitive for ONT reads
#   -y                     copy FASTQ comments to the SAM output
#   --junc-bed junctions   bonus for annotated splice junctions (GENCODE)
#   --secondary=no         no secondary alignments
# samtools view -F 3844    drop unmapped, secondary, QC-fail, duplicate and
#                          supplementary alignments
#
# Input : dRNA: ${READDIR}/dRNA/<dRNA_fastq>;  cDNA: output of step 02
# Output: ${WORKDIR}/${GENE}/flair_ALL/${TYPE}_<sample>-${GENE}.fastq  (one per sample)
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"
ALN="${WORKDIR}/genome_alignment/${TYPE}"
OUT="${WORKDIR}/${GENE}/flair_ALL"
mkdir -p "${ALN}" "${OUT}"

# annotated splice junctions in BED format (paftools.js is part of minimap2)
JUNC="${REFDIR}/junctions.bed"
[[ -s "${JUNC}" ]] || paftools.js gff2bed "${REF_ANNOT}" > "${JUNC}"

# target region of the gene
region=$(awk -v g="${GENE}" '$1==g {print $2":"$3"-"$4}' "${GENE_REGIONS}")
[[ -n "${region}" ]] || { echo "unknown gene ${GENE}"; exit 1; }

tail -n +2 "${CONFIG}/samples.tsv" | while IFS=$'\t' read -r sample _ _ _ _ _ drna_fq _ _ _; do
    if [[ "${TYPE}" == "dRNA" ]]; then
        fq="${READDIR}/dRNA/${drna_fq}"
    else
        fq="${WORKDIR}/pychopper/cDNA_${sample}_phmm_U.pychopper.fq.gz"
    fi
    bam="${ALN}/${TYPE}_${sample}.GRCh38.bam"

    # genome-wide spliced alignment
    minimap2 -t "${THREADS}" -ax splice -ub -k14 -w 4 -y --junc-bed "${JUNC}" --secondary=no \
        "${REF_GENOME}" "${fq}" | samtools sort -@ "${THREADS}" --write-index -o "${bam}"

    # primary alignments overlapping the gene region -> FASTQ with sample-tagged names
    samtools view -h -F 3844 "${bam}" "${region}" \
        | samtools fastq - \
        | seqkit replace -p '.+' -r "${TYPE}-${sample}-{nr}" \
        > "${OUT}/${TYPE}_${sample}-${GENE}.fastq"
done

seqkit stats "${OUT}"/${TYPE}_*-"${GENE}".fastq

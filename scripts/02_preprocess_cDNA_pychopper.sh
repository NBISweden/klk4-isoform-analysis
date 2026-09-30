#!/usr/bin/env bash
# =============================================================================
# Step 02 -- Identify, orient and trim full-length cDNA reads (cDNA only)
# -----------------------------------------------------------------------------
# Manuscript: sections 2.6.2 / 2.7
#
# Why: PCR-cDNA reads (SQK-PCB114.24) come in both orientations and carry
#      primer/adapter sequence. Pychopper detects the primers, keeps
#      full-length reads, trims the primers and orients all reads 5'->3'.
#      Direct RNA reads are always sequenced 3'->5' as native RNA and do NOT
#      need this step.
#
#   -k PCB114   primer configuration of the cDNA-PCR Barcoding Kit V14
#   -m phmm     profile-HMM based primer detection (more sensitive than edlib)
#   -U          report unclassified reads separately
#   -y          autotune the cutoff parameter
#
# Result: >94 % of the cDNA reads were classified as full length for all samples.
#
# Input : ${READDIR}/cDNA/<cDNA_fastq from config/samples.tsv>  (scripts/basecalling_dorado.sh)
# Output: ${WORKDIR}/pychopper/<sample>_phmm_U.pychopper.fq.gz
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"
OUT="${WORKDIR}/pychopper"
mkdir -p "${OUT}"

# columns of config/samples.tsv: 1 sample_id ... 8 cDNA_fastq
tail -n +2 "${CONFIG}/samples.tsv" | while IFS=$'\t' read -r sample _ _ _ _ _ _ cdna_fq _; do
    x="${READDIR}/cDNA/${cdna_fq}"
    n="${OUT}/cDNA_${sample}"
    pychopper -U -t "${THREADS}" -y -k PCB114 -m phmm \
        -S "${n}_phmm_U.stats.tsv" -A "${n}_phmm_U.aln.bed" -r "${n}_phmm_U.report.pdf" \
        -u "${n}_phmm_U.unclassified.fq" -w "${n}_phmm_U.rescued.fq" \
        "${x}" "${n}_phmm_U.pychopper.fq" > "${n}_phmm_U.log" 2>&1
    seqkit stats -a "${n}_phmm_U.pychopper.fq" > "${n}_phmm_U.pychopper.stats.txt"
    gzip "${n}_phmm_U.pychopper.fq"
done

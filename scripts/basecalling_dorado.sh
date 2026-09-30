#!/usr/bin/env bash
# =============================================================================
# Basecalling of the raw ONT signal data (POD5) with dorado
# -----------------------------------------------------------------------------
# Manuscript: sections 2.6.1 / 2.6.2
#
# The samples were sequenced in two batches (config/samples.tsv):
#
#  batch | type | samples                                      | dorado | model (command line)
#  ------|------|----------------------------------------------|--------|------------------------------
#  1     | cDNA | C4_10 (barcode01), C4_10W (barcode02)        | 0.7.1  | sup
#  1     | dRNA | C4_10W                                       | 0.7.1  | sup,m6A,pseU
#  2     | cDNA | C4_0 (03), C4_0W (06)  flow cell 1           | 0.9.0  | sup
#        |      | 22_0 (04), 22_0W (07)  flow cell 2           |        |
#        |      | 22_10 (05), 22_10W (08) flow cell 3          |        |
#  2     | dRNA | C4_0, C4_10, 22_0, 22_10, C4_0W, 22_0W, 22_10W (one flow cell each) | 0.9.1 | sup@v5.1.0,m5C,m6A_DRACH,pseU
#
# Models:
#   cDNA (FLO-PRO114M, SQK-PCB114-24):  dna_r10.4.1_e8.2_400bps_sup@v5.0.0
#   dRNA (FLO-PRO004RA, SQK-RNA004):    rna004_130bps_sup@v5.0.0 (batch 1),
#                                       rna004_130bps_sup@v5.1.0 (batch 2)
#   Base modifications were called but not used in this study.
#
# Options
#   -x cuda:all / cuda:0       GPU basecalling
#   --estimate-poly-a          estimate poly(A) tail length (pt:i tag)
#   --no-trim                  (cDNA) keep primers: Pychopper (step 02) uses them
#                              to identify and orient full-length reads
#   --kit-name SQK-PCB114-24   (cDNA) barcode classification
#   --batchsize 80             (dRNA batch 2) GPU memory setting
#
# Post-processing
#   cDNA : dorado demux --no-trim --kit-name SQK-PCB114-24 (per flow cell)
#   all  : keep reads with mean Q >= 10 (samtools view -e '[qs] >= 10') and
#          convert to FASTQ (samtools bam2fq)
#
# Usage: POD5_DIR=<dir with one sub-folder per flow cell> basecalling_dorado.sh <run>
#   <run> = cdna_batch2 | drna_batch2 | cdna_batch1 | drna_batch1
# =============================================================================
set -euo pipefail

POD5_DIR="${POD5_DIR:?set POD5_DIR}"
DORADO_071="${DORADO_071:-dorado-0.7.1-linux-x64/bin/dorado}"
DORADO_090="${DORADO_090:-dorado-0.9.0-linux-x64/bin/dorado}"
DORADO_091="${DORADO_091:-dorado-0.9.1-linux-x64/bin/dorado}"

# basecall all POD5 files of a folder into one BAM per POD5 file
basecall() {   # <dorado> <pod5 dir> <out dir> <dorado options...>
    local dorado="$1" in="$2" out="$3"; shift 3
    mkdir -p "${out}"
    for x in "${in}"/*.pod5; do
        "${dorado}" basecaller "$@" "${x}" > "${out}/$(basename "${x%.*}").bam"
    done
}

# Q >= 10 filter and FASTQ conversion
q10_fastq() {  # <bam> <fastq.gz>
    samtools view -e '[qs] >= 10' -@8 "$1" | samtools bam2fq -@8 -t | bgzip -@8 -c > "$2"
    seqkit stats -a -j 8 "$2" > "${2%.fastq.gz}.stats.txt"
}

# -----------------------------------------------------------------------------
# cDNA, batch 2 (dorado 0.9.0): three flow cells with two barcoded samples each
# -----------------------------------------------------------------------------
cdna_batch2() {
    local -A samples=( [fc1_pool2]="03:C4_0 06:C4_0W" [fc1_pool3]="04:22_0 07:22_0W" [fc1_pool4]="05:22_10 08:22_10W" )
    for fc in "${!samples[@]}"; do
        basecall "${DORADO_090}" "${POD5_DIR}/${fc}" "${fc}" \
            -x cuda:all -v --estimate-poly-a --no-trim --kit-name SQK-PCB114-24 sup
        "${DORADO_090}" demux -o "${fc}/demultiplex_noTrim" -v --no-trim --kit-name SQK-PCB114-24 -t 16 "${fc}"
        for s in ${samples[$fc]}; do
            bc=${s%%:*}; name=${s#*:}
            q10_fastq "$(ls "${fc}"/demultiplex_noTrim/*_SQK-PCB114-24_barcode"${bc}".bam)" \
                "${fc}_bc${bc#0}_${name}.pass-noTrim.fastq.gz"    # e.g. fc1_pool2_bc6_C4_0W.pass-noTrim.fastq.gz
        done
    done
}

# -----------------------------------------------------------------------------
# dRNA, batch 2 (dorado 0.9.1): one sample per flow cell
# -----------------------------------------------------------------------------
drna_batch2() {
    # POD5 folder name -> output prefix (as in config/samples.tsv)
    for d in C4_0_op_101_002 C4_10_op_101_003 22_0_op_101_004 22_10_op_101_005 \
             C4_0W_op_101_006 22_0W_op_101_007 22_10W_op_101_008; do
        basecall "${DORADO_091}" "${POD5_DIR}/${d}" "${d}" \
            -x cuda:0 --batchsize 80 -v --estimate-poly-a sup@v5.1.0,m5C,m6A_DRACH,pseU
        samtools merge -@8 -o "${d}.bam" "${d}"/*.bam
        q10_fastq "${d}.bam" "${d}_pass.fastq.gz"
    done
}

# -----------------------------------------------------------------------------
# cDNA, batch 1 (dorado 0.7.1): one flow cell, barcodes 01 (C4_10) and 02 (C4_10W)
# -----------------------------------------------------------------------------
cdna_batch1() {
    basecall "${DORADO_071}" "${POD5_DIR}/cDNA_batch1" cDNA_batch1 \
        -x cuda:all -v --estimate-poly-a --no-trim --kit-name SQK-PCB114-24 sup
    "${DORADO_071}" demux -o cDNA_batch1/demultiplex_noTrim -v --no-trim --kit-name SQK-PCB114-24 -t 8 cDNA_batch1
    q10_fastq cDNA_batch1/demultiplex_noTrim/*_barcode01.bam cDNA_bc1_C4_10.fastq.gz
    q10_fastq cDNA_batch1/demultiplex_noTrim/*_barcode02.bam cDNA_bc2_C4_10W.fastq.gz
}

# -----------------------------------------------------------------------------
# dRNA, batch 1 (dorado 0.7.1): sample C4_10W
# -----------------------------------------------------------------------------
drna_batch1() {
    basecall "${DORADO_071}" "${POD5_DIR}/dRNA_batch1" dRNA_batch1 \
        -x cuda:all -v --estimate-poly-a sup,m6A,pseU
    samtools merge -@8 -o dRNA_C4_10W.bam dRNA_batch1/*.bam
    q10_fastq dRNA_C4_10W.bam dRNA_C4_10W.fastq.gz
}

case "${1:-}" in
    cdna_batch2|drna_batch2|cdna_batch1|drna_batch1) "$1" ;;
    *) echo "usage: POD5_DIR=... $0 {cdna_batch2|drna_batch2|cdna_batch1|drna_batch1}"; exit 1 ;;
esac

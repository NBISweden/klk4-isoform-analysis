#!/usr/bin/env bash
# =============================================================================
# config.sh -- shared settings for all analysis steps
# -----------------------------------------------------------------------------
# Every step script sources this file:   source "$(dirname "$0")/config.sh"
# Adapt the paths in the "USER SETTINGS" block to your system before running.
# =============================================================================

# ---------------------------- USER SETTINGS ----------------------------------

# Root directory where all intermediate and final results are written.
export WORKDIR="${WORKDIR:-$PWD/analysis}"

# Directory holding the basecalled ONT reads, see config/samples.tsv.
export READDIR="${READDIR:-$PWD/data/reads}"

# Number of CPU threads used by multithreaded tools.
export THREADS="${THREADS:-20}"

# Read type to process: "dRNA" (direct RNA, SQK-RNA004) or "cDNA" (SQK-PCB114.24).
export TYPE="${TYPE:-dRNA}"

# Target gene. The workflow was developed for KLK4 and re-used for KLK2/KLK3.
export GENE="${GENE:-KLK4}"

# ---------------------------- DERIVED PATHS ----------------------------------

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export REPO_DIR
export BIN="${REPO_DIR}/scripts/bin"
export CONFIG="${REPO_DIR}/config"

# Reference files (created by step 00 and 01)
export REFDIR="${WORKDIR}/ref"
export REF_GENOME="${REFDIR}/GRCh38.p14.primary_assembly.genome.fa"  # GENCODE release 46 (GRCh38.p14)
export REF_ANNOT="${REFDIR}/gencode.v46.annotation.gtf"               # GENCODE release 46 annotation
export CHR19="${REFDIR}/chr19-cleaned.fa"                             # unmodified chr19, header "chr19"
export REF="${REFDIR}/chr19-cleaned-mod2.fa"                          # chr19 with the KLK4 exon-6 modification
export GTF_CHR19="${REFDIR}/chr19-cleaned-mod2-allKLKs_transcripts.fix.gtf"

# Gene regions used to subset reads from the genome-wide GRCh38.p14 alignments
# (1-based, inclusive), see config/gene_regions.tsv
export GENE_REGIONS="${CONFIG}/gene_regions.tsv"

# The eight samples. Naming: <cell line>-<nM DHT>[H]
#   22 = 22Rv1, C4 = C4-2B; 0/10 = nM DHT; suffix W = whole-cell lysate,
#   no suffix = pooled polysome fractions (see config/samples.tsv).
export SAMPLES="22_0W 22_0 22_10W 22_10 C4_0W C4_0 C4_10W C4_10"

# Gene id that flair appends to isoform names (Ensembl id of KLK4).
export GENE_ID_KLK4="ENSG00000167749.12"

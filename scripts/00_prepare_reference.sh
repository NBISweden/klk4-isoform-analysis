#!/usr/bin/env bash
# =============================================================================
# Step 00 -- Prepare the (modified) chr19 reference sequence
# -----------------------------------------------------------------------------
# Manuscript: section 2.7 "Data analysis"
#
# Why: Alignments of the KLK4 long reads revealed a 75 bp insertion and an
#      8 bp deletion in KLK4 exon 6 relative to GRCh38. Both variants are
#      supported by the vast majority of reads in both cell lines. To obtain
#      accurate spliced alignments and read assignments, the chr19 sequence
#      was edited accordingly.
#      The insertion occurs mostly as 75 bp, sporadically as 104/114/196 bp
#      tandem-repeat-like variants (see docs/workflow.md).
#
# Input : GENCODE GRCh38.p14 primary assembly (FASTA)
# Output: ${REFDIR}/chr19-cleaned.fa       chr19 with a clean FASTA header
#         ${REFDIR}/chr19-cleaned-mod2.fa  chr19 with the exon-6 modification
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"
mkdir -p "${REFDIR}"
cd "${REFDIR}"

# --- 1. Extract chr19 from the genome and simplify the header ----------------
samtools faidx "${REF_GENOME}" chr19 | sed 's/^>chr19.*/>chr19/' > chr19-cleaned.fa
samtools faidx chr19-cleaned.fa

# --- 2. Apply the 75 bp insertion and 8 bp deletion in KLK4 exon 6 -----------
python3 "${BIN}/modify_reference.py" chr19-cleaned.fa chr19-cleaned-mod2.fa
samtools faidx chr19-cleaned-mod2.fa

# --- 3. Sanity check: modified chr19 must be 58,617,683 bp (= +75 -8) --------
cut -f1,2 chr19-cleaned.fa.fai chr19-cleaned-mod2.fa.fai

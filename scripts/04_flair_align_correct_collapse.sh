#!/usr/bin/env bash
# =============================================================================
# Step 04 -- FLAIR: align, correct splice sites, collapse reads into isoforms
# -----------------------------------------------------------------------------
# Manuscript: section 2.7
#
# Reads of all 8 samples are pooled, so that one common set of isoform models
# is built; reads are assigned back to the samples in step 06.
#
# 4.1 flair align   -- minimap2 alignment to the MODIFIED chr19 + BED12 per read
#       --nvrna                native-RNA minimap2 settings (direct RNA)
#       --quality 20           minimum MAPQ 20
#       --remove_singleexon    drop unspliced reads (mostly degraded / fragments)
#       --filtertype separate  write chimeric alignments to a separate file
#
# 4.2 flair correct -- snap splice sites to the annotation (${GTF_CHR19});
#       reads whose junctions match annotated junctions -> *_all_corrected.bed
#       reads with unannotated junctions                -> *_all_inconsistent.bed
#
# 4.3 flair collapse -- build isoform models, run SEPARATELY on both read sets:
#       a) corrected    ("reference splice sites")  -> models named im*-a-*
#       b) inconsistent ("de novo splice sites")    -> models named im*-b-*
#       --stringent          supporting reads must cover 80 % of the model and
#                            span >= 25 bp of the first and last exon
#       --check_splice       4 of 6 bp around each splice site must match, no
#                            insertion > 3 bp at splice sites
#       --annotation_reliant generate  (set a only) match reads to annotated
#                            transcripts first, then detect novel ones
#       --generate_map       write the read -> isoform assignment
#
# Software: FLAIR 2.2.0 (development version, run from source)
#
# Input : ${WORKDIR}/${GENE}/flair_ALL/${TYPE}_<sample>-${GENE}.fastq (step 03)
# Output: ${WORKDIR}/${GENE}/flair_ALL/flair_${TYPE}-${GENE}/
#           ${TYPE}-${GENE}_noSE_all_corrected.isoforms.{bed,fa}
#           ${TYPE}-${GENE}_noSE_all_inconsistent.isoforms.{bed,fa}
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"
cd "${WORKDIR}/${GENE}/flair_ALL"

# merge the reads of all 22-* and C4-* libraries (read names keep the sample)
cat "${TYPE}"_*-"${GENE}".fastq > "${TYPE}-${GENE}.fastq"

x="${TYPE}-${GENE}.fastq"
pre="${x%.fastq}"               # e.g. dRNA-KLK4
out="flair_${pre}"
mkdir -p "${out}"

# 4.1 align
flair align -r "${x}" -g "${REF}" -f "${GTF_CHR19}" --nvrna --quality 20 \
    --remove_singleexon --filtertype separate --threads "${THREADS}" \
    -o "${out}/${pre}_noSE.align" > "${out}/align.out" 2> "${out}/align.err"

# 4.2 correct
flair correct -q "${out}/${pre}_noSE.align.bed" -f "${GTF_CHR19}" --nvrna --threads "${THREADS}" \
    -o "${out}/${pre}_noSE" > "${out}/correct.out" 2> "${out}/correct.err"

# 4.3a collapse reads with annotated splice sites
if [[ -s "${out}/${pre}_noSE_all_corrected.bed" ]]; then
    flair collapse -t "${THREADS}" -g "${REF}" --gtf "${GTF_CHR19}" -r "${x}" \
        -q "${out}/${pre}_noSE_all_corrected.bed" \
        --generate_map --check_splice --stringent --annotation_reliant generate \
        -o "${out}/${pre}_noSE_all_corrected" \
        > "${out}/collapse_corrected.out" 2> "${out}/collapse_corrected.err"
fi

# 4.3b collapse reads with novel ("inconsistent") splice sites
if [[ -s "${out}/${pre}_noSE_all_inconsistent.bed" ]]; then
    flair collapse -t "${THREADS}" -g "${REF}" --gtf "${GTF_CHR19}" -r "${x}" \
        -q "${out}/${pre}_noSE_all_inconsistent.bed" \
        --generate_map --stringent --check_splice \
        -o "${out}/${pre}_noSE_all_inconsistent" \
        > "${out}/collapse_inconsistent.out" 2> "${out}/collapse_inconsistent.err"
fi

# Result for dRNA KLK4 (17,780 reads):
#   flair align  : 14,509 reads kept (268 chimeric, 409 MAPQ<20, 2,594 single-exon removed)
#   flair correct: 6,926 corrected / 7,583 inconsistent reads

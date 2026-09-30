#!/usr/bin/env bash
# =============================================================================
# Step 01 -- Transfer the KLK gene annotation to the modified chr19
# -----------------------------------------------------------------------------
# Manuscript: section 2.7 "Data analysis"
#
# Why: The 75 bp insertion / 8 bp deletion (step 00) shift all coordinates
#      downstream of chr19:50,906,558 by +67 bp and change KLK4 exon 6.
#      Instead of shifting coordinates by hand, every reference transcript is
#      extracted as a sequence and re-aligned (spliced) to the modified chr19.
#      The alignment is then converted back to GTF:
#
#        GTF --gffread--> transcript FASTA --minimap2--> BAM --bedtools--> BED12
#            --bedToGenePred--> genePred --genePredToGtf--> GTF
#
#      The round trip loses the gene_id (every gene_id becomes the transcript
#      id), so gene_id is restored afterwards ("*.fix.gtf"). flair needs a
#      correct gene_id to name isoforms "<model>_<gene_id>".
#
# Input : ${CHR19}, ${REF}, chr19 subset of the GENCODE annotation (ano_chr19.gtf)
#         ref/KLK4T2.cdna.fa (KLK4T2 fusion transcript, Lundwall et al. 2021)
# Output: ${GTF_CHR19}  = merged annotation of 15 KLK genes + KLK4T2 on the
#         modified chr19 (input for flair in steps 04-06)
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/config.sh"
cd "${REFDIR}"

# chr19 subset of the GENCODE release 46 (GRCh38.p14) annotation.
ANNOT_CHR19="${REFDIR}/ano_chr19.gtf"
[[ -s "${ANNOT_CHR19}" ]] || awk '$1=="chr19"' "${REF_ANNOT}" > "${ANNOT_CHR19}"

# KLK genes of the locus used for FLAIR (KLK1 is not included).
GENES="KLK2 KLK3 KLK4 KLK5 KLK6 KLK7 KLK8 KLK9 KLK10 KLK11 KLK12 KLK13 KLK14 KLK15 KLKP1"

for gene in ${GENES}; do
    echo "== ${gene}"
    # 1. all annotation lines of the gene
    grep "gene_name \"${gene}\";" "${ANNOT_CHR19}" > "${gene}.gtf"

    # 2. transcript sequences from the UNMODIFIED chr19
    gffread "${gene}.gtf" -g "${CHR19}" -w "${gene}_transcripts.fasta"

    # 3. spliced alignment of the transcripts to the MODIFIED chr19
    #    -x splice:hq : high-quality spliced alignment (transcripts are error free)
    #    -uf          : use the transcript strand (forward) for splice motifs
    minimap2 -ax splice:hq -uf "${REF}" "${gene}_transcripts.fasta" \
        | samtools sort --write-index -o "chr19-cleaned-mod2-${gene}_transcripts.bam"

    # 4. BAM -> BED12 -> genePred -> GTF (UCSC tools)
    bedtools bamtobed -bed12 -i "chr19-cleaned-mod2-${gene}_transcripts.bam" \
        > "chr19-cleaned-mod2-${gene}_transcripts.bed"
    bedToGenePred "chr19-cleaned-mod2-${gene}_transcripts.bed" \
        "chr19-cleaned-mod2-${gene}_transcripts.genepred"
    genePredToGtf file "chr19-cleaned-mod2-${gene}_transcripts.genepred" \
        "chr19-cleaned-mod2-${gene}_transcripts.gtf"

    # 5. restore the Ensembl gene_id (one gene_id per gene)
    gene_id=$(awk '$3=="gene"' "${gene}.gtf" | grep -o 'gene_id "[^"]*"' | cut -d'"' -f2 | sort -u)
    [[ $(wc -w <<< "${gene_id}") -eq 1 ]] || { echo "ERROR: ${gene} has gene_ids: ${gene_id}"; exit 1; }
    sed "s/gene_id \"[^\"]*\"/gene_id \"${gene_id}\"/" \
        "chr19-cleaned-mod2-${gene}_transcripts.gtf" > "chr19-cleaned-mod2-${gene}_transcripts.fix.gtf"
done
# e.g. KLK4 ENSG00000167749.12, KLK3 ENSG00000142515.16, KLK15 ENSG00000174562.14

# --- KLK4T2 fusion transcript (not in GENCODE) -------------------------------
# KLK4T2 (KLK4 exons fused to KLKP1 sequence; Lundwall et al. 2021, IJMS 22:13023)
# is aligned from its cDNA sequence; its gene_id stays "KLK4T2".
minimap2 -ax splice:hq -uf "${REF}" "${REPO_DIR}/ref/KLK4T2.cdna.fa" \
    | samtools sort --write-index -o chr19-cleaned-mod2-KLK4T2_transcripts.bam
bedtools bamtobed -bed12 -i chr19-cleaned-mod2-KLK4T2_transcripts.bam > chr19-cleaned-mod2-KLK4T2_transcripts.bed
bedToGenePred chr19-cleaned-mod2-KLK4T2_transcripts.bed chr19-cleaned-mod2-KLK4T2_transcripts.genepred
genePredToGtf file chr19-cleaned-mod2-KLK4T2_transcripts.genepred chr19-cleaned-mod2-KLK4T2_transcripts.gtf

# --- Merge into one annotation for flair -------------------------------------
cat chr19-cleaned-mod2-KLK*_transcripts.fix.gtf chr19-cleaned-mod2-KLK4T2_transcripts.gtf \
    > "${GTF_CHR19}"

# BED12 of the KLK4 reference transcripts, used as "--ref_isoforms" when plotting.
# A copy of the result is provided in ref/KLK4_reference_transcripts.bed

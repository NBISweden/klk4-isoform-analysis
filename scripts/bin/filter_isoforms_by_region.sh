#!/usr/bin/env bash

# --- Usage check ---
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <isoforms.bed12> <regions.bed> <matched_out.bed12> <removed_out.bed12>"
    exit 1
fi

isoforms="$1"
regions="$2"
matched="$3"
removed="$4"

# --- Step 1: Extract exons from isoforms ---
bedtools bed12tobed6 -i "$isoforms" > tmp_exons.bed

# --- Step 2: Find which isoforms have any exon overlapping the regions ---
bedtools intersect -u -a tmp_exons.bed -b "$regions" \
    | cut -f4 \
    | sort -u \
    > tmp_matching_ids.txt

# --- Step 3: Filter isoforms into matched and removed ---
grep -F -w -f tmp_matching_ids.txt "$isoforms" > "$matched"
grep -F -w -v -f tmp_matching_ids.txt "$isoforms" > "$removed"

# --- Cleanup ---
rm tmp_exons.bed tmp_matching_ids.txt


#!/usr/bin/env bash

# --- Usage check ---
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <isoforms.bed12> <min_exon_length> <kept_out.bed12> <removed_out.bed12>"
    echo "Example: $0 isoforms.bed12 50 kept.bed12 removed.bed12"
    exit 1
fi

isoforms="$1"
minlen="$2"
kept="$3"
removed="$4"

# --- Step 1: Convert BED12 to exon-level BED6 ---
bedtools bed12tobed6 -i "$isoforms" > tmp_exons.bed

# --- Step 2: Collect isoform IDs that have at least one short exon ---
awk -v min="$minlen" '{ if (($3 - $2) < min) print $4 }' tmp_exons.bed \
    | sort -u \
    > tmp_short_ids.txt

# --- Step 3: Split original BED12 ---
grep -F -w -v -f tmp_short_ids.txt "$isoforms" > "$kept"
grep -F -w    -f tmp_short_ids.txt "$isoforms" > "$removed"

# --- Cleanup ---
rm tmp_exons.bed tmp_short_ids.txt


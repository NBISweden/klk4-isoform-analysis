#!/usr/bin/env python3
"""
rank_and_rename_models.py -- rank FLAIR isoform models by read support and rename them.

Models are sorted by their total read count over all samples (descending) and
renamed to

    im<rank>-<set>-R<reads>-P<percent>_<gene_id>      e.g. im1-a-R2344-P45.1_ENSG00000167749.12

  rank     running number, starting at --start (set b continues after set a)
  set      a = corrected (annotated splice sites), b = inconsistent (novel splice sites)
  reads    total read count (flair quantify)
  percent  share of the reads of all retained models of the set, truncated to
           one decimal, without leading zero (0.78 -> ".7", 0.04 -> "0")
  gene_id  the id of the target gene (--gene_id), for all models

Models with fewer than --min_reads reads are dropped
(dRNA KLK4: set a --min_reads 0, set b --min_reads 10).

Outputs (next to the counts file, <prefix> = counts file without ".counts.tsv"):
  <prefix>.counts_sorted.renamed.tsv   count table, sorted and renamed
  <prefix>.counts.replaceName.tsv      old name <tab> new name
  <prefix>.renamed.fa / .renamed.bed   renamed models (if --fasta / --bed given)

Usage:
  rank_and_rename_models.py <prefix>.counts.tsv --set a --start 1 --min_reads 0 \
      --fasta models.fa --bed models.bed
"""
import argparse
import re


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("counts", help="flair quantify counts table (*.counts.tsv)")
    ap.add_argument("--set", required=True, choices=["a", "b"])
    ap.add_argument("--start", type=int, default=1, help="rank of the first model")
    ap.add_argument("--min_reads", type=int, default=0)
    ap.add_argument("--gene_id", default="ENSG00000167749.12", help="gene id appended to all names (default: KLK4)")
    ap.add_argument("--fasta", help="model FASTA to rename")
    ap.add_argument("--bed", help="model BED12 to rename")
    a = ap.parse_args()

    with open(a.counts) as fh:
        header = fh.readline()
        rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    rows = [r for r in rows if sum(float(x) for x in r[1:]) >= a.min_reads]
    rows.sort(key=lambda r: -sum(float(x) for x in r[1:]))
    total = sum(sum(float(x) for x in r[1:]) for r in rows)

    rename = {}
    for i, r in enumerate(rows):
        old = r[0]
        model = old.rsplit("_", 1)[0]
        reads = int(round(sum(float(x) for x in r[1:])))
        pct = f"{int(1000 * reads / total) / 10:.1f}"
        pct = "0" if pct == "0.0" else (pct[1:] if pct.startswith("0") else pct)
        # keep the reference transcript / KLK4T2 name if the model matches one
        ref = "" if re.match(r"^(dRNA|cDNA)-", model) else f"{re.sub(r'_$', '', model)}-"
        rename[old] = f"im{a.start + i}-{a.set}-{ref}R{reads}-P{pct}_{a.gene_id}"

    prefix = re.sub(r"\.counts\.tsv$", "", a.counts)
    with open(prefix + ".counts_sorted.renamed.tsv", "w") as out:
        out.write(header)
        for r in rows:
            out.write("\t".join([rename[r[0]]] + r[1:]) + "\n")
    with open(prefix + ".counts.replaceName.tsv", "w") as out:
        for old, new in rename.items():
            out.write(f"{old}\t{new}\n")

    if a.fasta:
        with open(a.fasta) as fh, open(prefix + ".renamed.fa", "w") as out:
            keep = False
            for line in fh:
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    keep = name in rename
                    if keep:
                        out.write(f">{rename[name]}\n")
                elif keep:
                    out.write(line)
    if a.bed:
        with open(a.bed) as fh, open(prefix + ".renamed.bed", "w") as out:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[3] in rename:
                    f[3] = rename[f[3]]
                    out.write("\t".join(f) + "\n")
    print(f"{len(rename)} models renamed (im{a.start}..im{a.start + len(rename) - 1})")


if __name__ == "__main__":
    main()

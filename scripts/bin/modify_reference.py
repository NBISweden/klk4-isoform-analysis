#!/usr/bin/env python3
"""
modify_reference.py -- create the modified chr19 reference used for all RNA alignments.

The KLK4 long reads of both cell lines (22Rv1, C4-2B) consistently support two
differences to GRCh38 inside KLK4 exon 6:

  1. a 75 bp insertion directly after chr19:50,906,558
  2. an 8 bp deletion of chr19:50,906,576-50,906,583 (sequence TTGGGGGC)

Both coordinates are 1-based and refer to the *unmodified* GRCh38 chr19.
The modified chromosome is 67 bp longer (58,617,616 -> 58,617,683 bp).

Usage:
    python modify_reference.py chr19-cleaned.fa chr19-cleaned-mod2.fa
"""
import sys

INSERT_AFTER = 50906558            # 1-based position; insertion follows this base
INSERT_SEQ = ("GTCTGAGGGAGGAGGGGTTGGGGGCCTGGACCCCTGGGTCTGAGGG"
              "AGGAGGGGCTGGGGGATCTGGACTCCTGG")          # 75 bp
DEL_START, DEL_END = 50906576, 50906583                # 1-based, inclusive
DEL_SEQ = "TTGGGGGC"                                   # expected deleted bases
LINE_WIDTH = 60


def read_single_fasta(path):
    """Return (header, sequence) of a single-record FASTA file."""
    header, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    sys.exit(f"ERROR: {path} contains more than one sequence")
                header = line[1:].strip()
            else:
                chunks.append(line.strip())
    return header, "".join(chunks)


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    header, seq = read_single_fasta(sys.argv[1])

    assert len(INSERT_SEQ) == 75
    observed = seq[DEL_START - 1:DEL_END].upper()
    if observed != DEL_SEQ:
        sys.exit(f"ERROR: expected {DEL_SEQ} at {DEL_START}-{DEL_END}, found {observed}")

    # Python slices are 0-based / end-exclusive:
    #   seq[:50906558]           -> bases 1..50,906,558
    #   seq[50906558:50906575]   -> bases 50,906,559..50,906,575 (kept)
    #   seq[50906583:]           -> bases 50,906,584..end (after the deletion)
    mod = (seq[:INSERT_AFTER] + INSERT_SEQ
           + seq[INSERT_AFTER:DEL_START - 1] + seq[DEL_END:])

    with open(sys.argv[2], "w") as out:
        out.write(f">{header.split()[0]}\n")
        for i in range(0, len(mod), LINE_WIDTH):
            out.write(mod[i:i + LINE_WIDTH] + "\n")
    print(f"{header.split()[0]}: {len(seq)} bp -> {len(mod)} bp", file=sys.stderr)


if __name__ == "__main__":
    main()

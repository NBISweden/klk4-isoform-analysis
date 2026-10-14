#!/usr/bin/env python3
import pysam
import sys

if len(sys.argv) != 4:
    print("Usage: python match_splice_sites.py <isoforms.bam> <reads.bam> <output.tsv>")
    sys.exit(1)

iso_bam = sys.argv[1]
read_bam = sys.argv[2]
out_tsv = sys.argv[3]

def get_junctions(aln):
    juncs = []
    ref_pos = aln.reference_start
    for op, length in aln.cigartuples:
        if op in (0,7,8):  # M, =, X
            ref_pos += length
        elif op == 3:  # N
            donor = ref_pos
            ref_pos += length
            acceptor = ref_pos
            juncs.append((donor, acceptor))
        elif op == 2:  # D
            ref_pos += length
        elif op == 1:  # I does not consume ref
            pass
    return tuple(juncs)

def get_softclips(cigartuples):
    left = cigartuples[0][1] if cigartuples and cigartuples[0][0] == 4 else 0
    right = cigartuples[-1][1] if cigartuples and cigartuples[-1][0] == 4 else 0
    return left, right

def get_indels(cigartuples):
    ins = sum(l for op,l in cigartuples if op == 1)
    dels = sum(l for op,l in cigartuples if op == 2)
    return ins, dels

# --- Read isoform junctions ---
isoform_junctions = {}
isoform_bounds = {}
iso = pysam.AlignmentFile(iso_bam, "rb")
for aln in iso.fetch():
    if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
        continue
    jid = aln.query_name
    juncs = get_junctions(aln)
    strand = "-" if aln.is_reverse else "+"
    isoform_junctions[jid] = (aln.reference_name, strand, juncs)
    isoform_bounds[jid] = (aln.reference_start, aln.reference_end)
iso.close()

# --- Match reads ---
out = open(out_tsv, "w")
out.write("read_id\tisoform_id\tsoftclip_left\tsoftclip_right\tn_insertions\tn_deletions\tdelta_left\tdelta_right\tn_exons\n")

reads = pysam.AlignmentFile(read_bam, "rb")
for aln in reads.fetch():
    if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
        continue

    rchrom = aln.reference_name
    rstrand = "-" if aln.is_reverse else "+"
    rjuncs = get_junctions(aln)

    softL, softR = get_softclips(aln.cigartuples)
    ins, dels = get_indels(aln.cigartuples)
    nexons = len(rjuncs) + 1

    matched = None
    dleft = dright = "NA"

    for iid, (ichrom, istr, ijuncs) in isoform_junctions.items():
        if ichrom == rchrom and istr == rstrand and ijuncs == rjuncs:
            matched = iid
            istart, iend = isoform_bounds[iid]
            dleft = aln.reference_start - istart
            dright = aln.reference_end - iend
            break

    if matched:
        out.write(f"{aln.query_name}\t{matched}\t{softL}\t{softR}\t{ins}\t{dels}\t{dleft}\t{dright}\t{nexons}\n")

reads.close()
out.close()


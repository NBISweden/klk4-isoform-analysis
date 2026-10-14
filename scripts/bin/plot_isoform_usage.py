#!/usr/bin/env python3
import matplotlib

matplotlib.use("Agg")
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import seaborn as sns

parser = argparse.ArgumentParser(
    description="Plot isoform structures and usage, with optional reference isoforms."
)
parser.add_argument("isoforms", type=str, help="isoforms in BED format")
parser.add_argument("counts_matrix", type=str, help="counts matrix (TSV)")
parser.add_argument(
    "gene_name", type=str, help="Name of gene, must correspond with gene names in files"
)
parser.add_argument(
    "-o", dest="o", type=str, help="prefix used for output files (default=gene_name)"
)
parser.add_argument(
    "--min_reads",
    dest="min_reads",
    type=int,
    default=6,
    help="minimum supporting reads for an isoform (default=6)",
)
parser.add_argument(
    "--img_resolution",
    dest="img_resolution",
    type=int,
    default=600,
    help="image resolutions for isoform and usage plots (default=600)",
)
parser.add_argument(
    "--min_isoform_spacing",
    dest="min_isoform_spacing",
    type=int,
    default=5,
    help="Minimum distance (bp) required between isoforms to place them on the same row",
)
parser.add_argument(
    "--palette", dest="palette", help="palette file (one hex color per line)"
)
parser.add_argument(
    "--sort_by_reads",
    action="store_true",
    help="sort isoforms in plot by number of supporting reads (highest on TOP)",
)
parser.add_argument(
    "--ref_isoforms",
    dest="ref_isoforms",
    help="BED file with reference isoform annotations",
)
parser.add_argument(
    "--ref_on_top",
    action="store_true",
    help="draw reference isoforms above sample isoforms",
)
parser.add_argument(
    "--vcf",
    dest="vcf",
    help="Optional VCF file with variants to overlay (not required)",
)
parser.add_argument(
    "--xmargin_left_frac",
    type=float,
    default=0.05,
    help="Fraction of genomic width to add as margin left for labels (default=0.05)",
)
parser.add_argument(
    "--xmargin_right_frac",
    type=float,
    default=0.05,
    help="Fraction of genomic width to add as margin right for labels (default=0.05)",
)
parser.add_argument(
    "--label_padding_scale",
    type=float,
    default=1.0,
    help="Scale factor for estimated label width used during row packing (default=1.0)",
)
args = parser.parse_args()

# Open files
bedfh = open(args.isoforms)
counts_matrix = open(args.counts_matrix)
if not args.o:
    args.o = args.gene_name


def iso_label(s):
    idx = s.rfind("_")
    return s[:idx] if idx != -1 else s


def estimate_label_padding_bp(
    name,
    locus_width_bp,
    fig_width_in=9.0,
    panel_width_frac=0.95,
    fontsize_pt=8,
    char_width_factor=0.62,
):
    """
    Estimate label width in base pairs so row packing can reserve space for text.
    """
    if locus_width_bp <= 0:
        return 0
    label = iso_label(name)
    panel_width_in = max(fig_width_in * panel_width_frac, 1e-6)
    bp_per_in = float(locus_width_bp) / panel_width_in
    label_in = max(1, len(label)) * (fontsize_pt / 72.0) * char_width_factor
    return int(max(1.0, label_in * bp_per_in * args.label_padding_scale))


def parse_bed(bedfh, keepiso=set()):
    """
    Parse sample BED file.
    Returns:
      info: list of [blocksizes, blockstarts_abs, color, flag, name]
      lowbound, upbound, strand, usednames
    Note: blockstarts are absolute genomic positions (not shifted).
    """
    info = []
    usednames = []
    lowbound, upbound = 1e18, 0
    strand = None
    for line in bedfh:
        line = line.rstrip("\n").split("\t")
        if len(line) < 12:
            continue
        name = line[3]
        start = int(line[1])
        end = int(line[2])
        blocksizes = [int(n) for n in line[10].split(",")[:-1]]
        blockstarts = [
            int(n) + start for n in line[11].split(",")[:-1]
        ]  # absolute positions

        # productivity flags handling (keep original semantics)
        if "_" in name[-4:]:
            flag = name[name.rfind("_") + 1 :]
            if flag not in ["PRO", "PTC"]:
                flag = "Z" + flag
            name = name[: name.rfind("_")]
        elif "_PTC" in name:
            flag = "PTC"
            name = name.replace("_PTC", "")
        elif "_PRO" in name:
            flag = "PRO"
            name = name.replace("_PRO", "")
        elif "_NST" in name:
            flag = "Z"
            name = name.replace("_NST", "")
        elif "_NGO" in name:
            flag = "Z"
            name = name.replace("_NGO", "")
        else:
            flag = "PRO"

        # keep only isoforms with colors assigned in keepiso (mapping from counts)
        # keepiso keys are isoform names from counts file (line[0])
        try:
            keep_check = (name in keepiso) or (name[: name.rfind("_")] in keepiso)
        except Exception:
            keep_check = name in keepiso
        if not keep_check:
            continue

        strand = line[5]
        lowbound = min(lowbound, start)
        upbound = max(upbound, end)

        color = keepiso.get(name, keepiso.get(name[: name.rfind("_")], None))
        info.append([blocksizes, blockstarts, color, flag, name])
        usednames.append(name)

    if not info:
        return [], lowbound, upbound, strand, []

    return info, lowbound, upbound, strand, usednames


def parse_bed_ref(bedfh):
    """
    Parse reference BED: return list of [blocksizes, blockstarts_abs, color, 'REF', name]
    (all blockstarts are absolute genomic coordinates).
    """
    entries = []
    for line in bedfh:
        line = line.rstrip("\n").split("\t")
        if len(line) < 12:
            continue
        name = line[3]
        start = int(line[1])
        blocksizes = [int(n) for n in line[10].split(",")[:-1]]
        blockstarts = [int(n) + start for n in line[11].split(",")[:-1]]
        entries.append([blocksizes, blockstarts, "black", "REF", name])
    return entries


def pack(data, tosort=True, min_spacing=5, strand="+", label_padding_bp=None):
    """
    Pack isoforms into rows to avoid overlap.

    - data: list of isoform entries as before: [blocksizes, blockstarts, color, flag, name]
            blockstarts are absolute genome coordinates
    - tosort: if True sort isoforms by their leftmost start before packing
    - min_spacing: minimum distance (bp) required between ANY exon-span of candidate isoform
                   and ANY exon-span already present on a row
    - strand: used to reserve row-space on label side only
    - label_padding_bp: dict of isoform-name -> estimated label width in bp
    Returns: list of rows (each row is a list of isoform entries)
    """
    if not data:
        return []

    strand = strand if strand in ["+", "-"] else "+"

    # optionally sort by leftmost start (safer packing). If you want to preserve a
    # provided order (e.g. by read-count), pass tosort=False.
    if tosort:
        data = sorted(data, key=lambda d: min(d[1]))

    rows = []  # each row is a list of isoform entries
    row_intervals = []  # parallel list: for each row, list of (start, end) intervals (absolute coords)

    for iso in data:
        blocksizes, blockstarts = iso[0], iso[1]
        iso_start = min(blockstarts)
        iso_end = blockstarts[-1] + blocksizes[-1]  # last block start + last block size

        label_bp = 0
        if label_padding_bp is not None:
            label_bp = int(max(0, label_padding_bp.get(iso[4], 0)))

        interval_start = iso_start
        interval_end = max(iso_end, iso_start + label_bp)

        placed = False
        # try to place on an existing row
        for r_idx, intervals in enumerate(row_intervals):
            ok = True
            for s, e in intervals:
                # compute minimal distance between [interval_start, interval_end] and [s,e]
                if interval_end < s:
                    gap = s - interval_end
                elif e < interval_start:
                    gap = interval_start - e
                else:
                    gap = -1  # overlap (distance negative)
                if gap < min_spacing:
                    ok = False
                    break
            if ok:
                # place on this row
                rows[r_idx].append(iso)
                intervals.append((interval_start, interval_end))
                placed = True
                break

        if not placed:
            # create new row
            rows.append([iso])
            row_intervals.append([(interval_start, interval_end)])

    return rows


def plot_blocks(data, panel, lower_adj, upper_adj, strand, height=0.5, l=0.8):
    """
    Draw blocks. data contains entries with ABSOLUTE blockstarts.
    We transform to plotting coordinates by subtracting lower_adj.
    This function also sets x-limits with a small left margin to avoid label clipping.
    """
    strand = strand if strand in ["+", "-"] else "+"
    width = float(upper_adj - lower_adj) if (upper_adj - lower_adj) > 0 else 1.0
    left_margin = -args.xmargin_left_frac * width
    right_margin = width + args.xmargin_right_frac * width
    if strand == "-":
        panel.set_xlim(right_margin, left_margin)
    else:
        panel.set_xlim(left_margin, right_margin)

    panel.set_ylim(-0.6, len(data) * 2 - 0.4)
    panel.tick_params(
        axis="both",
        which="both",
        bottom=False,
        labelbottom=False,
        left=False,
        labelleft=False,
        right=False,
        labelright=False,
        top=False,
        labeltop=False,
    )

    di = 0
    ew = 1.5
    for i in range(0, len(data) * 2, 2):
        read = data[di]
        for j in range(len(read)):
            sizes, starts, color, flag, name_full = read[j]
            iso_name = iso_label(name_full)

            # absolute iso coordinates
            iso_start = starts[0]
            iso_end = starts[-1] + sizes[-1]

            # convert to plotting coordinates (relative to lower_adj)
            rel_iso_start = iso_start - lower_adj
            rel_iso_end = iso_end - lower_adj

            # Keep label start on the visual left side regardless of strand.
            # With a reversed x-axis ('-' strand), higher genomic coordinates are left.
            x = rel_iso_end if strand == "-" else rel_iso_start
            ha = "left"

            label_txt = iso_name
            panel.text(
                x,
                i - height / 2 + 1,
                label_txt,
                fontsize=8,
                ha=ha,
                va="center",
                color="black",
                clip_on=False,
                zorder=20,
            )

            for k in range(len(sizes)):
                xrect = starts[k] - lower_adj
                if flag == "REF":
                    rectangle = mplpatches.Rectangle(
                        [xrect, i - height / 2],
                        sizes[k],
                        height,
                        facecolor="none",
                        edgecolor=color,
                        linewidth=ew,
                        zorder=10,
                    )
                    panel.add_patch(rectangle)
                elif flag == "PTC":
                    rectangle = mplpatches.Rectangle(
                        [xrect, i - height / 2],
                        sizes[k],
                        height,
                        facecolor="none",
                        edgecolor=color,
                        linewidth=ew,
                        hatch="////",
                        zorder=10,
                    )
                    panel.add_patch(rectangle)
                    overlay = mplpatches.Rectangle(
                        [xrect, i - height / 2],
                        sizes[k],
                        height,
                        facecolor=color,
                        linewidth=0,
                        zorder=11,
                        alpha=0.2,
                    )
                    panel.add_patch(overlay)
                else:
                    rectangle = mplpatches.Rectangle(
                        [xrect, i - height / 2],
                        sizes[k],
                        height,
                        facecolor=color,
                        edgecolor=color,
                        linewidth=ew,
                        zorder=10,
                    )
                    panel.add_patch(rectangle)
                if k > 0:
                    prev_x = starts[k - 1] + sizes[k - 1] - lower_adj
                    curr_x = starts[k] - lower_adj
                    panel.plot([prev_x, curr_x], [i] * 2, "k-", lw=l, zorder=5)

        di += 1


# ---------- colors & usage ----------
hex_colors = ["#ba748a", "#3498db", "#34495e"]
name_colors = ["windows blue", "faded green", "dusty purple", "amber"]
color_palette = sns.xkcd_palette(name_colors) + sns.color_palette(hex_colors)
gray = sns.xkcd_palette(["greyish"])
if args.palette:
    # read hex colors directly from a palette file (one hex per line)
    color_palette = [line.strip() for line in open(args.palette) if line.strip()]

keepiso = {}
sample_ids = counts_matrix.readline().rstrip().split("\t")[1:]
proportions = []
totals = [0] * len(sample_ids)
gray_bar = [["lowexpr"] + [0] * len(sample_ids) + gray]
for line in counts_matrix:
    line = line.rstrip().split("\t")
    if args.gene_name not in line[0]:
        continue
    counts = [float(x) for x in line[1:]]
    total_counts = sum(counts)
    if total_counts < args.min_reads:
        for i in range(len(sample_ids)):
            gray_bar[0][i + 1] += counts[i]
            totals[i] += counts[i]
        continue
    proportions.append([line[0]] + counts + [total_counts])

proportions = sorted(proportions, key=lambda x: x[-1], reverse=True)
proportions_color = []
for idx, p in enumerate(proportions):
    counts = p[1:-1]
    for i in range(len(sample_ids)):
        totals[i] += counts[i]
    if idx >= len(color_palette):
        for i in range(len(sample_ids)):
            gray_bar[0][i + 1] += counts[i]
        continue
    proportions_color.append(p + [color_palette[idx]])
    keepiso[p[0]] = color_palette[idx]
proportions = [gray_bar[0]] + proportions_color
if len(proportions) == 1:
    sys.stderr.write("No isoforms above threshold\n")
    sys.exit(1)

# Usage plot (stacked bars)
figwidth = 1 + len(sample_ids) * 0.7
plt.figure(figsize=(figwidth, 6))
panel = plt.axes([0.7 / figwidth, 0.11, 1 - 0.7 / figwidth - 0.02, 0.88], frameon=False)
heights = [0] * len(sample_ids)
for iso in proportions:
    for i in range(len(sample_ids)):
        if totals[i] == 0:
            continue
        perc = iso[1 + i] / totals[i] * 100
        rect = mplpatches.Rectangle(
            [i + 1 - 0.4, heights[i]], 0.9, perc, facecolor=iso[-1], linewidth=0
        )
        panel.add_patch(rect)
        if perc >= 7.5:  # add usage percentage
            panel.text(
                i + 1 + 0.07,
                heights[i] + perc / 2,
                str(round(perc, 1)) + "%",
                fontsize=10,
                ha="center",
                va="center",
                color="white",
            )
        if perc >= 12.25:  # add read num
            if iso[1 + i] == 1:
                panel.text(
                    i + 1.45,
                    heights[i] + 1.75,
                    str(int(iso[1 + i])) + " read",
                    fontsize=6,
                    ha="right",
                    va="center",
                    color="white",
                )
            else:
                panel.text(
                    i + 1.44,
                    heights[i] + 1.75,
                    str(int(iso[1 + i])) + " reads",
                    fontsize=6,
                    ha="right",
                    va="center",
                    color="white",
                )

        # Add label if bar big enough

        # if perc > 5:  # threshold in percent
        #    panel.text(
        #        i + 1,
        #        heights[i] + perc / 2,
        #        f"{int(iso[1 + i])} reads",
        #        ha="left",
        #        va="bottom",
        #        fontsize=7,
        #        color="black",
        #        clip_on=False,
        #    )
        heights[i] += perc
panel.set_xlim(0.5, len(sample_ids) + 0.5)
panel.set_ylim(0, 100)
panel.set_xticks(np.arange(1, len(sample_ids) + 1))
panel.set_xticklabels(sample_ids, rotation=20, ha="right")
panel.set_ylabel("Percent Usage", fontsize=12)
plt.savefig(args.o + "_usage.png", dpi=args.img_resolution)

# ---------- isoform structures ----------
# parse sample isoforms (absolute coordinates)
bedfh.seek(0)
isoforms, low_s, up_s, strand, names = parse_bed(bedfh, keepiso=keepiso)
if not isoforms:
    sys.stderr.write("No sample isoforms to plot\n")
    sys.exit(1)

# optional: sort by reads so highest end up on top (we sort ascending so the highest are last -> top)
if args.sort_by_reads:
    iso_counts = {p[0]: p[-2] for p in proportions_color}
    isoforms = sorted(isoforms, key=lambda x: iso_counts.get(x[4], 0))

# parse reference isoforms (absolute coordinates)
ref_isoforms = []
if args.ref_isoforms:
    with open(args.ref_isoforms) as ref_bed:
        ref_all = parse_bed_ref(ref_bed)
    # try to match by gene_name first, otherwise overlap region
    matched_by_name = [r for r in ref_all if args.gene_name in r[4]]
    if matched_by_name:
        ref_isoforms = matched_by_name
    else:
        # select overlapping by absolute coords with sample span
        for r in ref_all:
            r_start = min(r[1])
            r_end = r[1][-1] + r[0][-1]
            if (r_end >= low_s) and (r_start <= up_s):
                ref_isoforms.append(r)
    if not ref_isoforms:
        sys.stderr.write(
            f"Warning: No reference isoforms matched for {args.gene_name}\n"
        )

# decide stacking order: ref_on_top True -> ref above samples (drawn last)
if args.ref_on_top:
    combined = isoforms + ref_isoforms
else:
    combined = ref_isoforms + isoforms

# compute plotting window to cover both sample and references
all_starts = [s for iso in combined for s in iso[1]]
all_ends = [s + size for iso in combined for s, size in zip(iso[1], iso[0])]
pad_bp = (
    max(50, int(0.02 * (max(all_ends) - min(all_starts) + 1))) if all_starts else 50
)
lower_adj = min(all_starts) - pad_bp if all_starts else low_s - pad_bp
upper_adj = max(all_ends) + pad_bp if all_ends else up_s + pad_bp

# dynamic height and bigger left margin for labels
fig_height = max(3, 0.4 * len(combined))
fig_0 = plt.figure(figsize=(9, fig_height))
# increase left margin slightly so long labels aren't clipped by the figure edge
panel = plt.axes([0.03, 0.03, 0.95, 0.94], frameon=True)


# Separate sample and reference isoforms
samples_only = [iso for iso in combined if iso[3] != "REF"]
refs_only = [iso for iso in combined if iso[3] == "REF"]

pack_strand = strand if strand in ["+", "-"] else "+"

locus_width_bp = max(1, int(upper_adj - lower_adj))
label_padding_bp = {
    iso[4]: estimate_label_padding_bp(iso[4], locus_width_bp) for iso in combined
}

# Pack each group independently
packed_samples = pack(
    samples_only,
    tosort=False,
    min_spacing=args.min_isoform_spacing,
    strand=pack_strand,
    label_padding_bp=label_padding_bp,
)
packed_refs = pack(
    refs_only,
    tosort=True,
    min_spacing=args.min_isoform_spacing,
    strand=pack_strand,
    label_padding_bp=label_padding_bp,
)

# Stack depending on --ref_on_top
if args.ref_on_top:
    packed_all = packed_samples + packed_refs
else:
    packed_all = packed_refs + packed_samples

# Pack WITHOUT extra sorting so order/top-bottom is preserved
# packed = pack(combined, tosort=False)
if not packed_all:
    sys.stderr.write("Nothing to plot after packing.\n")
    sys.exit(1)

# prepare iso_to_variant mapping using absolute coordinates (if provided)
iso_to_variant = {}
if args.vcf:
    try:
        with open(args.vcf) as vcf_fh:
            for line in vcf_fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if args.gene_name not in parts[-1]:
                    continue
                # variants: position is absolute (parts[1])
                variant_isos = "".join(parts[-1].split(":")[5:]).split(",")
                pos = int(parts[1])
                alt = parts[4] if len(parts) > 4 else ""
                # add to iso_to_variant if iso present in keepiso mapping
                for k in keepiso:
                    if "".join(k.split(":")) in variant_isos:
                        iso_name = k[: k.rfind("_")] if k.rfind("_") != -1 else k
                        if iso_name not in iso_to_variant:
                            iso_to_variant[iso_name] = []
                        iso_to_variant[iso_name].append(
                            (parts[0], pos, parts[3] if len(parts) > 3 else "", alt)
                        )
    except Exception as e:
        sys.stderr.write(f"Warning: Could not parse VCF: {e}\n")

# finally plot
plot_blocks(packed_all, panel, lower_adj, upper_adj, strand, l=1)

plt.savefig(args.o + "_isoforms.png", dpi=args.img_resolution)

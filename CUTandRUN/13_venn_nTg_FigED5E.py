#!/usr/bin/env python3
# ==============================================================================
# H3K4ac consensus peak overlap, nTg vehicle against nTg RA013915
# Extended Data Fig. 5e
#
# The counts are computed here rather than typed in, because the earlier
# version of this script carried them as literals:
#
#     NTG_VEH_ONLY = 122_922
#     OVERLAP      =  10_272
#     NTG_915_ONLY =  13_636
#
# and those three numbers cannot all be right. The vehicle circle sums to
# 133,194, which is the consensus total, but the treated circle sums to 23,908
# against a consensus total of 23,458. The reason is that bedtools' overlap
# count is not symmetric: the number of vehicle peaks that meet a treated peak
# is not the number of treated peaks that meet a vehicle peak, because peaks
# merge and split differently between the two calls. A two-circle Venn has one
# slot for the intersection and so cannot satisfy both totals at once.
#
# This version computes both intersection counts, draws the diagram with the
# vehicle-side value, and labels each circle with the consensus total it came
# from, so the two numbers inside a circle always add up to the number printed
# with it. The asymmetry is reported on standard output and should be stated in
# the legend rather than hidden by picking one number.
#
# Requires: bedtools on PATH, matplotlib
#
# Input : the two nTg consensus peak BED files
# Output: FigED5e_venn_nTg.pdf / .png, FigED5e_venn_nTg_counts.txt
# ==============================================================================
import os
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

PEAKS = os.environ.get("CUTRUN_PEAKS", "peaks_and_consensus_peaks")
VEH = os.path.join(PEAKS, "nTg-Veh-H3K4ac-consensus.peaks.bed")
D915 = os.path.join(PEAKS, "nTg-RA013915-H3K4ac-consensus.peaks.bed")
for f in (VEH, D915):
    if not os.path.exists(f):
        sys.exit(f"missing: {f}\nSet CUTRUN_PEAKS to the peak directory.")


def nlines(path):
    with open(path) as fh:
        return sum(1 for ln in fh if ln.strip())


def overlap(a, b):
    """Peaks of a that overlap at least one peak of b."""
    p = subprocess.run(["intersectBed", "-a", a, "-b", b, "-u"],
                       capture_output=True, text=True, check=True)
    return sum(1 for ln in p.stdout.splitlines() if ln.strip())


n_veh, n_915 = nlines(VEH), nlines(D915)
ov_veh = overlap(VEH, D915)      # vehicle peaks meeting a treated peak
ov_915 = overlap(D915, VEH)      # treated peaks meeting a vehicle peak
veh_only, d915_only = n_veh - ov_veh, n_915 - ov_915

report = [
    f"nTg Veh       consensus peaks {n_veh:>9,}   overlapping nTg RA013915 {ov_veh:>9,}   only {veh_only:>9,}",
    f"nTg RA013915  consensus peaks {n_915:>9,}   overlapping nTg Veh      {ov_915:>9,}   only {d915_only:>9,}",
    f"the two intersection counts differ by {abs(ov_veh - ov_915):,}; the diagram is drawn",
    "with the vehicle-side value and each circle is labelled with its own total.",
]
print("\n".join(report))
with open("FigED5e_venn_nTg_counts.txt", "w") as fh:
    fh.write("\n".join(report) + "\n")

# ---- the diagram -------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 8))
ax.set_xlim(0, 10)
ax.set_ylim(0, 8)
ax.axis("off")

VEH_COL, D915_COL = "#8C8C8C", "#056AE7"
# The vehicle circle was #CFCFCF, which is the panel's background grey and is
# illegible as a line on white. The darker grey keeps the same hue family.
for fc, ec, xy, w, h in ((VEH_COL, "none", (3.8, 4.0), 5.5, 6.0),
                         (D915_COL, "none", (6.2, 4.0), 3.5, 4.5)):
    ax.add_patch(Ellipse(xy, width=w, height=h, facecolor=fc, alpha=0.2, edgecolor=ec))
for ec, xy, w, h in ((VEH_COL, (3.8, 4.0), 5.5, 6.0),
                     (D915_COL, (6.2, 4.0), 3.5, 4.5)):
    ax.add_patch(Ellipse(xy, width=w, height=h, facecolor="none", edgecolor=ec, linewidth=3.5))

ax.text(1.4, 7.2, f"nTg Veh\n{n_veh:,} peaks", fontsize=13,
        fontweight="bold", color=VEH_COL, va="top")
ax.text(7.4, 7.2, f"nTg RA013915\n{n_915:,} peaks", fontsize=13,
        fontweight="bold", color=D915_COL, va="top")

ax.text(2.4, 4.0, f"{veh_only:,}", fontsize=13, ha="center", va="center")
ax.text(5.2, 4.0, f"{ov_veh:,}", fontsize=12, ha="center", va="center")
ax.text(7.9, 4.0, f"{d915_only:,}", fontsize=13, ha="center", va="center", color=D915_COL)

ax.set_title("H3K4ac consensus peaks\nnTg vehicle vs nTg RA013915", fontsize=16, pad=10)
plt.tight_layout()
plt.savefig("FigED5e_venn_nTg.pdf", dpi=300, bbox_inches="tight")
plt.savefig("FigED5e_venn_nTg.png", dpi=300, bbox_inches="tight")
print("\nwritten: FigED5e_venn_nTg.pdf / .png / _counts.txt")

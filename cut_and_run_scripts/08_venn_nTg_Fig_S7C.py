# ==============================================================================
# Script 08: nTg Vehicle vs nTg RA013915 H3K4ac peak Venn diagram
# Figure: S7C
#
# Overlap counts, generated with bedtools:
#
#   intersectBed -a nTg-Veh-H3K4ac-consensus.peaks.bed -b nTg-RA013915-H3K4ac-consensus.peaks.bed -v | wc -l   # nTg_Veh only
#   intersectBed -a nTg-Veh-H3K4ac-consensus.peaks.bed -b nTg-RA013915-H3K4ac-consensus.peaks.bed -u | wc -l   # overlap
#   intersectBed -a nTg-RA013915-H3K4ac-consensus.peaks.bed -b nTg-Veh-H3K4ac-consensus.peaks.bed -v | wc -l   # nTg_915 only
# ==============================================================================
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

fig, ax = plt.subplots(figsize=(10, 8))
ax.set_xlim(0, 10)
ax.set_ylim(0, 8)
ax.axis('off')

# Transparent fills
e1 = Ellipse((3.8, 4.0), width=5.5, height=6.0, facecolor='#CFCFCF', alpha=0.2, edgecolor='none')
e2 = Ellipse((6.2, 4.0), width=3.5, height=4.5, facecolor='#056AE7', alpha=0.2, edgecolor='none')
# Saturated borders
b1 = Ellipse((3.8, 4.0), width=5.5, height=6.0, facecolor='none', edgecolor='#CFCFCF', linewidth=3.5)
b2 = Ellipse((6.2, 4.0), width=3.5, height=4.5, facecolor='none', edgecolor='#056AE7', linewidth=3.5)
for e in [e1, e2, b1, b2]:
    ax.add_patch(e)

# Labels
ax.text(1.5, 7.2, 'nTg Veh', fontsize=14, fontweight='bold')
ax.text(7.5, 7.2, 'nTg 915', fontsize=14, fontweight='bold', color='#056AE7')

# Counts (see bedtools commands above)
NTG_VEH_ONLY = 122_922
OVERLAP = 10_272
NTG_915_ONLY = 13_636

ax.text(2.5, 4.0, f'{NTG_VEH_ONLY:,}', fontsize=13, ha='center', va='center')
ax.text(5.2, 4.0, f'{OVERLAP:,}', fontsize=12, ha='center', va='center')
ax.text(7.8, 4.0, f'{NTG_915_ONLY:,}', fontsize=13, ha='center', va='center', color='#056AE7')

ax.set_title('H3K4ac Consensus Peaks\nnTg Vehicle vs nTg 915', fontsize=16, pad=10)
plt.tight_layout()
plt.savefig('venn_nTg_relaxed.pdf', dpi=300, bbox_inches='tight')
plt.savefig('venn_nTg_relaxed.png', dpi=300, bbox_inches='tight')
print("Saved!")

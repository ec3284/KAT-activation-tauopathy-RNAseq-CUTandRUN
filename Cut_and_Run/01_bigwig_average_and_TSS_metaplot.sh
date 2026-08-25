#!/bin/bash
# ==============================================================================
# Script 01: BigWig averaging across replicates + TSS metaplot
# Figures: 4A, 4B, S7A
#
# Usage: run from a directory containing per-replicate spike-in-normalized
# bigWig files (same files deposited in GEO) and the mm10 GTF.
# ==============================================================================
set -euo pipefail

WORKDIR="${WORKDIR:-./cutrun_work}"
GTF="${GTF:-${WORKDIR}/Mus_musculus.GRCm38.102.gtf.gz}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

# ---- 1. Average bigWig tracks across biological replicates, per condition ----
bigwigAverage -b nTg-Veh-H3K4ac-rep1.bw nTg-Veh-H3K4ac-rep2.bw nTg-Veh-H3K4ac-rep3.bw \
    -o nTg-Veh-H3K4ac-avg.bw -p max
bigwigAverage -b nTg-RA013915-H3K4ac-rep1.bw nTg-RA013915-H3K4ac-rep2.bw nTg-RA013915-H3K4ac-rep3.bw \
    -o nTg-RA013915-H3K4ac-avg.bw -p max
bigwigAverage -b hTau-Veh-H3K4ac-rep1.bw hTau-Veh-H3K4ac-rep2.bw hTau-Veh-H3K4ac-rep3.bw \
    -o hTau-Veh-H3K4ac-avg.bw -p max
bigwigAverage -b hTau-RA013915-H3K4ac-rep1.bw hTau-RA013915-H3K4ac-rep2.bw hTau-RA013915-H3K4ac-rep3.bw \
    -o hTau-RA013915-H3K4ac-avg.bw -p max
# Output: nTg-Veh-H3K4ac-avg.bw, nTg-RA013915-H3K4ac-avg.bw,
#         hTau-Veh-H3K4ac-avg.bw, hTau-RA013915-H3K4ac-avg.bw

# ---- 2. computeMatrix at TSS +/- 3kb (Fig. 4A-B, S7A) ----
computeMatrix reference-point --referencePoint TSS -b 3000 -a 3000 \
    -R "$GTF" \
    -S nTg-Veh-H3K4ac-avg.bw hTau-Veh-H3K4ac-avg.bw hTau-RA013915-H3K4ac-avg.bw \
    --missingDataAsZero \
    -o matrix_TSS_avg3_zeros.gz -p max

# Metaplot, Figure 4A
plotProfile -m matrix_TSS_avg3_zeros.gz \
    --colors "#CFCFCF" "#FFE0C0" "#FF6A00" \
    --plotTitle "H3K4ac at TSS" \
    --samplesLabel "nTg Veh" "hTau Veh" "hTau 915" \
    --perGroup \
    -o metaplot_TSS_3cond.pdf \
    --plotWidth 12 --plotHeight 8 --plotFileFormat pdf --dpi 300

# ---- 3. computeMatrix at TSS +/- 3kb, all 4 conditions x 3 replicates ----
# (Fig. 4B, S7A: genome-wide promoter H3K4ac AUC per gene, both required
# for the nTg 915 vs nTg Veh comparison in S7A, which needs nTg 915 data
# not used in step 2 above.)
computeMatrix reference-point --referencePoint TSS -b 3000 -a 3000 \
    -R "$GTF" \
    -S WT_Veh_R1.bw WT_Veh_R2.bw WT_Veh_R3.bw \
       WT_915_R1.bw WT_915_R2.bw WT_915_R3.bw \
       hTau_Veh_R1.bw hTau_Veh_R2.bw hTau_Veh_R3.bw \
       hTau_915_R1.bw hTau_915_R2.bw hTau_915_R3.bw \
    --skipZeros \
    -o AUC_PROMOTER_PC.matrix_ALL12_PC.noSkip.gz -p max

# Sample names, in the exact order of the -S list above, one per line
printf "%s\n" \
    "WT_Veh_R1" "WT_Veh_R2" "WT_Veh_R3" \
    "WT_915_R1" "WT_915_R2" "WT_915_R3" \
    "hTau_Veh_R1" "hTau_Veh_R2" "hTau_Veh_R3" \
    "hTau_915_R1" "hTau_915_R2" "hTau_915_R3" \
    > AUC_PROMOTER_PC.samples.txt

# ---- 4. Build transcript_to_gene.tsv from the full GTF ($GTF, step 2) ----
python3 - <<'EOF'
import gzip

GTF = "Mus_musculus.GRCm38.102.gtf.gz"
transcript_to_gene = {}
with gzip.open(GTF, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        if '\ttranscript\t' not in line:
            continue
        if 'transcript_id' not in line:
            continue
        transcript_id = line.split('transcript_id "')[1].split('"')[0]
        gene_name = line.split('gene_name "')[1].split('"')[0] if 'gene_name' in line else 'NA'
        transcript_to_gene[transcript_id] = gene_name

print(f'Trascritti mappati: {len(transcript_to_gene)}')
print('Esempio:', list(transcript_to_gene.items())[:3])

import pandas as pd
pd.Series(transcript_to_gene).to_csv('transcript_to_gene.tsv', sep='\t', header=False)
EOF

# ---- 5. AUC quantification (+/-500bp of TSS) + paired Wilcoxon tests ----
# Figure 4B: nTg Veh vs hTau Veh vs hTau 915
# Figure S7A: nTg Veh vs nTg 915, paired two-sided Wilcoxon signed-rank
#             test, n = 51,352 genes (see legend)
# Uses transcript_to_gene.tsv generated in step 4 above.
python3 - <<'EOF'
import re, math, gzip
import pandas as pd
from scipy.stats import wilcoxon

WORKDIR = "."
OUTPREFIX = f"{WORKDIR}/AUC_PROMOTER_PC"
MATRIX_GZ = f"{WORKDIR}/AUC_PROMOTER_PC.matrix_ALL12_PC.noSkip.gz"
SAMPLES_TXT = f"{WORKDIR}/AUC_PROMOTER_PC.samples.txt"
TRANSCRIPT_MAP = f"{WORKDIR}/transcript_to_gene.tsv"

BIN_SIZE = 10
UPSTREAM = 3000
DOWNSTREAM = 3000
AUC_WINDOW = 500
PSEUDOCOUNT = 1e-6


def group_of(s):
    s = s.replace(" ", "")
    if re.search(r"WT.*Veh", s, re.I): return "nTg_Veh"
    if re.search(r"WT.*915", s, re.I): return "nTg_915"
    if re.search(r"hTau.*Veh", s, re.I): return "hTau_Veh"
    if re.search(r"hTau.*915", s, re.I): return "hTau_915"
    raise ValueError(f"Cannot assign group: {s}")


# Load transcript -> gene mapping
t2g = pd.read_csv(TRANSCRIPT_MAP, sep='\t', header=None, index_col=0)[1].to_dict()

with open(SAMPLES_TXT) as f:
    sample_names = [ln.strip() for ln in f if ln.strip()]

print(f"Samples: {sample_names}")

rows = []
with gzip.open(MATRIX_GZ, "rt") as fh:
    for line in fh:
        if line.startswith("@"):
            continue
        if not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        transcript = p[3] if len(p) > 3 else "NA"
        gene = t2g.get(transcript, transcript)
        vals = []
        for x in p[6:]:
            try:
                vals.append(float(x))
            except ValueError:
                vals.append(float('nan'))
        rows.append((gene, vals))

print(f"Genes parsed: {len(rows)}")

genes = [r[0] for r in rows]
vals_all = [r[1] for r in rows]

total_bp = UPSTREAM + DOWNSTREAM
nbins = total_bp // BIN_SIZE
start_bin = (UPSTREAM - AUC_WINDOW) // BIN_SIZE
end_bin = (UPSTREAM + AUC_WINDOW) // BIN_SIZE

auc_wide = pd.DataFrame({"Gene": genes})
for i, sname in enumerate(sample_names):
    offset = i * nbins
    auc = []
    for v in vals_all:
        window = v[offset + start_bin:offset + end_bin]
        window = [x for x in window if not math.isnan(x)]
        auc.append(sum(window) * BIN_SIZE if window else 0.0)
    auc_wide[sname] = auc

long = auc_wide.melt(id_vars="Gene", var_name="Sample", value_name="AUC")
long["Group"] = long["Sample"].apply(group_of)
long["log2AUC"] = (long["AUC"] + PSEUDOCOUNT).apply(lambda x: math.log2(x))

wide = long.groupby(["Gene", "Group"])["log2AUC"].median().unstack()
prism4 = wide[["nTg_Veh", "nTg_915", "hTau_Veh", "hTau_915"]].dropna()
prism4.to_csv(f"{OUTPREFIX}.PRISM_4conditions.csv")

print(f"Prism rows: {prism4.shape[0]}")

comparisons = [
    ("nTg_Veh_vs_hTau_Veh", ("nTg_Veh", "hTau_Veh")),      # Figure 4B
    ("hTau_Veh_vs_hTau_915", ("hTau_Veh", "hTau_915")),    # Figure 4B
    ("nTg_Veh_vs_hTau_915", ("nTg_Veh", "hTau_915")),      # Figure 4B
    ("nTg_Veh_vs_nTg_915", ("nTg_Veh", "nTg_915")),        # Figure S7A
]

rows_out = []
for name, (a, b) in comparisons:
    tmp = prism4[[a, b]].dropna()
    stat, p = wilcoxon(tmp[a], tmp[b])
    rows_out.append([name, tmp.shape[0], float((tmp[b] - tmp[a]).median()), float(stat), float(p)])

stats_df = pd.DataFrame(rows_out, columns=["Comparison", "N_genes", "Median_delta", "Wilcoxon_stat", "p_value"])
stats_df.to_csv(f"{OUTPREFIX}.WILCOXON_results.csv", index=False)
print(stats_df)
EOF

# The violin plots themselves (Fig. 4B: nTg Veh/hTau Veh/hTau 915; Fig.
# S7A: nTg Veh/nTg 915) were made manually in GraphPad Prism from
# AUC_PROMOTER_PC.PRISM_4conditions.csv, not scripted -- this is expected,
# not a missing step.

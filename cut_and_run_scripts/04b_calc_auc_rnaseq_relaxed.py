# ==============================================================================
# Script 04b: AUC quantification of H3K4ac signal at DEG-up peaks + statistics
# Figure: 4E
# Called by 04_RNAseq_integration_Fig4E.sh
# ==============================================================================
import math
import gzip
import pandas as pd
from scipy.stats import wilcoxon

MATRIX_GZ = "matrix_relaxed_DEGup_3cond.gz"
BIN_SIZE = 10
UPSTREAM = 2000
DOWNSTREAM = 2000
AUC_WINDOW = 500
PSEUDOCOUNT = 1e-6

sample_names = [
    'nTg-Veh-H3K4ac-rep1', 'nTg-Veh-H3K4ac-rep2', 'nTg-Veh-H3K4ac-rep3',
    'hTau-Veh-H3K4ac-rep1', 'hTau-Veh-H3K4ac-rep2', 'hTau-Veh-H3K4ac-rep3',
    'hTau-RA013915-H3K4ac-rep1', 'hTau-RA013915-H3K4ac-rep2', 'hTau-RA013915-H3K4ac-rep3',
]  # must match the order of the -S bigWig list in 04_RNAseq_integration_Fig4E.sh

rows = []
with gzip.open(MATRIX_GZ, "rt") as fh:
    for line in fh:
        if line.startswith("@") or not line.strip():
            continue
        p = line.rstrip("\n").split("\t")
        peak = f"{p[0]}:{p[1]}-{p[2]}"
        vals = []
        for x in p[6:]:
            try:
                vals.append(float(x))
            except ValueError:
                vals.append(0.0)
        rows.append((peak, vals))

print(f"Peaks parsed: {len(rows)}")

nbins = (UPSTREAM + DOWNSTREAM) // BIN_SIZE
start_bin = (UPSTREAM - AUC_WINDOW) // BIN_SIZE
end_bin = (UPSTREAM + AUC_WINDOW) // BIN_SIZE

auc_wide = pd.DataFrame({"Peak": [r[0] for r in rows]})
for i, sname in enumerate(sample_names):
    offset = i * nbins
    auc = [sum(v[offset + start_bin: offset + end_bin]) * BIN_SIZE for _, v in rows]
    auc_wide[sname] = auc

long = auc_wide.melt(id_vars="Peak", var_name="Sample", value_name="AUC")
long["Group"] = long["Sample"].apply(
    lambda s: "nTg_Veh" if "nTg-Veh" in s
    else "hTau_Veh" if "hTau-Veh" in s
    else "hTau_915"
)
long["log2AUC"] = (long["AUC"] + PSEUDOCOUNT).apply(lambda x: math.log2(x))

wide = long.groupby(["Peak", "Group"])["log2AUC"].median().unstack()
prism = wide[["nTg_Veh", "hTau_Veh", "hTau_915"]].dropna()
prism.to_csv("AUC_relaxed_DEGup_3cond.csv")
print(f"Peaks: {prism.shape[0]}")

comparisons = [
    ("nTg_Veh_vs_hTau_Veh", "nTg_Veh", "hTau_Veh"),
    ("hTau_Veh_vs_hTau_915", "hTau_Veh", "hTau_915"),
    ("nTg_Veh_vs_hTau_915", "nTg_Veh", "hTau_915"),
]

rows_out = []
for name, a, b in comparisons:
    stat, p = wilcoxon(prism[a], prism[b])
    delta = (prism[b] - prism[a]).median()
    print(f"\n{name}: delta={delta:.3f}, p={p:.2e}")
    rows_out.append([name, prism.shape[0], delta, stat, p])

pd.DataFrame(
    rows_out,
    columns=["Comparison", "N_peaks", "Median_delta", "Wilcoxon_stat", "p_value"]
).to_csv("WILCOXON_relaxed_DEGup_3cond.csv", index=False)

print("\nDone!")

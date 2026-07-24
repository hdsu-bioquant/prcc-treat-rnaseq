#!/usr/bin/env python3
"""Reproducible Corley-benchmark analysis (backs the numbers in
docs/Pipeline_and_Corley_Explained.md section 5.2).

Inputs (all already produced by the pipeline):
  - results/matrix/gene_counts_matrix.tsv          primary raw counts, 12 libraries
  - results/matrix/three_prime_restricted_matrix.tsv  FL reads in 3' windows
  - results/qc/multiqc_report_data/multiqc_star.txt   STAR mapping stats
  - samples.tsv                                    design (assay, condition, patient)
  - <pipeline>/resources/gdc/gencode.v36.gene_lengths.tsv   gene_id <-> gene_name

Outputs (analysis/tables/):
  qc_mapping_summary.tsv           per-library STAR mapping %
  spearman_correlation_12x12.tsv   full Spearman matrix (rank+Pearson; no scipy)
  correlation_summary.tsv          within-/cross-protocol means (protocol-axis test)
  isg_induction.tsv                per-ISG control->treated fold, both protocols
  isg_signature_per_sample.tsv     mean ISG log2CPM per library
  three_prime_restriction.tsv      FL<->QS agreement: full-gene vs 3'-restricted

Design: 6 patients S1-S6; control = S1/S3/S5, treated = S2/S4/S6; FL<->QS matched by patient.
Run with the snakemake env python (has pandas):  conda run -n snakemake-9.19.0 python corley_analysis.py
"""
import os
import numpy as np
import pandas as pd

# paths derived from this script's own location -> survives moving the pipeline
HERE   = os.path.dirname(os.path.abspath(__file__))   # <pipeline>/Corley_Benchmark/analysis
CORLEY = os.path.dirname(HERE)                        # <pipeline>/Corley_Benchmark
ROOT   = os.path.dirname(CORLEY)                      # <pipeline> root
OUT    = os.path.join(CORLEY, "analysis", "tables")
os.makedirs(OUT, exist_ok=True)

FL = [f"FL_S{i}" for i in range(1, 7)]
QS = [f"QS_S{i}" for i in range(1, 7)]

# ---------- load ----------------------------------------------------------- #
mat = pd.read_csv(os.path.join(CORLEY, "results/matrix/gene_counts_matrix.tsv"),
                  sep="\t", comment="#", index_col=0)
mat = mat[FL + QS]                                        # fix column order
samples = pd.read_csv(os.path.join(CORLEY, "samples.tsv"), sep="\t").set_index("sample")
cond = samples["condition"].to_dict()
gl = pd.read_csv(os.path.join(ROOT, "resources/gdc/gencode.v36.gene_lengths.tsv"), sep="\t")
name2id = {}
for gid, gn in zip(gl.gene_id, gl.gene_name):
    name2id.setdefault(gn, gid)

FL_ctrl = [c for c in FL if cond[c] == "control"]
FL_trt  = [c for c in FL if cond[c] == "treated"]
QS_ctrl = [c for c in QS if cond[c] == "control"]
QS_trt  = [c for c in QS if cond[c] == "treated"]

def cpm(df):
    return df / df.sum(axis=0) * 1e6

CPM = cpm(mat)

def spearman_vec(a, b):
    """Spearman rho of two 1-D arrays via rank + Pearson (no scipy)."""
    ra = pd.Series(np.asarray(a, float)).rank().to_numpy()
    rb = pd.Series(np.asarray(b, float)).rank().to_numpy()
    return float(np.corrcoef(ra, rb)[0, 1])

# ---------- 1. QC mapping summary ------------------------------------------ #
star = pd.read_csv(os.path.join(CORLEY, "results/qc/multiqc_report_data/multiqc_star.txt"), sep="\t")
star = star[~star["Sample"].astype(str).str.contains("_STARpass1")]     # keep final (pass-2) stats only
star["Sample"] = star["Sample"].str.replace(r"\.Log\.final$", "", regex=True).str.replace(r"\.$", "", regex=True)
keep = [c for c in ["Sample", "total_reads", "uniquely_mapped", "uniquely_mapped_percent",
                    "multimapped_percent", "unmapped_tooshort_percent"] if c in star.columns]
qc = star[keep].copy()
qc.to_csv(os.path.join(OUT, "qc_mapping_summary.tsv"), sep="\t", index=False)

# ---------- 2. Spearman correlation matrix (expressed genes) --------------- #
expressed = mat[mat.sum(axis=1) > 0]                     # drop all-zero genes
ranks = expressed.rank()                                 # rank each library across genes
spear = ranks.corr()                                     # Pearson on ranks == Spearman
spear.to_csv(os.path.join(OUT, "spearman_correlation_12x12.tsv"), sep="\t")

def mean_offdiag(block):
    v = block.to_numpy()
    iu = np.triu_indices_from(v, k=1)
    return float(v[iu].mean())

within_fl = mean_offdiag(spear.loc[FL, FL])
within_qs = mean_offdiag(spear.loc[QS, QS])
cross_matched = float(np.mean([spear.loc[f"FL_S{i}", f"QS_S{i}"] for i in range(1, 7)]))
cross_block = spear.loc[FL, QS].to_numpy()
matched_mask = np.eye(6, dtype=bool)
cross_unmatched = float(cross_block[~matched_mask].mean())

# nearest-neighbour test: is each FL sample's top correlate another FL sample?
nn_same_protocol = 0
for s in FL + QS:
    row = spear[s].drop(s)
    top = row.idxmax()
    same = (top in FL) == (s in FL)
    nn_same_protocol += int(same)

corr_summary = pd.DataFrame([
    ["within_FL_mean",        within_fl,      "mean Spearman among the 6 full-length libraries"],
    ["within_QS_mean",        within_qs,      "mean Spearman among the 6 QuantSeq libraries"],
    ["cross_matched_mean",    cross_matched,  "mean Spearman of FL_Sx vs its own QS_Sx (6 pairs)"],
    ["cross_unmatched_mean",  cross_unmatched,"mean Spearman of FL_Sx vs QS_Sy, x!=y (30 pairs)"],
    ["nn_same_protocol_frac", nn_same_protocol/12.0,
     "fraction of libraries whose nearest correlate is the SAME protocol (1.0 = protocol dominates)"],
], columns=["metric", "value", "description"])
corr_summary.to_csv(os.path.join(OUT, "correlation_summary.tsv"), sep="\t", index=False)

# ---------- 3. ISG signature (Poly(I:C) response) -------------------------- #
ISG = ["ISG15", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "IFIT1", "IFIT2", "IFIT3",
       "IFI6", "IFI27", "IFI44", "IFI44L", "RSAD2", "USP18", "IFI35", "IRF7", "STAT1",
       "STAT2", "DDX58", "IFIH1", "HERC5", "HERC6", "XAF1", "CMPK2", "EPSTI1", "LY6E",
       "BST2", "CXCL10"]
isg_ids = {g: name2id[g] for g in ISG if g in name2id and name2id[g] in CPM.index}

rows = []
for g, gid in isg_ids.items():
    fl_c, fl_t = CPM.loc[gid, FL_ctrl].mean(), CPM.loc[gid, FL_trt].mean()
    qs_c, qs_t = CPM.loc[gid, QS_ctrl].mean(), CPM.loc[gid, QS_trt].mean()
    rows.append({
        "gene": g, "gene_id": gid,
        "FL_ctrl_cpm": round(fl_c, 2), "FL_trt_cpm": round(fl_t, 2),
        "FL_fold": round((fl_t + 1) / (fl_c + 1), 2),
        "QS_ctrl_cpm": round(qs_c, 2), "QS_trt_cpm": round(qs_t, 2),
        "QS_fold": round((qs_t + 1) / (qs_c + 1), 2),
    })
isg = pd.DataFrame(rows).sort_values("FL_fold", ascending=False)
isg.to_csv(os.path.join(OUT, "isg_induction.tsv"), sep="\t", index=False)

# per-sample signature score = mean log2(CPM+1) over the ISG panel
logcpm = np.log2(CPM.loc[list(isg_ids.values())] + 1)
sig = pd.DataFrame({
    "sample": logcpm.columns,
    "assay": [samples.loc[s, "assay"] for s in logcpm.columns],
    "condition": [cond[s] for s in logcpm.columns],
    "isg_score_mean_log2cpm": logcpm.mean(axis=0).round(3).to_numpy(),
})
sig.to_csv(os.path.join(OUT, "isg_signature_per_sample.tsv"), sep="\t", index=False)

isg_fold_med_fl = float(isg["FL_fold"].median())
isg_fold_med_qs = float(isg["QS_fold"].median())

# ---------- 4. 3'-restriction: does it help FL<->QS agreement? ------------- #
threep = pd.read_csv(os.path.join(CORLEY, "results/matrix/three_prime_restricted_matrix.tsv"),
                     sep="\t", index_col=0)
threep = threep[[c for c in FL if c in threep.columns]]
common = mat.index.intersection(threep.index)
# genes expressed somewhere in the compared libraries (avoid all-zero ties)
expr = mat.loc[common, FL + QS].sum(axis=1) > 0
common = common[expr]

fl_full = cpm(mat.loc[common, FL])
qs_full = cpm(mat.loc[common, QS])
fl_3p   = cpm(threep.loc[common, FL])

rows = []
for i in range(1, 7):
    full = spearman_vec(fl_full[f"FL_S{i}"], qs_full[f"QS_S{i}"])
    thr  = spearman_vec(fl_3p[f"FL_S{i}"],   qs_full[f"QS_S{i}"])
    rows.append({"patient": f"S{i}", "condition": cond[f"FL_S{i}"],
                 "spearman_fullgene_FLvsQS": round(full, 3),
                 "spearman_3prime_FLvsQS": round(thr, 3),
                 "delta_3prime_minus_full": round(thr - full, 3)})
tp = pd.DataFrame(rows)
tp.loc[len(tp)] = ["MEAN", "", round(tp["spearman_fullgene_FLvsQS"].mean(), 3),
                   round(tp["spearman_3prime_FLvsQS"].mean(), 3),
                   round(tp["delta_3prime_minus_full"].mean(), 3)]
tp.to_csv(os.path.join(OUT, "three_prime_restriction.tsv"), sep="\t", index=False)

# ---------- console summary ------------------------------------------------ #
print("=== Corley reproducible analysis ===")
print(f"genes in matrix: {mat.shape[0]}   expressed (sum>0): {expressed.shape[0]}")
print(f"[correlation] within-FL={within_fl:.3f}  within-QS={within_qs:.3f}  "
      f"cross-matched={cross_matched:.3f}  cross-unmatched={cross_unmatched:.3f}")
print(f"[correlation] nearest-neighbour same-protocol: {nn_same_protocol}/12 "
      f"({nn_same_protocol/12*100:.0f}%)  -> protocol {'dominates' if nn_same_protocol==12 else 'partial'}")
print(f"[ISG] genes found: {len(isg_ids)}/{len(ISG)}  median fold  FL={isg_fold_med_fl:.1f}x  QS={isg_fold_med_qs:.1f}x")
print(f"[ISG] signature score control vs treated (FL): "
      f"{sig[(sig.assay.str.startswith('full')) & (sig.condition=='control')].isg_score_mean_log2cpm.mean():.2f}"
      f" -> {sig[(sig.assay.str.startswith('full')) & (sig.condition=='treated')].isg_score_mean_log2cpm.mean():.2f} (log2CPM)")
print(f"[3'-restriction] FL<->QS  full-gene={tp.iloc[-1]['spearman_fullgene_FLvsQS']}  "
      f"3'-restricted={tp.iloc[-1]['spearman_3prime_FLvsQS']}  "
      f"delta={tp.iloc[-1]['delta_3prime_minus_full']} (negative => 3'-restriction HURTS)")
print(f"tables -> {OUT}")

# ============================================================================
# 07_motif_combination_frequency.py
# ROLE   : A combinatorial (UpSet-style) motif frequency figure for an
#          8-area subset drawn entirely from this pipeline's alternate
#          PM(7)/AT(7) area definition (see ANALYSIS.md, common.py): RSP,
#          VIS, PTLp (PM); PERI, AI, TEa, ORB, PFC (AT). This mirrors a
#          similar combinatorial-motif figure done for the manuscript's
#          primary PM(7)/AT(13) analysis (a separate pipeline, not included
#          in this repo), whose AT picks for that figure (ACA, SS, MO) are
#          not modeled at all under this alternate definition, so they are
#          replaced here with PERI, AI, and TEa: three AT areas not already
#          used elsewhere in this figure's PM picks, chosen to span this
#          group's prevalence range (PERI 0.706 lowest, AI 0.731 mid, TEa
#          0.831 highest of the candidates) rather than clustering near one
#          end of it.
# INPUT  : data/master_detailed_comment.csv
# OUTPUT : outputs/motif_combination_frequency.{png,pdf}
# FIGURE : supplementary / exploratory, standalone in this repo
# NOTES  : PFC here is one of this definition's own 14 modeled areas (a
#          native `PFC_length`/`PFC_endpoint` column, not a derived PL+ILA
#          combination the model fits separately) -- so no OR-merge
#          special-casing is needed for either real or simulated data; all 8
#          areas are read the same generic way.
#          Duplicates script 04's refined-model constants/machinery rather
#          than importing script 04 (no __main__ guard, so importing it would
#          re-run its whole analysis and write files as an import side
#          effect). Re-sync from `fit_refined_params.py` if the model
#          is re-fit.
# ============================================================================

"""
motif_combination_frequency.py

Eight areas: RSP, VIS, PTLp (PM network), PERI, AI, TEa, ORB (AT network),
and PFC (AT network, a genuine modeled area under this definition -- see
NOTES). For each neuron, its "motif" is the exact 8-bit present/absent
pattern across these areas, restricted to just these 8 columns (the other 6
modeled areas are ignored for this classification). Real data's most
frequent motifs are compared against the refined model's predicted count for
those SAME motifs, from 2000 Monte Carlo replicate 160-neuron datasets
(median + 95% percentile band).

Layout (one figure, one shared x-axis of motifs, ordered by real-data
frequency descending):
  Panel 1 (top)    -- real data: neuron count per motif (plain bars).
  Panel 2 (middle) -- refined model: median predicted count + 95% band.
  Panel 3 (bottom) -- combination matrix: which of the 8 areas each motif
                      involves (filled dot = present, colored by PM/AT
                      network; a vertical line connects a motif's present
                      areas, the standard UpSet-plot convention).
"""

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.hermite import hermgauss

from common import (
    OUT_DIR, A, N, ispm, AREAS, sig, logit,
    Y as Yr, t, a as a0,
)

# =====================================================================
# 1. Refined-model machinery, duplicated from 04_refined_model_fit.py
#    (see module NOTES for why)
# =====================================================================
R_prev = Yr.mean(0)
A_PM, K_PM, C_PM = 5.869, 4.000, -0.009
A_AT, K_AT, C_AT = 8.316, 4.000, -0.421
SIGMA            = 1.513

_GHx, _GHw = hermgauss(16); _GH_W = _GHw / np.sqrt(np.pi)

def _gradient():
    f_pm = A_PM * sig(-K_PM * (t - C_PM))
    f_at = A_AT * sig( K_AT * (t - C_AT))
    return np.where(ispm[None, :], f_pm[:, None], f_at[:, None])
G_REF = _gradient()

def recalibrated_a():
    a = a0.copy()
    for _ in range(12):
        pmean = np.zeros(A)
        for q in range(len(_GHx)):
            pmean += _GH_W[q] * sig(a[None, :] + G_REF + np.sqrt(2) * SIGMA * _GHx[q]).mean(0)
        a += (logit(np.clip(R_prev, 0.02, 0.98)) - logit(np.clip(pmean, 0.02, 0.98))) * 0.7
    return a
A_CAL = recalibrated_a()

_check_pmean = np.zeros(A)
for _q in range(len(_GHx)):
    _check_pmean += _GH_W[_q] * sig(A_CAL[None, :] + G_REF + np.sqrt(2) * SIGMA * _GHx[_q]).mean(0)
assert np.max(np.abs(_check_pmean - R_prev)) < 0.01, (
    "Refined-model constants (duplicated from 04_refined_model_fit.py) no "
    "longer reproduce real per-area prevalence -- re-sync with script 04 and "
    "re-run scripts/fit_refined_params.py."
)

def simulate(rng):
    """One stochastic replicate Y (N,A) binary projection matrix, refined model."""
    delta = rng.normal(0.0, SIGMA, N)
    e = sig(A_CAL[None, :] + G_REF + delta[:, None])
    return (rng.random((N, A)) < e).astype(float)

# =====================================================================
# 2. The 8-area motif definition
# =====================================================================
# All 8 are ordinary modeled areas under this alternate definition -- unlike
# script 16, PFC needs no OR-merge special-casing (see module NOTES).
LABELS8 = ["RSP", "VIS", "PTLp", "PERI", "AI", "TEa", "ORB", "PFC"]
IDX8    = [AREAS.index(nm) for nm in LABELS8]
IS_PM8  = np.array([ispm[j] for j in IDX8])
POW2    = 2 ** np.arange(8)

def motif_codes(Y8):
    return (Y8.astype(int) * POW2).sum(1)

# =====================================================================
# 3. Real data: rank the 20 most frequent motifs
# =====================================================================
Y8_real = Yr[:, IDX8]
codes_real = motif_codes(Y8_real)
uniq, counts = np.unique(codes_real, return_counts=True)
# stable sort so ties (equal counts) break by ascending code -- reproducible,
# not an arbitrary/undocumented tie rule
order = np.argsort(-counts, kind="stable")
TOP_N = 20
top_codes = uniq[order][:TOP_N]
top_counts_real = counts[order][:TOP_N]

print(f"{len(uniq)} distinct motifs occur in real data (of 256 possible), N={N}")
print(f"Top {TOP_N} motifs by real-data frequency:")
for code, c in zip(top_codes, top_counts_real):
    bits = [(code >> k) & 1 for k in range(8)]
    present = "+".join(l for l, b in zip(LABELS8, bits) if b) or "(none)"
    print(f"  count={c:3d}  code={code:3d}  -> {present}")
print(f"(20th-place count={top_counts_real[-1]}; 21st-place count="
      f"{counts[order][TOP_N] if len(order) > TOP_N else 'n/a'})")

# =====================================================================
# 4. Refined model: Monte Carlo predicted count for the SAME 20 motifs
# =====================================================================
NREP = 2000
rng = np.random.default_rng(7)
per_rep_counts = np.zeros((NREP, TOP_N))
code_to_col = {code: i for i, code in enumerate(top_codes)}
for r in range(NREP):
    Ys = simulate(rng)
    codes_sim = motif_codes(Ys[:, IDX8])
    sim_uniq, sim_counts = np.unique(codes_sim, return_counts=True)
    for code, c in zip(sim_uniq, sim_counts):
        col = code_to_col.get(code)
        if col is not None:
            per_rep_counts[r, col] = c

band = np.percentile(per_rep_counts, [2.5, 50, 97.5], axis=0)   # (3, TOP_N)
print(f"\nRefined model ({NREP} replicates), median [2.5-97.5%] predicted count for the same {TOP_N} motifs:")
for i, code in enumerate(top_codes):
    bits = [(code >> k) & 1 for k in range(8)]
    present = "+".join(l for l, b in zip(LABELS8, bits) if b) or "(none)"
    print(f"  real={top_counts_real[i]:3d}  model=[{band[0,i]:.1f}, {band[1,i]:.1f}, {band[2,i]:.1f}]  -> {present}")

# =====================================================================
# 5. Figure
# =====================================================================
plt.rcParams.update({"font.family": "sans-serif", "font.size": 8, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5})

REAL, MOD, PMc, ATc = "#1a1a1a", "#2166AC", "#FF7800", "#3485DC"
FAINT = "#d9d9d9"

x = np.arange(TOP_N)
ymax = max(top_counts_real.max(), band[2].max()) * 1.12

fig = plt.figure(figsize=(9.5, 6.6))
gs = fig.add_gridspec(3, 1, height_ratios=[1.3, 1.3, 1.15], hspace=0.12)
ax_real, ax_mod, ax_mat = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])

# --- Panel 1: real data ---
ax_real.bar(x, top_counts_real, width=0.62, color=REAL)
ax_real.set_ylabel("Neuron count\n(real data)", fontsize=8)
ax_real.set_ylim(0, ymax)
ax_real.set_xticks(x); ax_real.set_xticklabels([])
ax_real.spines[["top", "right"]].set_visible(False)
ax_real.set_title(f"Top {TOP_N} motifs (of {len(uniq)} observed, 256 possible): "
                   f"real data vs. refined-model prediction (alt PM/AT definition)", fontsize=9.5)

# --- Panel 2: refined model, median + 95% band ---
ax_mod.bar(x, band[1], width=0.62, color=MOD, alpha=0.85,
           yerr=[band[1] - band[0], band[2] - band[1]], capsize=2.5,
           error_kw=dict(lw=0.9, ecolor="#0d3b66"))
ax_mod.set_ylabel("Neuron count\n(refined model,\nmedian & 95%)", fontsize=8)
ax_mod.set_ylim(0, ymax)
ax_mod.set_xticks(x); ax_mod.set_xticklabels([])
ax_mod.spines[["top", "right"]].set_visible(False)

# --- Panel 3: combination matrix (UpSet-style) ---
row_y = np.arange(8)
for col, code in enumerate(top_codes):
    bits = [(code >> k) & 1 for k in range(8)]
    present_rows = [ry for ry, b in zip(row_y, bits) if b]
    if len(present_rows) > 1:
        ax_mat.plot([col, col], [min(present_rows), max(present_rows)],
                   color="#4d4d4d", lw=1.4, zorder=1)
    for ry, b in zip(row_y, bits):
        color = (PMc if IS_PM8[ry] else ATc) if b else FAINT
        ax_mat.scatter(col, ry, s=48 if b else 26, color=color,
                       edgecolors="none" if b else "#bdbdbd", linewidths=0.6, zorder=2)

ax_mat.set_xlim(-0.6, TOP_N - 0.4)
ax_mat.set_ylim(-0.6, 7.6)
ax_mat.invert_yaxis()
ax_mat.set_yticks(row_y)
ax_mat.set_yticklabels(LABELS8, fontsize=8)
for tick, is_pm in zip(ax_mat.get_yticklabels(), IS_PM8):
    tick.set_color(PMc if is_pm else ATc)
ax_mat.set_xticks(x)
ax_mat.set_xticklabels([str(i + 1) for i in x], fontsize=7)
ax_mat.set_xlabel("Motif rank (by real-data frequency)", fontsize=8)
ax_mat.spines[["top", "right", "left"]].set_visible(False)
ax_mat.tick_params(left=False)

fig.align_ylabels([ax_real, ax_mod])
fig.savefig(OUT_DIR / "motif_combination_frequency.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT_DIR / "motif_combination_frequency.pdf", bbox_inches="tight")
print("\nSaved motif_combination_frequency.png / .pdf")

# ============================================================================
# 01_sim_models_fitted.py
# ROLE   : Fits 11 candidate models' free shape parameters to the real data
#          under this pipeline's alternate PM(7)/AT(7) area definition (see
#          ANALYSIS.md for the area lists and rationale), then simulates and
#          plots each fitted model's predictions. Mirrors the fit/simulate/
#          plot logic used by the manuscript's primary PM(7)/AT(13) analysis
#          (a separate pipeline, not included in this repo), applied here to
#          a smaller, more symmetric alternate area grouping.
# INPUT  : data/master_detailed_comment.csv  (160-neuron PMAT subset)
# OUTPUT : outputs/sim_models_fitted_grid.{png,pdf}
# NOTES  : This script's printed "4. Graded identity (1 axis)" fitted kappa
#          is the source for script 03's hardcoded KAPPA constant -- run this
#          before finalizing 03. Fits to aggregated PM/AT probability curves
#          (monotonic) -- see ANALYSIS.md for why this suppresses the
#          broadcast bump (script 03 fixes).
# ============================================================================

"""
sim_models_fitted.py

Fits each of the 11 models' free shape parameters to the real 160-neuron
PMAT-filtered data (14-area alternate definition), then simulates and plots
each fitted model's predictions. NO cross-model comparison/scoring is done
-- fit -> simulate -> plot, one model at a time.

FITTING METHOD
--------------
Per-area empirical intercepts a_j (logit of marginal projection prevalence)
are fixed exactly as before -- not re-fit here.

For each model we define a deterministic (no-RNG) prediction of two curves:
    meanPM(t)  = average projection probability across the 7 PM areas
    meanAT(t)  = average projection probability across the 7 AT areas
as a function of a small number of model-specific free shape parameters.
These are fit by least squares against smooth empirical target curves from
an aggregated binomial GLM of real PM/AT projection counts on real ML
position (logit link, statsmodels).

Models with zero identity-defining free parameters (null, area-intrinsic)
are not fit; they use a_j (and, for null, only the grand-mean intercept)
directly.

Structural/nuisance choices NOT part of the fit: m_geometry's random area
physical positions, m_discrete's intermediate-zone shift (fixed at 0),
m_choosegrow's growth-stage rates, and the shared Gamma length parameters
(LMAX, KSHAPE) -- the curve fit targets projection probability only.

BOUNDARY / DEGENERACY FIXES
----------------------------
(1) Developmental gradient: fitting scale_d and steepness independently
    is non-identifiable on these two curves -- as scale_d -> inf and
    steepness -> 0 with their product held fixed, the tanh curve becomes
    a straight line, i.e. the model's steep->0 limit IS the graded model.
    That produced the earlier boundary hit. Fixed by pinning the model's
    amplitude at t=+-1 to match the fitted graded model's amplitude
    (kappa_graded), leaving only the saturation steepness as a single
    free, identifiable parameter.

(2) Graded identity vs. choose-then-grow are exactly degenerate on every
    projection-probability statistic (choose-then-grow's growth stage
    never touches selection). Resolved with a genuinely discriminating
    statistic: conditional axon length (given the neuron projects to that
    area), regressed on sign-corrected ML position with area fixed
    effects (Gamma GLM, log link). Graded ties growth to the same
    engagement that drives selection -> nonzero slope; choose-then-grow's
    growth stage is area-random and ML-independent -> slope ~ 0. Reported
    as an explicit sixth column for every row, including real data.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.optimize import minimize, minimize_scalar
import statsmodels.api as sm

from common import (
    OUT_DIR, sub, PM, AT, AREAS, A, ispm, N, nPM, nAT,
    sig, logit, Y as Y_real, L as L_real, t, a, rsp_idx, orb_idx, other_idx,
    LMAX, KSHAPE,
)

rng = np.random.default_rng(7)

# =====================================================================
# 1. Real data
# =====================================================================
a_pm, a_at = a[ispm], a[~ispm]

# outcomes(e): stochastically draws whether a neuron projects to each area
# (Bernoulli(e), e = engagement probability) and, if so, its axon length
# (Gamma, mean set by elen -- defaults to e itself). This is the shared
# noise model every candidate model plugs its own engagement signal into.
def outcomes(e, elen=None):
    e = np.clip(e, 1e-4, 1 - 1e-4)
    if elen is None: elen = e
    Y = (rng.random(e.shape) < e).astype(float)
    mean = LMAX * np.clip(elen, 1e-4, None)
    L = rng.gamma(KSHAPE, mean / KSHAPE) * Y
    return Y, L

# =====================================================================
# 2. Empirical fitting targets: smooth PM/AT probability curves vs t
# =====================================================================
# The curve every model below is optimized to reproduce: mean PM / mean AT
# projection probability as a function of ML position, smoothed by a
# binomial GLM. This is a NETWORK-IDENTITY target (mirrors script 02's
# CV target) -- it says nothing about breadth or the broadcast phenotype,
# which is why fits here don't guarantee a match on those (see script 02).
pm_sum = Y_real[:, ispm].sum(axis=1)
at_sum = Y_real[:, ~ispm].sum(axis=1)
X = sm.add_constant(t)
glm_pm = sm.GLM(np.column_stack([pm_sum, nPM - pm_sum]), X, family=sm.families.Binomial()).fit()
glm_at = sm.GLM(np.column_stack([at_sum, nAT - at_sum]), X, family=sm.families.Binomial()).fit()
target_pm = glm_pm.predict(X)
target_at = glm_at.predict(X)

def loss(pred_pm, pred_at):
    return np.sum((pred_pm - target_pm)**2) + np.sum((pred_at - target_at)**2)

# =====================================================================
# 3. Per-model deterministic mean-PM(t)/mean-AT(t) predictors
# =====================================================================
# H: a single symmetric ML gradient sets PM-vs-AT engagement -- one
# continuous "identity" signal drives both curves in opposite directions.
def curves_graded(theta):
    kappa, = theta
    lam = -kappa * t
    return sig(a_pm.mean() + lam), sig(a_at.mean() - lam)

# H: neurons are discrete PM / AT / intermediate cell types (hard zones
# set by kappa_d, c), not points on a continuous gradient; "shift" sets
# each zone's identity strength.
def curves_discrete(theta):
    kappa_d, c, shift = theta
    zpm = sig(kappa_d * (-t - c)); zat = sig(kappa_d * (t - c))
    zint = np.clip(1 - zpm - zat, 0, 1)
    tot = zpm + zat + zint
    zpm, zat, zint = zpm/tot, zat/tot, zint/tot
    pm = zpm*sig(a_pm.mean()+shift) + zat*sig(a_pm.mean()-shift) + zint*sig(a_pm.mean())
    at = zat*sig(a_at.mean()+shift) + zpm*sig(a_at.mean()-shift) + zint*sig(a_at.mean())
    return pm, at

# H: each neuron allocates a fixed "budget" of N0 targets, traded off
# between PM and AT by position (f(t)) -- an allocation constraint, not an
# independent per-area probability.
def curves_tradeoff(theta):
    kappa_tr, N0 = theta
    f = sig(-kappa_tr * t)
    return np.clip(f*N0, 0, nPM)/nPM, np.clip(N0-f*N0, 0, nAT)/nAT

# H: same graded PM/AT gradient as curves_graded, plus a symmetric
# mid-axis "boost" (amp, width) meant to capture broadcast neurons -- but
# fitting only to these monotonic curves pushes the boost toward ~0 (see
# ANALYSIS.md / script 02's PPC for why this alone can't fix the mismatch).
def curves_broadcast(theta):
    kappa_b, amp, width = theta
    lam = -kappa_b * t
    boost = amp * np.exp(-t**2 / (2*width**2))
    return sig(a_pm.mean()+lam+boost), sig(a_at.mean()-lam+boost)

# H: identical selection curve to curves_graded (this and 4 are exactly
# degenerate on projection PROBABILITY) -- they differ only in axon-length
# growth (see m_choosegrow below and the conditional_length_slope
# discriminator further down).
def curves_choosegrow(theta):
    kappa_cg, = theta
    lam = -kappa_cg * t
    return sig(a_pm.mean()+lam), sig(a_at.mean()-lam)

# H: PM and AT engagement are two independent drives, each with its own
# slope/offset, rather than mirror images of one shared gradient.
def curves_twodrive(theta):
    slope_pm, off_pm, slope_at, off_at = theta
    return sig(slope_pm*(off_pm - t)), sig(slope_at*(t + off_at))

# H: neurons first choose a network (piPM/piAT, position-driven), then
# each area's probability within that network is set by the area's own
# intrinsic preference (scaled a) -- network choice and area choice are
# separate, hierarchical decisions.
def curves_hierarchical(theta):
    kappa_h, = theta
    piPM = sig(kappa_h*(-t)); piAT = sig(kappa_h*t)
    def scale(v): return 0.5 + 0.5*(v - v.min())/(v.max()-v.min()+1e-9)
    return piPM*scale(a_pm).mean(), piAT*scale(a_at).mean()

# H: projection probability follows physical proximity to each area's
# fixed anatomical position (pa_pm_mean/pa_at_mean) along the ML axis -- a
# wiring-cost/geometric account rather than a biological identity signal.
def curves_geometry(theta):
    coef, offset = theta
    pa_pm_mean, pa_at_mean = -0.65, 0.65
    d_pm = np.abs(t - pa_pm_mean); d_at = np.abs(t - pa_at_mean)
    return sig(a_pm.mean()+coef*(offset-d_pm)), sig(a_at.mean()+coef*(offset-d_at))

# H: PM/AT engagement saturates (tanh) rather than growing linearly with
# position, as if a developmental signal plateaus; amplitude is pinned to
# the fitted graded model's so only steepness is free (see BOUNDARY/
# DEGENERACY FIXES (1) in the module docstring).
def curves_developmental_matched(steep, lam_max):
    scale_d = lam_max / np.tanh(steep)
    lam2 = scale_d * np.tanh(steep * (-t))
    return sig(a_pm.mean()+lam2), sig(a_at.mean()-lam2)

SPECS = {
    "3. Discrete cell types":      (curves_discrete,     [3.0, 0.33, 1.7],   [(0.1,10),(-1,1),(0,4)]),
    "4. Graded identity (1 axis)": (curves_graded,       [2.5],              [(0.1, 10)]),
    "4a. Trade-off allocation":    (curves_tradeoff,     [2.5, 10.0],        [(0.1,10),(1,nPM+nAT)]),
    "4b. Broadcast / expanded":    (curves_broadcast,    [2.5, 2.5, 0.35],   [(0.1,10),(0,6),(0.05,1.5)]),
    "5. Choose-then-grow":         (curves_choosegrow,   [2.5],              [(0.1, 10)]),
    "6. Two independent drives":   (curves_twodrive,     [2.5, 0.3, 2.2, 0.3], [(0.1,10),(-1,1),(0.1,10),(-1,1)]),
    "7. Network-then-area":        (curves_hierarchical, [2.5],              [(0.1, 10)]),
    "8. Wiring cost / geometry":   (curves_geometry,     [4.0, 0.7],         [(0.1,10),(0.1,1.5)]),
}

fitted = {}
print(f"{'model':32s} {'fitted params'}")
for name, (fn, x0, bounds) in SPECS.items():
    res = minimize(lambda th: loss(*fn(th)), x0, bounds=bounds, method="L-BFGS-B")
    fitted[name] = res.x
    print(f"{name:32s} {np.round(res.x, 3)}   (loss {res.fun:.4f})")

lam_max = fitted["4. Graded identity (1 axis)"][0] * 1.0
res_dev = minimize_scalar(lambda s: loss(*curves_developmental_matched(s, lam_max)),
                           bounds=(0.05, 8), method="bounded")
scale_dev = lam_max / np.tanh(res_dev.x)
fitted["9. Developmental gradient"] = np.array([scale_dev, res_dev.x])
print(f"{'9. Developmental gradient (matched)':32s} "
      f"[scale={scale_dev:.3f} steep={res_dev.x:.3f}]   (loss {res_dev.fun:.4f})  "
      f"-- steep near lower bound => indistinguishable from linear/graded on these curves")

# =====================================================================
# 4. Full stochastic simulation per model using FITTED parameters
# =====================================================================
# No positional info at all: identical average engagement everywhere on
# the ML axis -- the "is there any position signal?" floor.
def m_null():
    return outcomes(sig(np.full((N, A), a.mean())))

# No positional info: each area keeps its own intrinsic prevalence, same
# for every neuron regardless of where it sits -- "areas just differ in
# popularity" floor.
def m_area():
    return outcomes(sig(np.tile(a, (N, 1))))

# Discrete cell types, simulated: each neuron is hard-assigned to a PM,
# AT, or intermediate zone by position, but with stochastic
# "misassignment" noise near the zone boundary (p_reassign) representing
# ambiguity at the boundary rather than a perfectly sharp cutoff.
def m_discrete(theta):
    kappa_d, c, shift = theta
    zone = np.where(t < -c, 0, np.where(t > c, 2, 1))
    p_reassign = 1 - np.abs(sig(kappa_d*(-t-c)) - sig(kappa_d*(t-c)))
    mis = rng.random(N) < np.clip(p_reassign*0.3, 0, 0.5)
    zone[mis] = rng.integers(0, 3, mis.sum())
    term = np.zeros((N, A))
    term[zone == 0] = np.where(ispm, shift, -shift)[None, :]
    term[zone == 2] = np.where(ispm, -shift, shift)[None, :]
    term[zone == 1] = 0.0
    return outcomes(sig(a[None, :] + term))

# Graded identity, simulated: one continuous engagement signal (lam) drives
# PM up and AT down (or vice versa) in lockstep -- the reference model
# every degeneracy check in this script is stated relative to.
def m_graded(theta):
    kappa, = theta
    lam = -kappa * t
    term = np.where(ispm[None, :], lam[:, None], -lam[:, None])
    return outcomes(sig(a[None, :] + term))

# Trade-off allocation, simulated: each neuron gets a hard integer budget
# of N0 targets, split between networks by f(t), and WHICH areas within
# each network get chosen is deterministic (rank by baseline preference)
# -- a fixed-budget allocation process rather than independent per-area
# probabilities.
def m_tradeoff(theta):
    kappa_tr, N0 = theta
    f = sig(-kappa_tr * t)
    qpm, qat = a_pm.copy(), a_at.copy()
    pm_idx = np.where(ispm)[0]; at_idx = np.where(~ispm)[0]
    Y = np.zeros((N, A)); L = np.zeros((N, A))
    for i in range(N):
        npm = int(np.clip(round(f[i]*N0), 0, nPM))
        nat = int(np.clip(N0 - npm, 0, nAT))
        cp = pm_idx[np.argsort(-qpm)[:npm]]
        ca = at_idx[np.argsort(-qat)[:nat]]
        if npm: L[i, cp] = f[i]*50.0/npm
        if nat: L[i, ca] = (1-f[i])*50.0/nat
        Y[i, np.r_[cp, ca]] = 1
    L = L * rng.gamma(KSHAPE, 1/KSHAPE, (N, A))
    return Y, L

# Broadcast/expanded, simulated: graded gradient plus the mid-axis boost
# -- if this boost survives fitting, neurons near the ML midpoint should
# show elevated breadth toward both networks at once (columns 2-4 below
# test whether it actually does).
def m_broadcast(theta):
    kappa_b, amp, width = theta
    lam = -kappa_b * t
    boost = amp * np.exp(-t**2 / (2*width**2))
    term = np.where(ispm[None, :], lam[:, None]+boost[:, None], -lam[:, None]+boost[:, None])
    return outcomes(sig(a[None, :] + term))

# Choose-then-grow, simulated: identical selection probability (e) to
# m_graded, but axon length is drawn using an area-random, ML-independent
# growth stage (mua) instead of e itself -- the point of this model is
# that growth is decoupled from the position-driven selection signal.
def m_choosegrow(theta):
    kappa_cg, = theta
    lam = -kappa_cg * t
    term = np.where(ispm[None, :], lam[:, None], -lam[:, None])
    e = sig(a[None, :] + term)
    mua = rng.uniform(0.4, 0.6, A)   # growth stage: not fit, area-random, ML-independent
    elen = np.tile(mua, (N, 1))
    return outcomes(e, elen=elen)

# Two independent drives, simulated: PM and AT engagement each follow
# their own fitted slope/offset rather than a single shared gradient.
def m_twodrive(theta):
    slope_pm, off_pm, slope_at, off_at = theta
    u = sig(slope_pm*(off_pm - t)); v = sig(slope_at*(t + off_at))
    term_u = logit(np.clip(u, 1e-4, 1-1e-4)) - a_pm.mean()
    term_v = logit(np.clip(v, 1e-4, 1-1e-4)) - a_at.mean()
    term = np.where(ispm[None, :], term_u[:, None], term_v[:, None])
    return outcomes(sig(a[None, :] + term))

# Network-then-area, simulated: position sets network choice (piPM/piAT),
# then each area's probability within that network is scaled by its own
# intrinsic preference (qpm/qat) -- a two-stage, hierarchical decision.
def m_hierarchical(theta):
    kappa_h, = theta
    piPM = sig(kappa_h*(-t)); piAT = sig(kappa_h*t)
    def scale(v): return 0.5 + 0.5*(v - v.min())/(v.max()-v.min()+1e-9)
    qpm, qat = scale(a_pm), scale(a_at)
    e = np.zeros((N, A))
    e[:, ispm] = piPM[:, None] * qpm[None, :]
    e[:, ~ispm] = piAT[:, None] * qat[None, :]
    return outcomes(e)

# Wiring cost/geometry, simulated: each area gets its own random fixed
# anatomical position along the ML axis (pa), and projection probability
# falls off with physical distance from it -- a spatial-proximity account
# with no notion of "identity" at all.
def m_geometry(theta):
    coef, offset = theta
    pa = np.where(ispm, rng.uniform(-0.9, -0.4, A), rng.uniform(0.4, 0.9, A))
    d = np.abs(t[:, None] - pa[None, :])
    term = coef*(offset - d) - a[None, :].mean()
    return outcomes(sig(a[None, :] + term))

# Developmental gradient, simulated: same shape as m_graded but saturating
# (tanh) rather than linear in t, as if guided by a developmental signal
# that plateaus toward the poles instead of growing without bound.
def m_developmental(theta):
    scale_d, steep = theta
    lam2 = scale_d * np.tanh(steep*(-t))
    term = np.where(ispm[None, :], lam2[:, None], -lam2[:, None])
    return outcomes(sig(a[None, :] + term))

ROWS = [
    ("REAL DATA",                    lambda: (Y_real, L_real)),
    ("1. Homogeneous (null)",        m_null),
    ("2. Area-intrinsic",            m_area),
    ("3. Discrete cell types",       lambda: m_discrete(fitted["3. Discrete cell types"])),
    ("4. Graded identity (1 axis)",  lambda: m_graded(fitted["4. Graded identity (1 axis)"])),
    ("4a. Trade-off allocation",     lambda: m_tradeoff(fitted["4a. Trade-off allocation"])),
    ("4b. Broadcast / expanded",     lambda: m_broadcast(fitted["4b. Broadcast / expanded"])),
    ("5. Choose-then-grow",          lambda: m_choosegrow(fitted["5. Choose-then-grow"])),
    ("6. Two independent drives",    lambda: m_twodrive(fitted["6. Two independent drives"])),
    ("7. Network-then-area",         lambda: m_hierarchical(fitted["7. Network-then-area"])),
    ("8. Wiring cost / geometry",    lambda: m_geometry(fitted["8. Wiring cost / geometry"])),
    ("9. Developmental gradient",    lambda: m_developmental(fitted["9. Developmental gradient"])),
]

# =====================================================================
# 5. Summary helpers, including the new conditional-length discriminator
# =====================================================================
# RSP (a PM area) and ORB (an AT area) as flagship network-engagement
# markers: splits neurons into PM-only, AT-only, and "both" (broadcast
# candidates) -- the same bridging-neuron grouping used in script 02.
def group_labels(Y):
    rsp = Y[:, rsp_idx] > 0; orb = Y[:, orb_idx] > 0
    g = np.full(N, "", dtype=object)
    g[rsp & ~orb] = "RSP+ ORB-"; g[orb & ~rsp] = "ORB+ RSP-"; g[rsp & orb] = "RSP+ ORB+"
    return g

# Per-neuron PM/AT identity: signed difference in how many PM vs AT areas
# a neuron projects to. This is the quantity whose relationship with ML
# position (rho) is the headline "position organizes identity" result.
def bias_score(Y):
    pm = Y[:, ispm].sum(axis=1); at = Y[:, ~ispm].sum(axis=1)
    return (pm - at) / np.clip(pm + at, 1, None)

def bin_means(x, y, nbins=6):
    edges = np.linspace(x.min(), x.max(), nbins+1)
    idx = np.clip(np.digitize(x, edges[1:-1]), 0, nbins-1)
    cx, cm, cs = [], [], []
    for b in range(nbins):
        sel = idx == b
        if sel.sum() == 0: continue
        cx.append(x[sel].mean()); cm.append(y[sel].mean())
        cs.append(y[sel].std(ddof=1)/np.sqrt(sel.sum()) if sel.sum() > 1 else 0.0)
    return np.array(cx), np.array(cm), np.array(cs)

def conditional_length_slope(Y, L):
    """Gamma GLM (log link) of conditional log-length on sign-corrected
    ML position, with area fixed effects. Returns (coef, se) or (nan, nan)
    if too few conditional observations to fit."""
    rows = []
    for j, area in enumerate(AREAS):
        proj = Y[:, j] > 0
        if proj.sum() < 5: continue
        signed_t = t[proj] * (1 if not ispm[j] else -1)
        for st, ln in zip(signed_t, L[proj, j]):
            if ln > 0:
                rows.append((area, st, np.log(ln)))
    if len(rows) < 30:
        return np.nan, np.nan
    gdf = pd.DataFrame(rows, columns=["area", "signed_t", "logL"])
    dummies = pd.get_dummies(gdf["area"], drop_first=True).astype(float)
    Xg = sm.add_constant(pd.concat([gdf["signed_t"], dummies], axis=1))
    try:
        with np.errstate(all="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fit = sm.GLM(gdf["logL"], Xg, family=sm.families.Gamma(link=sm.families.links.Log())).fit()
        c, s = fit.params["signed_t"], fit.bse["signed_t"]
        if not np.isfinite(c):
            return np.nan, np.nan
        return c, s
    except Exception:
        return np.nan, np.nan

def plot_length_slope_panel(ax, coef, se):
    """Draws the conditional-log-length-vs-signed-ML-slope panel from a
    (coef, se) pair returned by conditional_length_slope(). Kept for reuse
    (e.g. a standalone supplemental panel) -- not called in the saved grid
    below, which only shows columns 1-5."""
    if np.isnan(coef):
        ax.text(0.5, 0.5, "n/a\n(too few obs)", ha="center", va="center", fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])
    else:
        ax.errorbar([0], [coef], yerr=[1.96*se], fmt="o", color="#55A868", ms=8, capsize=4)
        ax.axhline(0, color="grey", lw=0.5)
        ax.set_xlim(-1, 1); ax.set_xticks([])
        ax.set_ylabel("conditional log-length\nslope (signed ML)")

# =====================================================================
# 6. Plot: five columns (conditional-length column computed/printed below
#    but not plotted -- see plot_length_slope_panel() to re-enable it)
# =====================================================================
fig, axes = plt.subplots(len(ROWS), 5, figsize=(19, 3.0*len(ROWS)))
groups = ["RSP+ ORB-", "ORB+ RSP-", "RSP+ ORB+"]

print("\nconditional log-length ~ signed ML slope (area fixed effects):")
for row, (name, fn) in enumerate(ROWS):
    Y, L = fn()
    b = bias_score(Y); br = Y.sum(axis=1)
    rho, _ = spearmanr(t, b)
    ax0, ax1, ax2, ax3, ax4 = axes[row]

    # Col 1: does ML position organize PM/AT identity? (the headline
    # rho ~ -0.83 result -- real data's per-neuron scatter plus rank corr.)
    ax0.scatter(t, b, s=8, alpha=0.35, color="#4C72B0")
    ax0.set_title(f"{name}\nrho={rho:.2f}", fontsize=8)
    ax0.set_xlabel("ML (t)"); ax0.set_ylabel("PM/AT bias"); ax0.axhline(0, color="grey", lw=0.5)

    # Col 2: does total breadth vary with position? Real data peaks
    # mid-axis (a hump); check whether each model reproduces that shape or
    # only a monotonic trend (script 02's rho PPC quantifies this).
    cx, cm, cs = bin_means(t, br, 6)
    ax1.errorbar(cx, cm, yerr=cs, fmt="o-", color="#C44E52", ms=4, lw=1)
    ax1.set_xlabel("ML (binned)"); ax1.set_ylabel("# areas (all 14)")

    # Col 3: mean per-area projection probability by network -- the
    # crossing PM/AT curves every model's curves_* function was fit to.
    cx2, pmm, pms = bin_means(t, Y[:, ispm].mean(axis=1), 6)
    _, atm, ats = bin_means(t, Y[:, ~ispm].mean(axis=1), 6)
    ax2.errorbar(cx2, pmm, yerr=pms, fmt="o-", color="#FF7800", ms=4, lw=1, label="PM")
    ax2.errorbar(cx2, atm, yerr=ats, fmt="o-", color="#3485DC", ms=4, lw=1, label="AT")
    ax2.set_xlabel("ML (binned)"); ax2.set_ylabel("mean proj. prob.")
    if row == 0: ax2.legend(fontsize=6, frameon=False)

    # Cols 4-5: the broadcast-phenotype group contrast, per model -- do
    # neurons bridging both network anchors (RSP+ORB) also project more
    # widely / with more total length across the other 12 areas, versus
    # single-network specialists?
    g = group_labels(Y)
    breadth_other = Y[:, other_idx].sum(axis=1)
    length_other = L[:, other_idx].sum(axis=1)
    ax3.boxplot([breadth_other[g==grp] for grp in groups], tick_labels=[x.replace(" ","\n") for x in groups])
    ax3.set_ylabel("# other targets (12)"); ax3.tick_params(axis="x", labelsize=6)
    ax4.boxplot([length_other[g==grp] for grp in groups], tick_labels=[x.replace(" ","\n") for x in groups])
    ax4.set_ylabel("total length, other 12"); ax4.tick_params(axis="x", labelsize=6)

    # Graded-vs-choose-then-grow degeneracy breaker (see BOUNDARY/
    # DEGENERACY FIXES (2) above): graded ties growth to the same signal
    # that drives selection (expect nonzero slope); choose-then-grow's
    # growth stage is ML-independent (expect slope ~ 0).
    coef, se = conditional_length_slope(Y, L)
    print(f"  {name:32s} coef={coef if np.isnan(coef) else round(coef,3)}"
          f"{'' if np.isnan(coef) else f'  se={se:.3f}'}")

plt.tight_layout()
plt.savefig(OUT_DIR / "sim_models_fitted_grid.png", dpi=140, bbox_inches="tight")
plt.savefig(OUT_DIR / "sim_models_fitted_grid.pdf", bbox_inches="tight")
print("\nSaved sim_models_fitted_grid.png / .pdf")

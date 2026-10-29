# ============================================================================
# 02_sim_models_comparison.py
# ROLE   : Compares the 11 candidate models fit by script 01, under this
#          pipeline's alternate PM(7)/AT(7) area definition (see
#          ANALYSIS.md). Mirrors the comparison logic used by the manuscript's
#          primary PM(7)/AT(13) analysis (a separate pipeline, not included
#          in this repo).
# INPUT  : data/master_detailed_comment.csv
# OUTPUT : outputs/sim_models_comparison.{png,pdf}
# NOTES  : CV LL is an aggregate (mean-PM/mean-AT) pseudo-likelihood -- ranks
#          curve shape, not per-area structure.
# ============================================================================

"""
sim_models_comparison.py

Compares the 11 models using two complementary strategies (from the
earlier "how should we compare fits" discussion):

(A) 5-FOLD CROSS-VALIDATED LOG-LIKELIHOOD (primary quantitative score)
    For each fold: intercepts a_j AND each model's shape parameters are
    RE-FIT on the training neurons only (80%), then evaluated as held-out
    binomial log-likelihood on the untouched test neurons (20%). This is
    done at the same aggregate level the curves were originally fit at --
    PM count out of 7 ~ Binomial(7, p_PM(t)), AT count out of 13 ~
    Binomial(13, p_AT(t)) -- so scoring is on exactly what was fit, and
    extra free parameters are only rewarded if they generalize, not just
    fit training noise. Null/area-intrinsic have no shape parameters to
    refit (trivial "fit").

(B) POSTERIOR-PREDICTIVE CHECKS on statistics NEVER used in fitting
    Using each model's FULL-DATA fitted parameters, 300 stochastic
    replicate datasets are simulated and used to build a predictive
    distribution for:
      1. Spearman rho(breadth, ML)      -- the ~22-vs-~13 breadth signature
      2. RSP/ORB "both" vs "single" group breadth contrast
    The real value's percentile within each model's predictive
    distribution is reported: percentiles near 50 indicate the model's
    predictive spread comfortably contains the real effect size;
    percentiles near 0/100 indicate the model over/under-shoots it.
    Conditional-length slope (graded vs. choose-then-grow discriminator)
    is carried over as a single-realization point estimate with 95% CI,
    as established previously.

NO single number is treated as the final word -- (A) is the primary
generalization score, (B) is qualitative discriminating evidence on
independent statistics.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.optimize import minimize, minimize_scalar
import statsmodels.api as sm
import warnings

from common import (
    OUT_DIR, PM, AT, AREAS, A, ispm, N, nPM, nAT,
    sig, logit, Y as Y_real, L as L_real, t as t_all, a as a_full_base,
    rsp_idx, orb_idx, other_idx, LMAX, KSHAPE,
)

rng_global = np.random.default_rng(7)

# =====================================================================
# 1. Real data
# =====================================================================
# How many of each neuron's 7 PM / 13 AT areas it projects to -- the
# aggregate readout of NETWORK ENGAGEMENT that every model below is fit
# and CV-scored against. This is about identity (PM vs AT), not breadth.
pm_sum_all = Y_real[:, ispm].sum(axis=1)
at_sum_all = Y_real[:, ~ispm].sum(axis=1)

def outcomes(e, elen=None, rng=rng_global):
    e = np.clip(e, 1e-4, 1 - 1e-4)
    if elen is None: elen = e
    Y = (rng.random(e.shape) < e).astype(float)
    mean = LMAX * np.clip(elen, 1e-4, None)
    L = rng.gamma(KSHAPE, mean / KSHAPE) * Y
    return Y, L

# =====================================================================
# 2. Curve functions, parametrized by (theta, a_pm_mean, a_at_mean, t)
# so they can be refit/evaluated on arbitrary train/test splits.
# =====================================================================
# H: a single symmetric ML gradient sets PM-vs-AT engagement -- one
# continuous "identity" signal drives both curves in opposite directions.
def curves_graded(theta, apm, aat, tt):
    kappa, = theta
    lam = -kappa * tt
    return sig(apm + lam), sig(aat - lam)

# H: neurons are discrete PM / AT / intermediate cell types (hard zones
# set by kappa_d, c), not points on a continuous gradient; "shift" sets
# each zone's identity strength.
def curves_discrete(theta, apm, aat, tt):
    kappa_d, c, shift = theta
    zpm = sig(kappa_d * (-tt - c)); zat = sig(kappa_d * (tt - c))
    zint = np.clip(1 - zpm - zat, 0, 1)
    tot = zpm + zat + zint
    zpm, zat, zint = zpm/tot, zat/tot, zint/tot
    pm = zpm*sig(apm+shift) + zat*sig(apm-shift) + zint*sig(apm)
    at = zat*sig(aat+shift) + zpm*sig(aat-shift) + zint*sig(aat)
    return pm, at

# H: each neuron allocates a fixed "budget" of N0 targets, traded off
# between PM and AT by position (f(t)) -- an allocation constraint, not an
# independent per-area probability.
def curves_tradeoff(theta, apm, aat, tt):
    kappa_tr, N0 = theta
    f = sig(-kappa_tr * tt)
    return np.clip(f*N0, 0, nPM)/nPM, np.clip(N0-f*N0, 0, nAT)/nAT

# H: same graded PM/AT gradient as curves_graded, plus a symmetric
# mid-axis "boost" (amp, width) meant to capture broadcast neurons -- the
# boost this model needs is exactly what the PPC in section 5 tests for.
def curves_broadcast(theta, apm, aat, tt):
    kappa_b, amp, width = theta
    lam = -kappa_b * tt
    boost = amp * np.exp(-tt**2 / (2*width**2))
    return sig(apm+lam+boost), sig(aat-lam+boost)

# H: identical selection curve to curves_graded (this and 4 are exactly
# degenerate on projection PROBABILITY) -- they differ only in axon-length
# growth, which is what script 01's conditional-length-slope discriminates.
def curves_choosegrow(theta, apm, aat, tt):
    kappa_cg, = theta
    lam = -kappa_cg * tt
    return sig(apm+lam), sig(aat-lam)

# H: PM and AT engagement are two independent drives, each with its own
# slope/offset, rather than mirror images of one shared gradient.
def curves_twodrive(theta, apm, aat, tt):
    slope_pm, off_pm, slope_at, off_at = theta
    return sig(slope_pm*(off_pm - tt)), sig(slope_at*(tt + off_at))

# H: neurons first choose a network (piPM/piAT, position-driven), then
# each area's probability within that network is set by the area's own
# intrinsic preference (QPM_BAR/QAT_BAR) -- network choice and area choice
# are separate, hierarchical decisions.
def curves_hierarchical(theta, apm, aat, tt):
    kappa_h, = theta
    piPM = sig(kappa_h*(-tt)); piAT = sig(kappa_h*tt)
    # within-network preference (scaled a) doesn't depend on the ML split,
    # so it's held at its full-data value for simplicity across folds
    return piPM*QPM_BAR, piAT*QAT_BAR

# H: projection probability follows physical proximity to each area's
# fixed anatomical position (pa_pm_mean/pa_at_mean) along the ML axis -- a
# wiring-cost/geometric account rather than a biological identity signal.
def curves_geometry(theta, apm, aat, tt):
    coef, offset = theta
    pa_pm_mean, pa_at_mean = -0.65, 0.65
    d_pm = np.abs(tt - pa_pm_mean); d_at = np.abs(tt - pa_at_mean)
    return sig(apm+coef*(offset-d_pm)), sig(aat+coef*(offset-d_at))

# H: PM/AT engagement saturates (tanh) rather than growing linearly with
# position, as if a developmental signal plateaus; steepness is the only
# free parameter once amplitude is pinned to the graded model's (see
# script 01 NOTES on why scale/steepness aren't jointly identifiable).
def curves_developmental(theta, apm, aat, tt):
    scale_d, steep = theta
    lam2 = scale_d * np.tanh(steep * (-tt))
    return sig(apm+lam2), sig(aat-lam2)

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
MODEL_NAMES = ["1. Homogeneous (null)", "2. Area-intrinsic"] + list(SPECS.keys()) + ["9. Developmental gradient"]

# fixed nuisance scaling for hierarchical (doesn't depend on train/test split)
a_full = a_full_base
def _scale(v): return 0.5 + 0.5*(v - v.min())/(v.max()-v.min()+1e-9)
QPM_BAR = _scale(a_full[ispm]).mean()
QAT_BAR = _scale(a_full[~ispm]).mean()

def binom_loglik(y_succ, n_trials, p):
    p = np.clip(p, 1e-6, 1-1e-6)
    return np.sum(y_succ*np.log(p) + (n_trials-y_succ)*np.log(1-p))

def fit_curve_params(fn, x0, bounds, apm, aat, tt, pm_target, at_target):
    def obj(th):
        pm, at = fn(th, apm, aat, tt)
        return np.sum((pm-pm_target)**2) + np.sum((at-at_target)**2)
    res = minimize(obj, x0, bounds=bounds, method="L-BFGS-B")
    return res.x

# =====================================================================
# 3. 5-fold CV log-likelihood
# =====================================================================
# Biological question: does ML position generalize, out-of-sample, as a
# predictor of NETWORK IDENTITY (PM vs AT engagement)? This section says
# nothing about breadth or the broadcast phenotype -- that's section 5.
n_folds = 5
fold_id = np.tile(np.arange(n_folds), N//n_folds + 1)[:N]
rng_fold = np.random.default_rng(0)
rng_fold.shuffle(fold_id)

cv_ll = {name: 0.0 for name in MODEL_NAMES}
in_sample_loss = {}

for k in range(n_folds):
    train = fold_id != k
    test = fold_id == k
    t_tr, t_te = t_all[train], t_all[test]

    p_hat_tr = np.clip(Y_real[train].mean(axis=0), 0.02, 0.98)
    a_tr = logit(p_hat_tr)
    apm_tr, aat_tr = a_tr[ispm].mean(), a_tr[~ispm].mean()

    pm_tr = Y_real[train][:, ispm].sum(axis=1); at_tr = Y_real[train][:, ~ispm].sum(axis=1)
    X_tr = sm.add_constant(t_tr)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        g_pm = sm.GLM(np.column_stack([pm_tr, nPM-pm_tr]), X_tr, family=sm.families.Binomial()).fit()
        g_at = sm.GLM(np.column_stack([at_tr, nAT-at_tr]), X_tr, family=sm.families.Binomial()).fit()
    pm_target_tr = g_pm.predict(X_tr); at_target_tr = g_at.predict(X_tr)

    pm_te = Y_real[test][:, ispm].sum(axis=1); at_te = Y_real[test][:, ~ispm].sum(axis=1)

    # null: constant p = train grand-mean intercept (same for PM and AT)
    # -- is there ANY position signal at all, or does one number explain
    # engagement equally well everywhere along the ML axis?
    p_null = sig(a_tr.mean())
    cv_ll["1. Homogeneous (null)"] += binom_loglik(pm_te, nPM, p_null) + binom_loglik(at_te, nAT, p_null)

    # area-intrinsic: constant per-network mean of per-area sigmoids (no t
    # dependence) -- do areas just differ in overall popularity, with no
    # role for the neuron's position at all?
    p_area_pm = sig(a_tr[ispm]).mean(); p_area_at = sig(a_tr[~ispm]).mean()
    cv_ll["2. Area-intrinsic"] += binom_loglik(pm_te, nPM, p_area_pm) + binom_loglik(at_te, nAT, p_area_at)

    for name, (fn, x0, bounds) in SPECS.items():
        theta = fit_curve_params(fn, x0, bounds, apm_tr, aat_tr, t_tr, pm_target_tr, at_target_tr)
        pm_pred, at_pred = fn(theta, apm_tr, aat_tr, t_te)
        cv_ll[name] += binom_loglik(pm_te, nPM, pm_pred) + binom_loglik(at_te, nAT, at_pred)

    # developmental: amplitude-matched to this fold's graded fit
    kappa_g = fit_curve_params(curves_graded, [2.5], [(0.1,10)], apm_tr, aat_tr, t_tr, pm_target_tr, at_target_tr)[0]
    lam_max_tr = kappa_g * 1.0
    res_dev = minimize_scalar(
        lambda s: np.sum((curves_developmental([lam_max_tr/np.tanh(s), s], apm_tr, aat_tr, t_tr)[0]-pm_target_tr)**2)
                + np.sum((curves_developmental([lam_max_tr/np.tanh(s), s], apm_tr, aat_tr, t_tr)[1]-at_target_tr)**2),
        bounds=(0.05, 8), method="bounded")
    theta_dev = [lam_max_tr/np.tanh(res_dev.x), res_dev.x]
    pm_pred, at_pred = curves_developmental(theta_dev, apm_tr, aat_tr, t_te)
    cv_ll["9. Developmental gradient"] += binom_loglik(pm_te, nPM, pm_pred) + binom_loglik(at_te, nAT, at_pred)

print("5-fold CV held-out log-likelihood (higher = better generalization):")
for name in sorted(cv_ll, key=cv_ll.get, reverse=True):
    print(f"  {name:32s} {cv_ll[name]:9.2f}")

# =====================================================================
# 4. Full-data fits (for posterior predictive checks)
# =====================================================================
X_all = sm.add_constant(t_all)
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    glm_pm = sm.GLM(np.column_stack([pm_sum_all, nPM-pm_sum_all]), X_all, family=sm.families.Binomial()).fit()
    glm_at = sm.GLM(np.column_stack([at_sum_all, nAT-at_sum_all]), X_all, family=sm.families.Binomial()).fit()
target_pm_all, target_at_all = glm_pm.predict(X_all), glm_at.predict(X_all)
apm_full, aat_full = a_full[ispm].mean(), a_full[~ispm].mean()

fitted_full = {}
for name, (fn, x0, bounds) in SPECS.items():
    fitted_full[name] = fit_curve_params(fn, x0, bounds, apm_full, aat_full, t_all, target_pm_all, target_at_all)
kappa_g_full = fitted_full["4. Graded identity (1 axis)"][0]
lam_max_full = kappa_g_full * 1.0
res_dev_full = minimize_scalar(
    lambda s: np.sum((curves_developmental([lam_max_full/np.tanh(s), s], apm_full, aat_full, t_all)[0]-target_pm_all)**2)
            + np.sum((curves_developmental([lam_max_full/np.tanh(s), s], apm_full, aat_full, t_all)[1]-target_at_all)**2),
    bounds=(0.05, 8), method="bounded")
fitted_full["9. Developmental gradient"] = np.array([lam_max_full/np.tanh(res_dev_full.x), res_dev_full.x])

# full stochastic simulators (area-level), reusing fitted_full theta
def sim_null(rng): return outcomes(sig(np.full((N, A), a_full.mean())), rng=rng)
def sim_area(rng): return outcomes(sig(np.tile(a_full, (N, 1))), rng=rng)

def sim_discrete(rng):
    kappa_d, c, shift = fitted_full["3. Discrete cell types"]
    zone = np.where(t_all < -c, 0, np.where(t_all > c, 2, 1))
    p_reassign = 1 - np.abs(sig(kappa_d*(-t_all-c)) - sig(kappa_d*(t_all-c)))
    mis = rng.random(N) < np.clip(p_reassign*0.3, 0, 0.5)
    zone[mis] = rng.integers(0, 3, mis.sum())
    term = np.zeros((N, A))
    term[zone == 0] = np.where(ispm, shift, -shift)[None, :]
    term[zone == 2] = np.where(ispm, -shift, shift)[None, :]
    term[zone == 1] = 0.0
    return outcomes(sig(a_full[None, :] + term), rng=rng)

def sim_graded(rng):
    kappa, = fitted_full["4. Graded identity (1 axis)"]
    lam = -kappa * t_all
    term = np.where(ispm[None, :], lam[:, None], -lam[:, None])
    return outcomes(sig(a_full[None, :] + term), rng=rng)

def sim_tradeoff(rng):
    kappa_tr, N0 = fitted_full["4a. Trade-off allocation"]
    f = sig(-kappa_tr * t_all)
    qpm, qat = a_full[ispm].copy(), a_full[~ispm].copy()
    pm_idx = np.where(ispm)[0]; at_idx = np.where(~ispm)[0]
    Y = np.zeros((N, A)); L = np.zeros((N, A))
    for i in range(N):
        npm = int(np.clip(round(f[i]*N0), 0, nPM))
        nat = int(np.clip(N0 - npm, 0, nAT))
        cp = pm_idx[np.argsort(-qpm)[:npm]]; ca = at_idx[np.argsort(-qat)[:nat]]
        if npm: L[i, cp] = f[i]*50.0/npm
        if nat: L[i, ca] = (1-f[i])*50.0/nat
        Y[i, np.r_[cp, ca]] = 1
    L = L * rng.gamma(KSHAPE, 1/KSHAPE, (N, A))
    return Y, L

def sim_broadcast(rng):
    kappa_b, amp, width = fitted_full["4b. Broadcast / expanded"]
    lam = -kappa_b * t_all
    boost = amp * np.exp(-t_all**2 / (2*width**2))
    term = np.where(ispm[None, :], lam[:, None]+boost[:, None], -lam[:, None]+boost[:, None])
    return outcomes(sig(a_full[None, :] + term), rng=rng)

def sim_choosegrow(rng):
    kappa_cg, = fitted_full["5. Choose-then-grow"]
    lam = -kappa_cg * t_all
    term = np.where(ispm[None, :], lam[:, None], -lam[:, None])
    e = sig(a_full[None, :] + term)
    mua = rng.uniform(0.4, 0.6, A)
    elen = np.tile(mua, (N, 1))
    return outcomes(e, elen=elen, rng=rng)

def sim_twodrive(rng):
    slope_pm, off_pm, slope_at, off_at = fitted_full["6. Two independent drives"]
    u = sig(slope_pm*(off_pm - t_all)); v = sig(slope_at*(t_all + off_at))
    term_u = logit(np.clip(u, 1e-4, 1-1e-4)) - a_full[ispm].mean()
    term_v = logit(np.clip(v, 1e-4, 1-1e-4)) - a_full[~ispm].mean()
    term = np.where(ispm[None, :], term_u[:, None], term_v[:, None])
    return outcomes(sig(a_full[None, :] + term), rng=rng)

def sim_hierarchical(rng):
    kappa_h, = fitted_full["7. Network-then-area"]
    piPM = sig(kappa_h*(-t_all)); piAT = sig(kappa_h*t_all)
    qpm, qat = _scale(a_full[ispm]), _scale(a_full[~ispm])
    e = np.zeros((N, A))
    e[:, ispm] = piPM[:, None] * qpm[None, :]
    e[:, ~ispm] = piAT[:, None] * qat[None, :]
    return outcomes(e, rng=rng)

def sim_geometry(rng):
    coef, offset = fitted_full["8. Wiring cost / geometry"]
    pa = np.where(ispm, rng.uniform(-0.9, -0.4, A), rng.uniform(0.4, 0.9, A))
    d = np.abs(t_all[:, None] - pa[None, :])
    term = coef*(offset - d) - a_full[None, :].mean()
    return outcomes(sig(a_full[None, :] + term), rng=rng)

def sim_developmental(rng):
    scale_d, steep = fitted_full["9. Developmental gradient"]
    lam2 = scale_d * np.tanh(steep*(-t_all))
    term = np.where(ispm[None, :], lam2[:, None], -lam2[:, None])
    return outcomes(sig(a_full[None, :] + term), rng=rng)

SIMULATORS = {
    "1. Homogeneous (null)": sim_null, "2. Area-intrinsic": sim_area,
    "3. Discrete cell types": sim_discrete, "4. Graded identity (1 axis)": sim_graded,
    "4a. Trade-off allocation": sim_tradeoff, "4b. Broadcast / expanded": sim_broadcast,
    "5. Choose-then-grow": sim_choosegrow, "6. Two independent drives": sim_twodrive,
    "7. Network-then-area": sim_hierarchical, "8. Wiring cost / geometry": sim_geometry,
    "9. Developmental gradient": sim_developmental,
}

# =====================================================================
# 5. Posterior predictive checks: breadth-ML rho, RSP/ORB group contrast
# =====================================================================
# Both statistics below were NEVER used to fit any model above -- genuine
# out-of-sample checks of the broadcast phenomenon this pipeline cares
# about, not just curve-fit quality.

# Bridging-neuron hypothesis: RSP (a PM area) and ORB (an AT area) are
# used as flagship markers of network engagement. Neurons that reach BOTH
# (the broadcast candidates) are compared to neurons committed to only
# one -- do the "both" neurons also project more widely across the other
# 18 areas, or do they just happen to hit these two extra targets?
def group_contrast(Y):
    rsp = Y[:, rsp_idx] > 0; orb = Y[:, orb_idx] > 0
    both = rsp & orb; single = (rsp & ~orb) | (orb & ~rsp)
    breadth_other = Y[:, other_idx].sum(axis=1)
    if both.sum() < 2 or single.sum() < 2: return np.nan
    return breadth_other[both].mean() - breadth_other[single].mean()

# NB: Spearman rho only detects MONOTONIC ML-breadth relationships. Real
# breadth peaks mid-axis (a hump -- see script 01's "#areas vs ML" column),
# not a monotonic trend, so a modest real rho does not mean "no position
# effect" -- it means the effect isn't monotonic. Models below impose
# smooth monotonic curves and so structurally cannot reproduce that hump;
# expect them to overshoot rho rather than match it (script 05 does the
# direct, shape-aware spatial test). Also note this rho is computed over
# ALL 14 areas (unlike group_contrast's "other 12").
real_rho, _ = spearmanr(t_all, Y_real.sum(axis=1))
real_contrast = group_contrast(Y_real)

N_REP = 300
ppc_rows = []
for name, sim_fn in SIMULATORS.items():
    rhos, contrasts = [], []
    r = np.random.default_rng(123)
    for rep in range(N_REP):
        Y, L = sim_fn(r)
        rho, _ = spearmanr(t_all, Y.sum(axis=1))
        rhos.append(rho)
        contrasts.append(group_contrast(Y))
    rhos = np.array(rhos); contrasts = np.array(contrasts)
    pct_rho = (rhos < real_rho).mean() * 100
    pct_contrast = np.nanmean(contrasts < real_contrast) * 100
    ppc_rows.append((name, rhos, contrasts, pct_rho, pct_contrast))

print(f"\nPosterior-predictive checks ({N_REP} reps), real rho={real_rho:.3f}, real group contrast={real_contrast:.2f}:")
for name, rhos, contrasts, pct_rho, pct_contrast in ppc_rows:
    print(f"  {name:32s} rho pct={pct_rho:5.1f}   contrast pct={pct_contrast:5.1f}")

# =====================================================================
# 6. Summary figure
# =====================================================================
# Panel A: does ML position generalize as a predictor of network identity?
# Panel B: does ML position predict overall breadth on its own (it
#          shouldn't, if the real relationship is a mid-axis hump rather
#          than a monotonic trend -- see the NB above group_contrast)?
# Panel C: does bridging both network anchors (RSP+ORB) predict extra
#          breadth -- the broadcast-phenotype signature this pipeline
#          is chasing?
fig, axes = plt.subplots(1, 3, figsize=(20, 5))

names_sorted = sorted(cv_ll, key=cv_ll.get, reverse=True)
axes[0].barh(names_sorted, [cv_ll[n] for n in names_sorted], color="#4C72B0")
axes[0].set_xlabel("5-fold CV held-out log-likelihood (higher better)")
axes[0].invert_yaxis()

for name, rhos, contrasts, pct_rho, pct_contrast in ppc_rows:
    axes[1].violinplot(rhos, positions=[list(SIMULATORS.keys()).index(name)], widths=0.7, showextrema=False)
axes[1].axhline(real_rho, color="red", lw=1.5, label="real data")
axes[1].set_xticks(range(len(SIMULATORS)))
axes[1].set_xticklabels(list(SIMULATORS.keys()), rotation=90, fontsize=7)
axes[1].set_ylabel("Spearman rho(breadth, ML)"); axes[1].legend(fontsize=8)

for name, rhos, contrasts, pct_rho, pct_contrast in ppc_rows:
    c = contrasts[~np.isnan(contrasts)]
    if len(c) > 1:
        axes[2].violinplot(c, positions=[list(SIMULATORS.keys()).index(name)], widths=0.7, showextrema=False)
axes[2].axhline(real_contrast, color="red", lw=1.5, label="real data")
axes[2].set_xticks(range(len(SIMULATORS)))
axes[2].set_xticklabels(list(SIMULATORS.keys()), rotation=90, fontsize=7)
axes[2].set_ylabel("breadth contrast: both - single group"); axes[2].legend(fontsize=8)

plt.tight_layout()
plt.savefig(OUT_DIR / "sim_models_comparison.png", dpi=140, bbox_inches="tight")
plt.savefig(OUT_DIR / "sim_models_comparison.pdf", bbox_inches="tight")
print("\nSaved sim_models_comparison.png / .pdf")

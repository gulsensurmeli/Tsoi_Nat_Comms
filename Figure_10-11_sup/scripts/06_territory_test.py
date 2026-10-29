# ============================================================================
# 06_territory_test.py
# ROLE   : Directly tests the PM/AT ML-territory asymmetry the refined model
#          (script 04) uses, under this pipeline's alternate PM(7)/AT(7) area
#          definition (see ANALYSIS.md). Mirrors the territory-test logic
#          used by the manuscript's primary PM(7)/AT(13) analysis (a separate
#          pipeline, not included in this repo).
# INPUT  : data/master_detailed_comment.csv
# OUTPUT : outputs/territory_test.{png,pdf}, territory_asymmetric_fit.{png,pdf}
# NOTES  : Fits its own territory parameters live at runtime (no dependency on
#          script 04's hardcoded constants). Fixed-effects logistic (delta_i
#          affects dispersion, not the mean ML-response tested here; cf.
#          script 05). Overlap is measured as the CO-OCCUPANCY zone (both
#          networks' mean prob > 0.5), NOT the difference of fitted logistic
#          midpoints -- AT is a broad medial-to-lateral ramp whose midpoint
#          extrapolates outside the sampled range, so the midpoint difference
#          is not a clean overlap measure.
# ============================================================================

"""
territory_test.py

The refined model (script 04) gives PM and AT their own ML "territories":
  f_PM(t) = A_PM * sigmoid(-k_PM * (t - c_PM))     high medially, soft edge
  f_AT(t) = A_AT * sigmoid( k_AT * (t - c_AT))     rises toward lateral
The central breadth peak requires the two networks to be jointly reachable over
a central ML band. This script asks, directly from the data, (1) whether a
per-network (asymmetric) territory beats a reflection-symmetric one, (2) whether
a central co-occupancy zone (both networks' mean prob > 0.5) robustly exists,
and (3) how much the asymmetry helps out-of-sample.

Reflection-symmetric null: A_PM=A_AT=A, k_PM=k_AT=k, c_PM=-m, c_AT=+m (edges
mirror-imaged about the centre). Nested in the full 6-parameter territory model
via 3 equality constraints -> LRT with 3 df.
"""

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2

from common import OUT_DIR, ispm, A, N, Y, t, a as a0, sig, logit

pm_j = ispm; at_j = ~ispm
rng_master = np.random.default_rng(7)         # seeds only the bootstrap in part (2)
TG = np.linspace(-1, 1, 201)                 # dense ML grid for occupancy curves

# ---- binning (per-network occupancy profile, for the figure) ----
# Same 6-bin scheme used elsewhere in this pipeline (scripts 01/04): bin real
# neurons by ML position, then for each bin take the mean fraction projecting
# to PM areas / AT areas. These are the plotted data points (circles) that the
# fitted curves below are compared against -- not used in any fit.
edges = np.linspace(t.min(), t.max(), 7)
idx   = np.clip(np.digitize(t, edges[1:-1]), 0, 5)
binmask = [idx == b for b in range(6)]
bt   = np.array([t[m].mean() for m in binmask])
occ_pm_dat = np.array([Y[m][:, pm_j].mean() for m in binmask])
occ_at_dat = np.array([Y[m][:, at_j].mean() for m in binmask])

# ---- gradient builders (return (n, A) logit shift) ----
# kbound: reparametrizes the logistic sharpness k as a squashed function of an
# unconstrained "raw" value, so the optimizer below (which works in
# unconstrained space) can never wander into a degenerate near-step (k too
# large) or near-flat (k too small) territory boundary.
def kbound(raw):  return 0.3 + 11.7 / (1.0 + np.exp(-raw))    # k in (0.3, 12)

# Reflection-symmetric territory (the null/comparison model): ONE shared
# amplitude A_ and sharpness k for both networks, with edges forced to mirror
# each other about the ML midpoint (-m for PM, +m for AT). Only 3 free
# territory parameters -- this is what "symmetric" means in this script, not
# to be confused with script 01/04's single-shared-gradient model (that one
# has no separate territory edges at all).
def grad_sym(par, tv):
    A_, k, m = par[0], kbound(par[1]), par[2]
    fpm = A_ * sig(-k * (tv + m)); fat = A_ * sig(k * (tv - m))
    return np.where(pm_j[None, :], fpm[:, None], fat[:, None])

# Per-network territory (the model script 04 actually uses): PM and AT each
# get their own amplitude, sharpness, and edge -- 6 free territory parameters,
# fit directly here (NOT the hardcoded values from script 04's separate fit).
# grad_sym is nested inside this via 3 equality constraints (A1=A2, k1=k2,
# c1=-c2), which is what makes the LRT below valid.
def grad_pn(par, tv):
    A1, k1, c1, A2, k2, c2 = par[0], kbound(par[1]), par[2], par[3], kbound(par[4]), par[5]
    fpm = A1 * sig(-k1 * (tv - c1)); fat = A2 * sig(k2 * (tv - c2))
    G = np.where(pm_j[None, :], fpm[:, None], fat[:, None])
    return G, (c1, c2)

def recal(G, Yset, iters=12):
    """Profile-MLE per-area intercepts: match predicted mean to observed mean.

    For a FIXED territory shift G (whatever candidate parameters the outer
    optimizer in fit() is currently trying), this finds the per-area
    intercepts a_j whose predicted mean matches the area's empirical mean --
    a simpler, no-random-effect cousin of script 04's recalibrated_a() (no
    Gauss-Hermite integration needed here since there's no delta_i to
    marginalize; see the module header NOTES on why delta_i is out of scope
    for this script). Same fixed-point idea: nudge `a` by a fraction (0.8) of
    the gap between target and current predicted mean, iterate.
    """
    emp = logit(np.clip(Yset.mean(0), .02, .98)); a = emp.copy()
    for _ in range(iters):
        pj = sig(a[None, :] + G).mean(0)
        a += (emp - logit(np.clip(pj, .02, .98))) * 0.8
    return a

# Bernoulli log-likelihood of Yset under logit z = a_j + G, written via
# logaddexp for numerical stability (avoids overflow in sig()/log() for
# extreme z, equivalent to sum(Y*log(p) + (1-Y)*log(1-p))).
def loglik(G, Yset, a):
    z = a[None, :] + G
    return (Yset * (-np.logaddexp(0, -z)) + (1 - Yset) * (-np.logaddexp(0, z))).sum()

def fit(kind, tv, Yset, restarts=6, x0=None):
    """Profile-likelihood fit of one territory model ('sym' or 'pn') to (tv, Yset).

    For any candidate territory-parameter vector p, nll(p) builds the
    territory shift G, then calls recal() to re-optimize the per-area
    intercepts for THAT G before scoring the log-likelihood -- i.e. the
    intercepts are profiled out, and Nelder-Mead only has to search the
    territory parameters (3 for 'sym', 6 for 'pn'). Because Nelder-Mead can
    get stuck in a local optimum from a bad start, `restarts` different
    starting points are tried (the first is the caller-supplied x0 or a fixed
    base guess; the rest are that base guess jittered by a fixed cosine
    pattern, not random, so results stay fully reproducible) and the
    best-scoring one is kept.
    """
    build = (lambda p: (grad_sym(p, tv), None)) if kind == 'sym' else (lambda p: grad_pn(p, tv))
    def nll(p):
        try:
            G = build(p)[0]; a = recal(G, Yset); return -loglik(G, Yset, a)
        except Exception:
            return 1e9
    base = ([5.0, 0.0, 0.3] if kind == 'sym'
            else [7.0, 0.0, -0.13, 6.0, 0.0, -0.45])
    starts = [np.array(base if x0 is None else x0, float)]
    for r in range(1, restarts):
        starts.append(np.array(base, float) + 0.3 * np.cos(np.arange(len(base)) * 1.7 + r))
    best = None
    for s in starts:
        res = minimize(nll, s, method='Nelder-Mead',
                       options=dict(maxiter=6000, xatol=1e-6, fatol=1e-9))
        if best is None or res.fun < best.fun:
            best = res
    G = build(best.x)[0]; a = recal(G, Yset)
    return dict(par=best.x, ll=loglik(G, Yset, a), a=a)

def occupancy(par, a, tv=TG):
    """Mean per-network projection probability along the ML grid.

    Always uses the per-network gradient builder (grad_pn) regardless of
    which model's parameters are passed in -- callers pass either the
    per-network fit's own parameters, or (for the symmetric model, see the
    figure section) reconstruct the equivalent per-network-shaped call
    themselves. This is the curve actually plotted in panel A.
    """
    Gp, edges_ = grad_pn(par, tv); P = sig(a[None, :] + Gp)
    return P[:, pm_j].mean(1), P[:, at_j].mean(1), edges_

def cooccupancy_width(par, a, thresh=0.5, tv=TG):
    """Direct, extrapolation-safe measure of central territory overlap.

    Rather than comparing the two fitted logistic midpoints (c_PM vs c_AT --
    flagged in the module NOTES as unreliable, since AT's midpoint can
    extrapolate outside the sampled ML range), this scans the dense ML grid
    and finds where BOTH networks' mean projection probability exceeds
    `thresh` (0.5) at the same time. `width` converts the fraction of the
    grid in that zone into ML units; `peak` is the highest value the SMALLER
    of the two probabilities reaches in that zone (how strongly "joint" the
    overlap gets, not just how wide it is); `both` is the boolean mask over
    `tv` used later to shade the zone in the figure.
    """
    opm, oat, _ = occupancy(par, a, tv)
    both = (opm > thresh) & (oat > thresh)
    width = both.mean() * (tv[-1] - tv[0])           # in ML (t) units
    peak = np.minimum(opm, oat).max()                # max joint reachability
    return width, peak, both

# ======================================================================
# (1) ASYMMETRY LRT on the full data
# ======================================================================
# Fit both models once on all 160 neurons. Since 'sym' (3 territory params)
# is nested inside 'pn' (6 territory params) via 3 equality constraints, twice
# the log-likelihood gap is chi2-distributed with 3 df under the null that the
# extra per-network flexibility doesn't help -- a standard nested-model LRT.
# k_sym/k_pn count ALL free parameters (territory params + one intercept per
# area, A=14) for AIC/BIC, which penalize model size more conservatively than
# the LRT's p-value alone (relevant given nobs = N*A is large).
fs = fit('sym', t, Y); fp = fit('pn', t, Y)
lrt = 2 * (fp['ll'] - fs['ll']); df = 3; p_lrt = chi2.sf(lrt, df)
nobs = N * A
k_sym, k_pn = 3 + A, 6 + A
def aic(ll, k): return 2 * k - 2 * ll
def bic(ll, k): return k * np.log(nobs) - 2 * ll
cpm_hat, cat_hat = fp['par'][2], fp['par'][5]
# both_hat: boolean mask over TG marking the full-data co-occupancy zone --
# reused unchanged in the figure section below to shade panel A.
w_hat, peak_hat, both_hat = cooccupancy_width(fp['par'], fp['a'])

print("=" * 72)
print("(1) ASYMMETRY: per-network territory vs reflection-symmetric")
print("=" * 72)
print(f"  logLik  symmetric = {fs['ll']:.2f}  (k={k_sym})")
print(f"  logLik  per-net   = {fp['ll']:.2f}  (k={k_pn})")
print(f"  LRT chi2({df}) = {lrt:.2f}   p = {p_lrt:.2e}")
print(f"  AIC  sym={aic(fs['ll'],k_sym):.1f}  pn={aic(fp['ll'],k_pn):.1f}   "
      f"(dAIC={aic(fs['ll'],k_sym)-aic(fp['ll'],k_pn):+.1f}, +ve favours per-net)")
print(f"  BIC  sym={bic(fs['ll'],k_sym):.1f}  pn={bic(fp['ll'],k_pn):.1f}   "
      f"(dBIC={bic(fs['ll'],k_sym)-bic(fp['ll'],k_pn):+.1f})")
print(f"  fitted logistic midpoints (note c_AT extrapolates < -1; use co-occupancy below):")
print(f"    c_PM={cpm_hat:+.3f}  c_AT={cat_hat:+.3f}")
print(f"  full-data co-occupancy width (both mean prob>0.5) = {w_hat:.3f} ML units; "
      f"peak joint prob = {peak_hat:.2f}")

# ======================================================================
# (2) BOOTSTRAP the central co-occupancy zone (neuron resampling)
# ======================================================================
# Is the co-occupancy zone found above a robust feature of the data, or just
# an artifact of this one point estimate? Resample neurons with replacement
# (rows of t/Y, so each resample is still a valid 160-neuron dataset), refit
# the per-network model on each resample, and recompute its co-occupancy
# width/peak. restarts=1 (vs. 6 for the full-data fits above) and warm-starting
# from the full-data optimum (x0=fp['par']) trade a little per-resample fit
# quality for 600x fewer optimizer restarts -- reasonable since each resample
# is a small perturbation of the data the full fit already solved well.
NBOOT = 600
w_bs = np.empty(NBOOT); peak_bs = np.empty(NBOOT)
for b in range(NBOOT):
    ridx = rng_master.integers(0, N, N)
    fb = fit('pn', t[ridx], Y[ridx], restarts=1, x0=fp['par'])
    w_bs[b], peak_bs[b], _ = cooccupancy_width(fb['par'], fb['a'])
def ci(x): return np.percentile(x, [2.5, 50, 97.5])
w_ci, peak_ci = ci(w_bs), ci(peak_bs)
# Fraction of resamples where a co-occupancy zone exists at all (width > 0) --
# reported alongside the CI as a direct "how often is there any central
# overlap zone" answer, since a CI alone doesn't rule out some resamples
# having none.
frac_pos = float((w_bs > 0).mean())

print("\n" + "=" * 72)
print(f"(2) BOOTSTRAP central co-occupancy ({NBOOT} neuron-resamples)")
print("=" * 72)
print(f"  co-occupancy width  median {w_ci[1]:.3f}   95% CI [{w_ci[0]:.3f}, {w_ci[2]:.3f}] ML units")
print(f"  peak joint prob     median {peak_ci[1]:.2f}   95% CI [{peak_ci[0]:.2f}, {peak_ci[2]:.2f}]")
print(f"  P(a central both->0.5 zone exists) = {frac_pos:.3f}")

# ======================================================================
# (3) 5-fold held-out gradient log-likelihood (symmetric vs per-network)
# ======================================================================
# Does the per-network model's extra flexibility generalize, or does it just
# fit the training data closer? Split neurons into 5 folds; for each fold,
# refit BOTH models (warm-started from their own full-data optimum) on the
# other 4/5, then score per-neuron log-likelihood on the held-out fold using
# the training fold's fitted territory params and intercepts -- genuinely
# out-of-sample, since neither the territory shape nor the intercepts saw
# the held-out neurons during fitting. d_ho is the paired per-neuron
# difference (pn minus sym), which is what lets the SE below be a paired SE
# rather than a two-sample one.
perm = np.random.default_rng(0).permutation(N)
folds = np.array_split(perm, 5)
ho = {'sym': np.zeros(N), 'pn': np.zeros(N)}
for fold in folds:
    ev = fold; tr = np.setdiff1d(perm, fold)
    for kind in ('sym', 'pn'):
        f = fit(kind, t[tr], Y[tr], restarts=3, x0=(fs['par'] if kind == 'sym' else fp['par']))
        G_ev = (grad_sym(f['par'], t[ev]) if kind == 'sym' else grad_pn(f['par'], t[ev])[0])
        z = f['a'][None, :] + G_ev
        ho[kind][ev] = -(Y[ev] * np.logaddexp(0, -z) + (1 - Y[ev]) * np.logaddexp(0, z)).sum(1)
d_ho = ho['pn'] - ho['sym']
print("\n" + "=" * 72)
print("(3) 5-FOLD held-out per-neuron log-likelihood")
print("=" * 72)
print(f"  symmetric   = {ho['sym'].mean():.3f}")
print(f"  per-network = {ho['pn'].mean():.3f}")
print(f"  paired dLL/neuron (pn - sym) = {d_ho.mean():+.3f}  SE {d_ho.std(ddof=1)/np.sqrt(N):.3f}")
print("  (small by design -- the asymmetry is justified by the breadth-vs-ML")
print("   profile/both-group it reproduces, not per-neuron prediction; see script 04.)")

# ======================================================================
# FIGURE
# ======================================================================
plt.rcParams.update({"font.family": "sans-serif", "font.size": 8, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5})
PMc, ATc, OV, GREY = "#FF7800", "#3485DC", "#6E5A9E", "#9a9a9a"

# --- Panel A content, factored out so it can be drawn both inside the
# composite 3-panel figure (with the symmetric-model dashed curves for
# comparison) and as its own standalone figure (data + asymmetric fit only,
# no symmetric curves) without duplicating the plotting logic. ---
opm_pn, oat_pn, _ = occupancy(fp['par'], fp['a'])
opm_sym = sig(fs['a'][None, :] + grad_sym(fs['par'], TG))[:, pm_j].mean(1)
oat_sym = sig(fs['a'][None, :] + grad_sym(fs['par'], TG))[:, at_j].mean(1)
zone = TG[both_hat]

def draw_panelA(ax, show_symmetric, title, fontsize_title=7.3):
    if zone.size:
        ax.axvspan(zone.min(), zone.max(), color=OV, alpha=0.14, lw=0)
    ax.axhline(0.5, color=GREY, lw=0.5, ls=":")
    ax.plot(TG, opm_pn, color=PMc, lw=1.8, label="PM"); ax.plot(TG, oat_pn, color=ATc, lw=1.8, label="AT")
    if show_symmetric:
        ax.plot(TG, opm_sym, color=PMc, lw=1.0, ls=(0, (4, 2)), alpha=0.6)
        ax.plot(TG, oat_sym, color=ATc, lw=1.0, ls=(0, (4, 2)), alpha=0.6)
    ax.plot(bt, occ_pm_dat, 'o', color=PMc, ms=3.5); ax.plot(bt, occ_at_dat, 'o', color=ATc, ms=3.5)
    ax.set_xlabel("ML position (medial → lateral)"); ax.set_ylabel("Projection probability")
    ax.set_title(title, fontsize=fontsize_title)
    ax.legend(fontsize=7, frameon=False, loc="center right"); ax.spines[["top", "right"]].set_visible(False)

fig, ax = plt.subplots(1, 3, figsize=(8.6, 2.7))

# --- Panel A: fitted territories + data + co-occupancy zone ---
draw_panelA(ax[0], show_symmetric=True,
            title="Per-network ML territories\n(solid=per-net, dashed=symmetric, band=both>0.5)")

# --- Panel B: bootstrap distribution of co-occupancy width (part 2's w_bs) ---
ax[1].hist(w_bs, bins=30, color=OV, alpha=0.8, edgecolor="white", lw=0.3)
ax[1].axvline(0, color="k", lw=1.0)
ax[1].axvline(w_ci[0], color=GREY, lw=0.8, ls="--"); ax[1].axvline(w_ci[2], color=GREY, lw=0.8, ls="--")
ax[1].set_xlabel("co-occupancy width (ML units, both>0.5)"); ax[1].set_ylabel("bootstrap count")
ax[1].set_title(f"Central both-zone robust: P>0={frac_pos:.2f}\n95% CI [{w_ci[0]:.2f}, {w_ci[2]:.2f}]", fontsize=7.3)
ax[1].spines[["top", "right"]].set_visible(False)

# --- Panel C: held-out LL, symmetric vs per-network (part 3's ho/d_ho) ---
means = [ho['sym'].mean(), ho['pn'].mean()]
ses = [ho['sym'].std(ddof=1)/np.sqrt(N), ho['pn'].std(ddof=1)/np.sqrt(N)]
ax[2].bar([0, 1], means, yerr=ses, color=[GREY, OV], width=0.6, capsize=3)
ax[2].set_xticks([0, 1]); ax[2].set_xticklabels(["symmetric", "per-network"], fontsize=7)
ax[2].set_ylabel("held-out logLik / neuron")
ax[2].set_title(f"Out-of-sample: Δ={d_ho.mean():+.2f}/neuron\n(profile-driven, not prediction)", fontsize=7.3)
ax[2].set_ylim(min(means) - 0.6, max(means) + 0.3)
ax[2].spines[["top", "right"]].set_visible(False)

fig.suptitle("Are the PM and AT ML-response functions genuinely different?", fontsize=9, y=1.04)
fig.tight_layout(w_pad=1.4)
fig.savefig(OUT_DIR / "territory_test.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT_DIR / "territory_test.pdf", bbox_inches="tight")
print("\nsaved territory_test.png / .pdf")

# --- Standalone figure: data + asymmetric (per-network) fit only, no
# symmetric-model comparison curves -- same panel A content, on its own. ---
figA, axA = plt.subplots(figsize=(4.3, 3.4))
draw_panelA(axA, show_symmetric=False,
            title="Per-network ML territories\n(data vs. asymmetric fit; band = both>0.5)",
            fontsize_title=8.5)
figA.tight_layout()
figA.savefig(OUT_DIR / "territory_asymmetric_fit.svg", dpi=300, bbox_inches="tight")
figA.savefig(OUT_DIR / "territory_asymmetric_fit.pdf", bbox_inches="tight")
print("saved territory_asymmetric_fit.svg / .pdf")

# ============================================================================
# fit_refined_params.py   (derivation utility -- NOT a numbered pipeline step)
# ROLE   : Re-derive, end-to-end from the data, the per-network territory +
#          continuous-propensity parameters that are HARDCODED in
#          04_refined_model_fit.py, and check they still match. Run this
#          after any change to the refined model's mechanism, or to the
#          area definition in common.py, and paste the printed block back
#          into script 04.
# INPUT  : data/master_detailed_comment.csv  (via common.py)
# OUTPUT : stdout only (no figure). Prints a copy-paste-ready parameter block
#          and a PASS/DIFFERS check against script 04's current constants
#          (skipped on the very first run -- see CURRENT below).
# NOTES  : Fits to the FIGURE's summary statistics (breadth-vs-ML profile, RSP/
#          ORB group target numbers, breadth mean & variance) -- NOT the full
#          per-neuron likelihood (that objective is used for model *selection*
#          in scripts 02/06, which is a different question; see ANALYSIS.md).
#          delta_i is integrated out with Gauss-Hermite quadrature; per-area
#          intercepts are profiled out by recalibration; the 7 structural
#          parameters are fit by Nelder-Mead. Deterministic (no RNG) --
#          reproducible to optimizer tolerance. Runtime several minutes
#          (6 sharpness caps x 14 restarts); run with `python -u` for live
#          progress (stdout is block-buffered).
# ============================================================================
"""
fit_refined_params.py

Model (identical to 04_refined_model_fit.py):
    e_ij = sigmoid( a_j + f_net(j)(t_i) + delta_i )
      f_PM(t) = A_PM * sigmoid(-k_PM * (t - c_PM))
      f_AT(t) = A_AT * sigmoid( k_AT * (t - c_AT))
      delta_i ~ Normal(0, sigma^2)

Fit target: the summary statistics the refined-model figure reproduces --
  - breadth (# target areas) vs ML position, 6 ML bins;
  - number of other-12 targets for the RSP-only / ORB-only / both groups;
  - breadth mean and variance across neurons.

Boundary sharpness k is capped: an unconstrained fit can drive k_PM to a
near-step that overfits the coarse 6-bin profile. We sweep an upper cap KMAX
and take the smallest cap whose fitted profile is still within tolerance.
See ANALYSIS.md ("hardcoded fitted parameters") for why these live as pasted
constants in script 04 rather than being fit at import time.
"""

import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import minimize

from common import (
    ispm, A, N, Y as Yr, t, a as a0, sig, logit,
    rsp_idx, orb_idx, other_idx,
)

pm_j = ispm; at_j = ~ispm
R_prev = Yr.mean(0)                         # per-area marginal prevalence (recal target)

# 16-node Gauss-Hermite quadrature to average over delta_i ~ N(0, sigma^2)
_GHx, _GHw = hermgauss(16); _LOGW = np.log(_GHw / np.sqrt(np.pi))

# ---- fit targets: real summary statistics ----
edges = np.linspace(t.min(), t.max(), 7)
idx   = np.clip(np.digitize(t, edges[1:-1]), 0, 5)
binmask = [idx == b for b in range(6)]
bt   = np.array([t[m].mean() for m in binmask])
rs = Yr[:, rsp_idx] > 0; ob = Yr[:, orb_idx] > 0
real_group = np.array([Yr[rs & ~ob][:, other_idx].sum(1).mean(),
                       Yr[~rs & ob][:, other_idx].sum(1).mean(),
                       Yr[rs & ob][:, other_idx].sum(1).mean()])
real_bml   = np.array([Yr[m].sum(1).mean() for m in binmask])
real_bmean = Yr.sum(1).mean()
real_bvar  = Yr.sum(1).var(ddof=1)

# ---- model pieces ----
def grad_pn(A_pm, k_pm, c_pm, A_at, k_at, c_at):
    """Per-network ML logit shift, shape (N, A)."""
    fpm = A_pm * sig(-k_pm * (t - c_pm)); fat = A_at * sig(k_at * (t - c_at))
    return np.where(pm_j[None, :], fpm[:, None], fat[:, None])

def comps_delta(sigma):
    """Quadrature mixture components (log-weight, delta shift) for N(0, sigma^2)."""
    return [(_LOGW[q], np.sqrt(2) * sigma * _GHx[q]) for q in range(len(_GHx))]

def recalibrate(G, comps, iters=12):
    """Profile-MLE per-area intercepts: match model marginal prevalence to data."""
    emp = logit(np.clip(R_prev, .02, .98)); a = a0.copy()
    for _ in range(iters):
        pj = np.zeros(A)
        for lw, sh in comps:
            pj += (np.exp(lw) * sig(a[None, :] + G + sh)).mean(0)
        a += (emp - logit(np.clip(pj, .02, .98))) * 0.7
    return a

def stats(G, comps):
    """Analytic summary statistics of the model (delta integrated out)."""
    a = recalibrate(G, comps)
    m = np.zeros(N); m2 = np.zeros(N); vb = np.zeros(N)
    numg = np.zeros(3); deng = np.zeros(3)
    for lw, sh in comps:
        w = np.exp(lw); p = sig(a[None, :] + G + sh); M = p.sum(1)
        m += w * M; m2 += w * M**2; vb += w * (p * (1 - p)).sum(1)
        pr, po = p[:, rsp_idx], p[:, orb_idx]; oth = p[:, other_idx].sum(1)
        for gi, want in enumerate([(1, 0), (0, 1), (1, 1)]):
            mb = (pr if want[0] else 1 - pr) * (po if want[1] else 1 - po)
            numg[gi] += (w * mb * oth).sum(); deng[gi] += (w * mb).sum()
    var = m.var(ddof=0) + (vb + m2 - m**2).mean()
    return dict(bml=np.array([m[k].mean() for k in binmask]),
                groups=numg / deng, bmean=m.mean(), bvar=var)

def loss(s):
    """Weighted summary-statistic objective (primary = breadth-vs-ML profile)."""
    return (np.mean((s['bml'] - real_bml)**2)          # profile (6 bins)
            + np.mean((s['groups'] - real_group)**2)    # RSP/ORB group targets
            + 0.5 * (s['bmean'] - real_bmean)**2         # breadth mean
            + 0.06 * (s['bvar'] - real_bvar)**2)         # breadth variance

# ---- fit at a given boundary-sharpness cap KMAX ----
def fit_capped(KMAX, restarts=14):
    # unconstrained params: [A_pm, k_pm_raw, c_pm, A_at, k_at_raw, c_at, log_sigma]
    # k = 0.3 + (KMAX-0.3)*sigmoid(raw) keeps 0.3 < k < KMAX
    def unpack(x):
        A_pm = x[0]; k_pm = 0.3 + (KMAX - 0.3) * sig(x[1]); c_pm = x[2]
        A_at = x[3]; k_at = 0.3 + (KMAX - 0.3) * sig(x[4]); c_at = x[5]
        sigma = np.exp(x[6])
        return A_pm, k_pm, c_pm, A_at, k_at, c_at, sigma
    def obj(x):
        try:
            A_pm, k_pm, c_pm, A_at, k_at, c_at, sigma = unpack(x)
            G = grad_pn(A_pm, k_pm, c_pm, A_at, k_at, c_at)
            return loss(stats(G, comps_delta(sigma))) + 0.002 * (A_pm**2 + A_at**2)  # ridge on amps
        except Exception:
            return 1e9
    base = np.array([5.0, 1.0, 0.1, 5.0, 0.5, -0.4, np.log(1.4)])
    best = None
    for r in range(restarts):
        x0 = base + (0 if r == 0 else 0.3 * np.cos(np.arange(7) * 1.7 + r))   # deterministic jitter
        res = minimize(obj, x0, method='Nelder-Mead',
                       options=dict(maxiter=14000, xatol=1e-7, fatol=1e-10))
        if best is None or res.fun < best.fun:
            best = res
    par = unpack(best.x)
    return par, stats(grad_pn(*par[:6]), comps_delta(par[6]))

# ---- sweep the cap; pick the smallest that keeps the profile within tolerance ----
TOL = 0.55       # max allowed |per-bin breadth residual|
print("Sweeping boundary-sharpness cap KMAX (smallest within tolerance is chosen):")
print(f"{'KMAX':>5} {'k_pm':>6} {'k_at':>6} {'c_pm':>7} {'c_at':>7} {'sigma':>6} "
      f"{'max|resid|':>10} {'both':>6}")
sweep = {}
for KMAX in [4, 6, 8, 10, 15, 60]:
    par, s = fit_capped(KMAX); sweep[KMAX] = (par, s)
    A_pm, k_pm, c_pm, A_at, k_at, c_at, sigma = par
    mr = np.abs(s['bml'] - real_bml).max()
    print(f"{KMAX:>5} {k_pm:6.2f} {k_at:6.2f} {c_pm:+7.3f} {c_at:+7.3f} {sigma:6.3f} "
          f"{mr:10.2f} {s['groups'][2]:6.2f}", flush=True)

chosen = next((K for K in [4, 6, 8, 10, 15, 60]
               if np.abs(sweep[K][1]['bml'] - real_bml).max() <= TOL), 60)
par, s = sweep[chosen]
A_pm, k_pm, c_pm, A_at, k_at, c_at, sigma = par

print(f"\nChosen KMAX = {chosen} (smallest with max|profile residual| <= {TOL})")
print("Fit diagnostics vs data:")
print(f"  breadth profile : " + " ".join(f"{v:5.2f}" for v in s['bml']))
print(f"  real profile    : " + " ".join(f"{v:5.2f}" for v in real_bml))
print(f"  residual        : " + " ".join(f"{v:+5.2f}" for v in s['bml'] - real_bml))
print(f"  groups (RSP/ORB/both) model {np.round(s['groups'],2).tolist()}  "
      f"real {np.round(real_group,2).tolist()}")
print(f"  breadth mean {s['bmean']:.2f} (real {real_bmean:.2f})  "
      f"var {s['bvar']:.2f} (real {real_bvar:.2f})")

print("\n" + "=" * 70)
print("Copy-paste block for 04_refined_model_fit.py:")
print("=" * 70)
print(f"A_PM, K_PM, C_PM = {A_pm:.3f}, {k_pm:.3f}, {c_pm:.3f}")
print(f"A_AT, K_AT, C_AT = {A_at:.3f}, {k_at:.3f}, {c_at:.3f}")
print(f"SIGMA            = {sigma:.3f}")

# ---- check against the constants currently pasted in script 04 ----
# (kept here as a literal so this utility needs no import of script 04, which
#  would execute that module's figure code. Update if you re-paste new values.)
# Baseline established by the first fit of this area definition, now pasted
# into script 04 (04_refined_model_fit.py) -- kept here as a literal so this
# utility needs no import of script 04, which would execute that module's
# figure code. Update if you re-paste new values.
CURRENT = dict(A_PM=5.869, K_PM=4.000, C_PM=-0.009,
               A_AT=8.316, K_AT=4.000, C_AT=-0.421, SIGMA=1.513)
got = dict(A_PM=A_pm, K_PM=k_pm, C_PM=c_pm, A_AT=A_at, K_AT=k_at, C_AT=c_at, SIGMA=sigma)
print("\n" + "=" * 70)
if CURRENT is None:
    print("No prior baseline yet -- this run establishes it.")
    print("Paste the block above into script 04, then set CURRENT here to the same")
    print("values so future re-fits can PASS/DIFFERS-check against it.")
else:
    worst = max(abs(got[k] - CURRENT[k]) for k in CURRENT)
    if worst <= 0.01:
        print(f"PASS: re-derived params match script 04 (max abs diff {worst:.4f} <= 0.01).")
    else:
        print(f"DIFFERS from script 04 (max abs diff {worst:.4f} > 0.01). If the model or")
        print("objective changed intentionally, paste the block above into script 04.")
        for k in CURRENT:
            print(f"    {k:5s} script04={CURRENT[k]:+.3f}  re-derived={got[k]:+.3f}")

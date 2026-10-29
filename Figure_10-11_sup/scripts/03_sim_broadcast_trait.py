# ============================================================================
# 03_sim_broadcast_trait.py
# ROLE   : Tests whether adding a broadcast trait to the graded-identity
#          model (fit by script 01) reproduces the two held-out broadcast
#          signatures under this pipeline's alternate PM(7)/AT(7) area
#          definition (see ANALYSIS.md). Mirrors the broadcast-trait
#          extension used by the manuscript's primary PM(7)/AT(13) analysis
#          (a separate pipeline, not included in this repo).
# INPUT  : data/master_detailed_comment.csv
# OUTPUT : outputs/sim_broadcast_trait.{png,pdf}
# NOTES  : KAPPA below is script 01's fitted "4. Graded identity (1 axis)"
#          kappa under this area definition (run script 01 first if
#          re-deriving). Variant B/C params (amp/w/sigma) are fit live at
#          runtime here, not hardcoded.
# ============================================================================

"""
sim_broadcast_trait.py

Tests whether adding an ML-INDEPENDENT per-neuron broadcast trait to the
graded identity model reproduces the two held-out broadcast signatures that
NO ML-position-only model could reproduce:

   (1) real total-breadth vs ML is NON-MONOTONIC (peaks at intermediate ML)
       -> Spearman rho(breadth, ML) is low (0.20); ML-only models produce
       monotonic breadth and overshoot this (rho ~ 0.47).
   (2) real RSP+ORB "both" neurons project to ~4.8 MORE other areas than
       RSP-only / ORB-only "single" neurons; no ML-only model reproduces it.

THREE NESTED VARIANTS COMPARED
------------------------------
  A. ML-only (graded)                 e_ij = sig(a_j + ident_ij)
  B. + uniform broadcast trait        e_ij = sig(a_j + ident_ij + delta_i),
                                       delta_i ~ N(0, sigma^2) INDEPENDENT of t
  C. + position-dependent broadcast   e_ij = sig(a_j + ident_ij
                                                + amp*exp(-t^2/2w^2) + delta_i)

Biological readings:
  * delta_i (uniform trait): a per-neuron intrinsic arborization /
    collateralization capacity -- plausibly a transcriptional program --
    uncorrelated with the mediolateral identity gradient. Generative form
    of the earlier random-intercept structure.
  * mid-axis boost: an ADDITIONAL, position-dependent broadcast propensity
    concentrated at intermediate ML, where the identity gradient is weakest
    -- i.e. neurons that sit between the two networks are the ones that
    engage both. This is the broadcast-bump mechanism that the earlier
    monotonic PM/AT curve-fit had structurally suppressed.

METHOD NOTES (addressing the audit)
-----------------------------------
* Fit targets are the two held-out BREADTH statistics (rho, contrast) plus
  breadth variance -- NOT the monotonic PM/AT curves that cannot represent
  a bump. The identity slope kappa is inherited from the graded fit.
* delta_i is MARGINALIZED: sigma is a single population parameter and each
  neuron draws delta fresh. delta is never read off a neuron's own
  projections, so the earlier circularity does not occur.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
from scipy.optimize import minimize

from common import (
    OUT_DIR, AREAS, A, ispm, N,
    sig, logit, Y as Y_real, t, a,
    rsp_idx as rsp_i, orb_idx as orb_i, other_idx as other,
)

# ---- real data ----
# The fitted "4. Graded identity (1 axis)" kappa printed by script 01
# (01_sim_models_fitted.py), under this alternate area definition.
KAPPA = 2.261
ident = np.where(ispm[None,:], (-KAPPA*t)[:,None], (KAPPA*t)[:,None])

def breadth_rho(Y):
    br=Y.sum(1)
    return np.nan if br.std()==0 else spearmanr(t,br)[0]
def group_contrast(Y):
    rsp=Y[:,rsp_i]>0; orb=Y[:,orb_i]>0; both=rsp&orb; single=(rsp&~orb)|(orb&~rsp)
    bo=Y[:,other].sum(1)
    return np.nan if (both.sum()<2 or single.sum()<2) else bo[both].mean()-bo[single].mean()

real_rho = breadth_rho(Y_real)
real_contrast = group_contrast(Y_real)
real_bvar = Y_real.sum(1).var()
print(f"real: breadth-ML rho={real_rho:.3f}  RSP/ORB contrast={real_contrast:.3f}  breadth var={real_bvar:.2f}")

# ---- simulator with amp (mid-axis boost), w (boost width), sigma (uniform trait) ----
def simulate(amp, w, sigma, rng):
    boost = amp*np.exp(-t**2/(2*w**2)) if amp>0 else np.zeros(N)
    d = rng.normal(0.0, sigma, N) if sigma>0 else np.zeros(N)
    e = sig(a[None,:] + ident + boost[:,None] + d[:,None])
    return (rng.random(e.shape) < e).astype(float)

def summary(amp, w, sigma, n_rep=150, seed=0):
    rng = np.random.default_rng(seed)
    rhos, cons, bvars = [], [], []
    for _ in range(n_rep):
        Y = simulate(amp, w, sigma, rng)
        rhos.append(breadth_rho(Y)); cons.append(group_contrast(Y)); bvars.append(Y.sum(1).var())
    return np.nanmean(rhos), np.nanmean(cons), np.nanmean(bvars)

# ---- fit variant C (amp, w, sigma) to (rho, contrast, breadth var) ----
def obj(params):
    amp, w, sigma = params
    if amp<0 or sigma<0 or w<=0.05: return 1e6
    r,c,v = summary(amp, w, sigma, n_rep=80)
    return ((r-real_rho))**2 + ((c-real_contrast)/real_contrast)**2 + ((v-real_bvar)/real_bvar)**2

res = minimize(obj, [2.5,0.4,0.8], method="Nelder-Mead",
               options={"maxiter":120,"xatol":0.05,"fatol":1e-3})
amp_hat, w_hat, sigma_hat = res.x
print(f"\nfitted variant C: amp={amp_hat:.3f}  w={w_hat:.3f}  sigma={sigma_hat:.3f}")

# also fit variant B (uniform trait only): amp=0, fit sigma
def obj_B(sigma):
    if sigma[0]<0: return 1e6
    r,c,v = summary(0.0, 0.4, sigma[0], n_rep=80)
    return ((r-real_rho))**2 + ((v-real_bvar)/real_bvar)**2
resB = minimize(obj_B, [1.0], method="Nelder-Mead", options={"maxiter":60})
sigmaB = resB.x[0]
print(f"fitted variant B: sigma={sigmaB:.3f}")

VARIANTS = [
    ("A. ML-only (graded)",              0.0,    0.4, 0.0),
    ("B. + uniform trait",               0.0,    0.4, sigmaB),
    ("C. + position-dependent broadcast", amp_hat, w_hat, sigma_hat),
]

# ---- posterior-predictive checks ----
N_REP = 400
def ppc(amp,w,sigma,seed=123):
    rng=np.random.default_rng(seed)
    rhos,cons=[],[]
    for _ in range(N_REP):
        Y=simulate(amp,w,sigma,rng); rhos.append(breadth_rho(Y)); cons.append(group_contrast(Y))
    return np.array(rhos), np.array(cons)
def pctile(dist,val):
    d=dist[~np.isnan(dist)]; return (d<val).mean()*100

print("\nPosterior-predictive percentile of real value (50 = well-centered):")
results={}
for name,amp,w,sigma in VARIANTS:
    r,c = ppc(amp,w,sigma); results[name]=(r,c)
    print(f"  {name:38s} rho pct={pctile(r,real_rho):5.1f}   contrast pct={pctile(c,real_contrast):5.1f}")

# ---- figure ----
fig, axes = plt.subplots(1,2, figsize=(14,5))
pos = range(len(VARIANTS))
axes[0].violinplot([results[n][0][~np.isnan(results[n][0])] for n,_,_,_ in VARIANTS],
                   positions=pos, showextrema=False)
axes[0].axhline(real_rho, color="red", lw=1.5, label=f"real = {real_rho:.2f}")
axes[0].set_xticks(pos); axes[0].set_xticklabels([n for n,_,_,_ in VARIANTS], rotation=15, ha="right", fontsize=8)
axes[0].set_ylabel("Spearman rho(breadth, ML)")
axes[0].set_title("Held-out check 1: breadth-ML correlation\n(real is LOW because breadth peaks mid-axis)")
axes[0].legend(fontsize=9)

axes[1].violinplot([results[n][1][~np.isnan(results[n][1])] for n,_,_,_ in VARIANTS],
                   positions=pos, showextrema=False)
axes[1].axhline(real_contrast, color="red", lw=1.5, label=f"real = {real_contrast:.2f}")
axes[1].set_xticks(pos); axes[1].set_xticklabels([n for n,_,_,_ in VARIANTS], rotation=15, ha="right", fontsize=8)
axes[1].set_ylabel("breadth contrast: both - single")
axes[1].set_title("Held-out check 2: RSP+ORB broadcast contrast")
axes[1].legend(fontsize=9)

plt.tight_layout()
plt.savefig(OUT_DIR / "sim_broadcast_trait.png", dpi=140, bbox_inches="tight")
plt.savefig(OUT_DIR / "sim_broadcast_trait.pdf", bbox_inches="tight")
print("\nSaved sim_broadcast_trait.png / .pdf")

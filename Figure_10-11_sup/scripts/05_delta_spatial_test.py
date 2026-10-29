# ============================================================================
# 05_delta_spatial_test.py
# ROLE   : Tests whether the refined model's per-neuron breadth propensity
#          (delta_i, from script 04) carries genuine spatial structure along
#          the ML axis, under this pipeline's alternate PM(7)/AT(7) area
#          definition (see ANALYSIS.md). Mirrors the spatial-test logic used
#          by the manuscript's primary PM(7)/AT(13) analysis (a separate
#          pipeline, not included in this repo).
# INPUT  : data/master_detailed_comment.csv
# OUTPUT : outputs/delta_spatial_test.{png,pdf}
# NOTES  : Uses statsmodels BinomialBayesMixedGLM. No hardcoded fit-specific
#          constants in this script -- fully self-contained per run.
# ============================================================================

"""
delta_spatial_test.py

Per-neuron broadcast-propensity model with the guardrails discussed:

  logit P(project_ij) = a_j
                        + [PM polynomial in t]*I(j in PM)
                        + [AT polynomial in t]*I(j in AT)      (directional ML gradient)
                        + delta_i                              (per-neuron breadth propensity)

delta_i is a partially-pooled (shrunk) random intercept, estimated by
variational Bayes (BinomialBayesMixedGLM). We then ask whether the fitted
delta_i has any SPATIAL distribution along the ML axis:

  * monotonic trend      : Spearman(delta_i, t_i)
  * central concentration: correlation(delta_i, -t_i^2)  (intermediate-high)

GUARDRAIL: the only way a mis-specified ML gradient could leak spatial
structure into delta is through the *non-directional* (breadth) component of
its residual. We therefore re-fit with the directional gradient at linear,
quadratic and cubic flexibility and require any delta-vs-t relationship to be
STABLE across them. A relationship that appears only at low flexibility, or
moves as the gradient is enriched, is gradient misfit, not spatial broadcasting.
delta_i never sees t_i during fitting (t enters only the shared gradient), so
the test is not self-predictive; shrinkage biases the test toward the null,
so a surviving relationship is conservative.
"""

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM
from scipy.stats import spearmanr, pearsonr
import warnings; warnings.filterwarnings("ignore")

from common import (
    OUT_DIR, sub, AREAS, ispm, A, N,
    Y, t,
)

breadth = Y.sum(1)

rows=[]
for i in range(N):
    for j in range(A):
        rows.append(dict(y=Y[i,j], neuron=f"n{i:03d}", area=AREAS[j],
                         pm=int(ispm[j]), at=int(not ispm[j]), t=t[i]))
D=pd.DataFrame(rows)
for k in [1,2,3]:
    D[f'pm_t{k}']=D['pm']*D['t']**k; D[f'at_t{k}']=D['at']*D['t']**k

def fit_delta(Kpoly):
    grad=[f'pm_t{k}' for k in range(1,Kpoly+1)]+[f'at_t{k}' for k in range(1,Kpoly+1)]
    fe="y ~ C(area) + " + " + ".join(grad)
    res=BinomialBayesMixedGLM.from_formula(fe, {"neuron":"0 + C(neuron)"}, D).fit_vb()
    delta=np.asarray(res.vc_mean).ravel()            # (N,) posterior-mean random intercepts
    sd_delta=np.exp(res.vcp_mean[0])
    return delta, sd_delta

results={}
print(f"{'grad':>6} | {'SD(delta)':>9} | {'rho(delta,t)':>13} p | {'corr(delta,-t^2)':>16} p | {'corr(delta,breadth)':>18}")
for K,label in [(1,"linear"),(2,"quad"),(3,"cubic")]:
    delta,sdd=fit_delta(K)
    rho_t,p_t=spearmanr(t,delta)
    r_c,p_c=pearsonr(-t**2,delta)            # central concentration: high delta at small |t|
    r_b,_=pearsonr(delta,breadth)            # sanity: delta should track breadth
    results[label]=(delta,sdd,rho_t,p_t,r_c,p_c)
    print(f"{label:>6} | {sdd:9.3f} | {rho_t:+.3f} {p_t:6.3f} | {r_c:+.3f} {p_c:6.3f} | {r_b:+.3f}")

# figure: delta vs t for each gradient form, with binned means
fig,axes=plt.subplots(1,3,figsize=(13,4),sharey=True)
for ax,(label,(delta,sdd,rho_t,p_t,r_c,p_c)) in zip(axes,results.items()):
    ax.scatter(t,delta,s=12,alpha=0.4,color="#4C72B0")
    edges=np.linspace(t.min(),t.max(),7); idx=np.clip(np.digitize(t,edges[1:-1]),0,5)
    bx=[t[idx==b].mean() for b in range(6)]; bm=[delta[idx==b].mean() for b in range(6)]
    bs=[delta[idx==b].std(ddof=1)/np.sqrt((idx==b).sum()) for b in range(6)]
    ax.errorbar(bx,bm,yerr=bs,fmt="o-",color="#C44E52",ms=5,lw=1.5,zorder=3)
    ax.axhline(0,color="grey",lw=0.5)
    ax.set_title(f"{label} gradient\nrho(delta,t)={rho_t:+.2f} (p={p_t:.2f})\ncorr(delta,-t^2)={r_c:+.2f} (p={p_c:.2f})",fontsize=9)
    ax.set_xlabel("ML position (medial \u2192 lateral)")
axes[0].set_ylabel("Estimated broadcast propensity  \u03b4\u1d62")
plt.tight_layout()
plt.savefig(OUT_DIR / "delta_spatial_test.png", dpi=150, bbox_inches="tight")
plt.savefig(OUT_DIR / "delta_spatial_test.pdf", bbox_inches="tight")
print("\nSaved delta_spatial_test.png / .pdf")

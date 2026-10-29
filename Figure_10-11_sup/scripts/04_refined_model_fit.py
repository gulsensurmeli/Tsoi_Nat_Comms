# ============================================================================
# 04_refined_model_fit.py
# ROLE   : Fits and evaluates the refined generative model (per-network ML
#          territories + continuous per-neuron propensity) under this
#          pipeline's alternate PM(7)/AT(7) area definition (see
#          ANALYSIS.md). Mirrors the refined-model logic used by the
#          manuscript's primary PM(7)/AT(13) analysis (a separate pipeline,
#          not included in this repo).
# INPUT  : data/master_detailed_comment.csv
# OUTPUT : outputs/refined_model_fit.{png,pdf}
# NOTES  : Refined-model params (per-network logistic ML territories + continuous
#          propensity sigma) are HARDCODED from a separate fit -- see ANALYSIS.md
#          and scripts/fit_refined_params.py. The ML-only reference
#          model's KAPPA is also script 01's fitted graded-model kappa under
#          this area definition.
# ============================================================================

"""
refined_model_fit.py  (refined model, alternate PM/AT definition)

Panel G: a refined generative model reproduces the target-number, axon-length,
and breadth observations that position-only models could not.

REFINED MODEL
-------------
  e_ij = sigmoid( a_j + f_net(j)(t_i) + delta_i )
    a_j       per-area intercept (recalibrated so marginal projection
              prevalence per area matches the data)
    f_PM(t) = A_PM * sigmoid(-k_PM * (t - c_PM))    PM "territory": full-strength
    f_AT(t) = A_AT * sigmoid( k_AT * (t - c_AT))    medially, soft ML boundary.
              Each network occupies an ML territory with its own edge (c) and
              boundary sharpness (k). The AT edge sits medial of the PM edge
              (c_AT < c_PM), so the territories OVERLAP centrally: only central
              neurons fall within both and project to both networks -> broadest.
    delta_i   per-neuron breadth propensity, delta_i ~ Normal(0, sigma^2): a
              continuous over-dispersion term (a subset of neurons project more
              widely). NOT spatially organised along ML (script 05).

This replaces the earlier symmetric-gradient + binary-broadcaster model. The
symmetric gradient forced PM and AT to mirror each other and could not hold PM's
plateau into the centre, so it systematically UNDER-predicted central breadth;
the binary broadcaster was the only dispersion source and, strong enough to
match the variance, produced a bimodal breadth distribution (the data are
unimodal). Per-network territories reproduce the breadth-vs-ML profile, and the
continuous propensity matches breadth variance while keeping breadth unimodal.

INTERPRETATION NOTE (shown, not hidden)
---------------------------------------
The central "both RSP and ORB" group is broadest because RSP and ORB sit at
opposite ML poles and only central neurons fall within both network territories
-- an anchor / territory-overlap effect, NOT position-dependent broadcasting.
The propensity delta_i is spatially uniform; script 05 shows breadth propensity
is not organised along the ML axis.

NOTE ON FRAMING: broad projectors are the upper tail of a continuous propensity,
not a discrete "broadcaster" cell type -- consistent with the unimodal breadth
distribution. Describe them as "a subset of neurons that project more widely",
not as a distinct subpopulation.
"""

import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import numpy as np
from numpy.polynomial.hermite import hermgauss
from matplotlib.patches import Patch
from scipy.stats import spearmanr

from common import (
    OUT_DIR, AREAS, A, ispm, N,
    sig, logit, Y as Yr, L as Lr, t, a as a0,
    rsp_idx as rsp_i, orb_idx as orb_i, other_idx as other,
    KSHAPE,
)

# sgn: per-area sign used only by the "ML-only" comparison model in simulate()
# below (+1 for AT areas, -1 for PM areas), so one shared coefficient moves PM
# and AT probability in opposite directions from a single symmetric gradient.
# R_prev: real per-area marginal projection prevalence (fraction of the 160
# neurons projecting to each area) -- the target recalibrated_a() below tries
# to match once delta_i's effect on the mean is accounted for.
sgn    = np.where(ispm, -1.0, 1.0)
R_prev = Yr.mean(0)

# Fitted refined-model params (HARDCODED from a separate fit -- see ANALYSIS.md).
# Per-network logistic ML "territories": each network is full-strength on its
# own side of a soft ML edge (c) with boundary sharpness (k). k is capped at 4
# here (the smallest cap within tolerance for this 14-area definition -- see
# fit_refined_params.py's KMAX sweep) to avoid an overfit near-step
# boundary.
# CAUTION: these 6 numbers are NOT derived from any code in this file -- they
# were pasted in from a separate offline fit. If you change _gradient(),
# recalibrated_a(), or simulate() below, these will silently go stale; you'd
# need to re-fit and re-paste, not assume they still apply. Re-derive with
# `python scripts/fit_refined_params.py`, which re-runs the fit
# end-to-end and prints a copy-paste block plus a PASS/DIFFERS check against
# these values.
A_PM, K_PM, C_PM = 5.869, 4.000, -0.009   # PM territory: strength, sharpness, ML edge
A_AT, K_AT, C_AT = 8.316, 4.000, -0.421   # AT edge sits medial of PM's -> central overlap
SIGMA            = 1.513                   # SD of per-neuron breadth propensity delta_i

# Gauss-Hermite quadrature nodes/weights: lets recalibrated_a() below compute
# the exact average of a function over delta_i ~ Normal(0, SIGMA^2) using 16
# fixed points, instead of drawing thousands of random delta_i samples --
# noise-free, and needed repeatedly inside that function's iteration loop.
_GHx, _GHw = hermgauss(16); _GH_W = _GHw / np.sqrt(np.pi)

def _gradient():
    """Per-network ML logit shift f_net(t); shape (N, A).

    f_pm(t) is near its max (A_PM) for medial t (well below C_PM) and falls off
    past C_PM; f_at(t) is the mirror shape, near its max (A_AT) for t above
    C_AT. Because C_AT < C_PM, there is a middle band of t where BOTH curves
    are simultaneously near their maxima -- this is the territory-overlap zone
    (see module docstring). np.where then broadcasts f_pm(t)/f_at(t) (each
    length N) across all 20 areas: PM areas get f_pm(t_i) added to their logit
    for neuron i, AT areas get f_at(t_i). Depends only on real, fixed ML
    positions (t) and the hardcoded params above -- so G_REF itself never
    changes across simulate() calls or Monte Carlo reps; only the stochastic
    pieces (Y draws, delta_i) do.
    """
    f_pm = A_PM * sig(-K_PM * (t - C_PM))
    f_at = A_AT * sig( K_AT * (t - C_AT))
    return np.where(ispm[None, :], f_pm[:, None], f_at[:, None])
G_REF = _gradient()

def recalibrated_a():
    """Per-area intercepts so marginal prevalence matches data (delta integrated out).

    Why this is needed: simulate() adds a per-neuron random propensity delta_i
    to every area's logit. Because sigmoid is nonlinear, the AVERAGE (over
    delta_i) projection probability for an area is not the same as sigmoid of
    the average logit -- so plugging the original intercepts a0 straight in
    would make each area's simulated prevalence drift away from R_prev. This
    loop is a fixed-point correction: for each area, average sig(a + G_REF +
    noise) over the 16 Gauss-Hermite nodes (standing in for delta_i ~
    Normal(0, SIGMA^2)) and over all neurons (.mean(0)); compare that average
    to the real prevalence R_prev in logit space, and nudge `a` a fraction
    (0.7) of the gap. 12 iterations converges empirically; the step size and
    iteration count aren't derived analytically, just tuned to converge
    cleanly for these parameter values.
    """
    a = a0.copy()
    for _ in range(12):
        pmean = np.zeros(A)
        for q in range(len(_GHx)):
            pmean += _GH_W[q] * sig(a[None, :] + G_REF + np.sqrt(2) * SIGMA * _GHx[q]).mean(0)
        a += (logit(np.clip(R_prev, 0.02, 0.98)) - logit(np.clip(pmean, 0.02, 0.98))) * 0.7
    return a
A_CAL = recalibrated_a()

def simulate(refined, LMAX, rng):
    """Draw one stochastic replicate dataset (Y = who-projects-where, L = axon lengths).

    refined=True:  the full model -- fixed per-network territory shift (G_REF)
                   + recalibrated intercepts (A_CAL) + a FRESH random draw of
                   this replicate's own per-neuron propensity delta_i.
    refined=False: the OLD "ML-only" comparison model -- one symmetric graded
                   gradient (kappa=2.261, pinned from the earlier fit in scripts
                   17/18) with the ORIGINAL, uncalibrated intercepts a0, and no
                   propensity term at all. Plotted in every panel purely as the
                   pre-refinement baseline, not as a serious rival hypothesis.
    LMAX sets the scale of the Gamma-distributed axon length (fit separately
    per variant by calib() below, so both variants' total length matches real
    data and only the SHAPE of the group/ML profile is being compared).
    """
    if refined:
        delta = rng.normal(0.0, SIGMA, N)                      # per-neuron propensity
        e = sig(A_CAL[None, :] + G_REF + delta[:, None])
    else:  # ML-only reference (symmetric graded gradient, base intercepts)
        e = sig(a0[None, :] + (2.261 * t)[:, None] * sgn[None, :])
    Y = (rng.random((N, A)) < e).astype(float)
    L = rng.gamma(KSHAPE, (LMAX * np.clip(e, 1e-4, None)) / KSHAPE) * Y
    return Y, L

# Splits neurons into the 3 RSP/ORB groups used throughout this pipeline
# (RSP-only, ORB-only, both) -- same split as script 01/02's group_labels(),
# just returned as boolean masks here instead of string labels.
def groups(Y): rs=Y[:,rsp_i]>0; ob=Y[:,orb_i]>0; return [rs&~ob,ob&~rs,rs&ob]

real_overall_len=Lr[:,other].sum(1).mean()
# calib(): finds the single multiplicative length scale (LMAX) that makes this
# model variant's average total axon length (over the "other 12" areas,
# averaged over 30 stochastic reps to smooth out noise) equal the real data's
# average total length. Run once per variant (refined vs ML-only) so both are
# compared to data on the same length scale -- this is a scale calibration,
# not something the model is being tested on.
def calib(refined):
    rng=np.random.default_rng(1); base=np.mean([simulate(refined,1.0,rng)[1][:,other].sum(1).mean() for _ in range(30)])
    return real_overall_len/base
Lref,Lml=calib(True),calib(False)

# Monte Carlo engine for this figure's uncertainty bands: simulates `nrep`
# (600) independent replicate datasets of one model variant and summarizes
# them into the three quantities each figure panel needs.
#   NT[g] / TL[g]: for RSP/ORB group g, one number per replicate -- that
#     replicate's mean "other 12" target count / total length for neurons in
#     that group (same groups() split as above, now applied to simulated Y/L).
#   BR: one array per replicate of every neuron's total breadth (all 14
#     areas), later binned by ML position.
#   bd(): collapses one group's 600 replicate means into a [2.5th, 50th,
#     97.5th] percentile triple -- this is what gets drawn as the shaded
#     box/line in panels 0-1 of the figure.
#   bml: for each of 6 ML bins, breadth is averaged first within a replicate's
#     matching neurons, then averaged again across the 600 replicates -- the
#     model's predicted breadth-vs-ML curve, using the same binning as the
#     real-data curve computed below so the two overlay directly (panel 2).
def band(refined,LMAX,nrep=600):
    rng=np.random.default_rng(7); NT=[[],[],[]]; TL=[[],[],[]]; BR=[]
    for _ in range(nrep):
        Y,L=simulate(refined,LMAX,rng); gm=groups(Y)
        for g in range(3):
            NT[g].append(Y[gm[g]][:,other].sum(1).mean() if gm[g].sum() else np.nan)
            TL[g].append(L[gm[g]][:,other].sum(1).mean() if gm[g].sum() else np.nan)
        BR.append(Y.sum(1))
    def bd(x): x=np.array(x); x=x[~np.isnan(x)]; return np.percentile(x,[2.5,50,97.5])
    BR=np.array(BR); edges=np.linspace(t.min(),t.max(),7); idx=np.clip(np.digitize(t,edges[1:-1]),0,5)
    bml=[np.nanmean([BR[r][idx==b].mean() for r in range(nrep)]) for b in range(6)]
    return [bd(NT[g]) for g in range(3)],[bd(TL[g]) for g in range(3)],np.array(bml)
# r = refined-model bands, m = ML-only reference-model bands; each variant
# uses its own length-scale calibration (Lref / Lml) from calib() above.
nt_r,tl_r,bml_r=band(True,Lref); nt_m,tl_m,bml_m=band(False,Lml)

# real: the same statistics computed directly on the real data (not
# simulated), for the black data points/error bars plotted in every panel.
gm_r=groups(Yr)
def mci(v): v=np.asarray(v,float); return v.mean(),1.96*v.std(ddof=1)/np.sqrt(len(v))  # mean, 95% CI half-width
rnt=[mci(Yr[m][:,other].sum(1)) for m in gm_r]; rtl=[mci(Lr[m][:,other].sum(1)) for m in gm_r]
edges=np.linspace(t.min(),t.max(),7); idx=np.clip(np.digitize(t,edges[1:-1]),0,5)  # same 6 ML bins as band() above
bctr=[t[idx==b].mean() for b in range(6)]; rbml=[Yr[idx==b].sum(1).mean() for b in range(6)]
rbsem=[Yr[idx==b].sum(1).std(ddof=1)/np.sqrt((idx==b).sum()) for b in range(6)]

plt.rcParams.update({"font.family":"sans-serif","font.size":8,"axes.linewidth":0.6,
    "xtick.major.width":0.6,"ytick.major.width":0.6,"xtick.major.size":2.5,"ytick.major.size":2.5})
REAL="#1a1a1a"; MOD="#2166AC"; ML="#9a9a9a"; grp=["RSP\nonly","ORB\nonly","both"]
fig,axes=plt.subplots(1,3,figsize=(7.4,2.5))

# Draws one of the two group-comparison panels (target number or length):
# rvals = real data's [(mean, 95%CI), ...] per group; bmod = this variant's
# [(2.5,50,97.5), ...] percentile band per group (drawn as a shaded box + a
# median line); bml_line = the OTHER (ML-only) variant's band, drawn as just a
# flat grey line at its median so both models' predictions sit on one axis
# alongside the real data point. `scale` divides values (used to convert
# length from um to mm in panel 1 below).
def gpanel(ax,rvals,bmod,bml_line,ylab,scale=1.0):
    for g in range(3): ax.plot([g-0.30,g+0.30],[bml_line[g][1]*scale]*2,color=ML,lw=1.2,zorder=1)
    for g in range(3):
        lo,md,hi=[v*scale for v in bmod[g]]
        ax.add_patch(plt.Rectangle((g-0.20,lo),0.40,hi-lo,facecolor=MOD,alpha=0.22,edgecolor="none",zorder=2))
        ax.plot([g-0.20,g+0.20],[md,md],color=MOD,lw=1.4,zorder=3)
    for g in range(3):
        m,ci=[v*scale for v in rvals[g]]; ax.errorbar(g,m,yerr=ci,fmt="o",color=REAL,ms=4,capsize=2.2,lw=1,zorder=4)
    ax.set_xticks([0,1,2]); ax.set_xticklabels(grp,fontsize=7); ax.set_ylabel(ylab,fontsize=8)
    ax.spines[["top","right"]].set_visible(False); ax.margins(x=0.16)
# Panel 0: target number by RSP/ORB group -- data vs refined-model band vs
# ML-only-model line. This is the panel that shows the refined model's
# both-group undershoot (see the printed summary at the bottom of this file).
gpanel(axes[0],rnt,nt_r,nt_m,"Number of other targets")
# Panel 1: same comparison for total axon length; um -> mm via the /1e3 scale.
gpanel(axes[1],[(m/1e3,c/1e3) for m,c in rtl],[b/1e3 for b in tl_r],[b/1e3 for b in tl_m],"Total axon length (mm)")

# Panel 2: breadth-vs-ML-position curve -- ML-only model, refined model, and
# real data (binned mean +/- SEM) overlaid on one axis. This is the panel that
# shows the refined model reproducing the mid-axis breadth profile that the
# ML-only model (a flat/monotonic single gradient) cannot.
ax=axes[2]
ax.plot(bctr,bml_m,color=ML,lw=1.5,label="ML-only",zorder=1)
ax.plot(bctr,bml_r,color=MOD,lw=1.5,label="Refined",zorder=2)
ax.errorbar(bctr,rbml,yerr=rbsem,fmt="o",color=REAL,ms=4,capsize=2.2,lw=1,label="Data",zorder=4)
ax.set_xlabel("ML position (medial \u2192 lateral)",fontsize=8); ax.set_ylabel("Breadth (# targets)",fontsize=8)
ax.spines[["top","right"]].set_visible(False)
ax.legend(fontsize=6.5,frameon=False,loc="upper center",ncol=3,handlelength=1.2,columnspacing=0.9)

# Manual legend (placed on panel 0 only) spelling out what the box/line/marker
# styles mean, since gpanel() itself draws them without per-call labels.
leg=[plt.Line2D([0],[0],marker="o",color=REAL,lw=0,ms=4,label="Data (mean \u00b1 95% CI)"),
     Patch(facecolor=MOD,alpha=0.22,label="Refined model (95%)"),
     plt.Line2D([0],[0],color=ML,lw=1.2,label="ML-only model")]
axes[0].legend(handles=leg,fontsize=5.8,frameon=False,loc="upper left",handlelength=1.1,labelspacing=0.3)

fig.tight_layout(w_pad=1.3)
fig.savefig(OUT_DIR / "refined_model_fit.png", dpi=300, bbox_inches="tight")
fig.savefig(OUT_DIR / "refined_model_fit.pdf", bbox_inches="tight")
# Print the both-group numbers explicitly on every run so the documented
# undershoot (refined model ~13.4 vs real ~14.1, see module docstring) stays
# visible rather than getting silently lost in the figure alone.
print("saved. group n-targets [RSP-only, ORB-only, both] --",
      "data:", [round(float(m), 1) for m, _ in rnt],
      "refined:", [round(float(b[1]), 1) for b in nt_r],
      "(both-group is a documented undershoot)")

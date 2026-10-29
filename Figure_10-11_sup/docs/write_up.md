# Alternate PM/AT area-group definition: results write-up

## Introduction

The manuscript's primary projection-modeling analysis (a separate pipeline,
not included in this repository) splits ECL5a projection neurons into
posteromedial (PM, 7 areas) and anterior-temporal (AT, 13 areas) networks,
using RSP and ORB as flagship anchors for the "broadcast" phenotype. This
repository asks a scoping question: do the same conclusions hold under a
different, smaller, and more symmetric curated PM/AT grouping? PM is
unchanged as a set (VISpor, RSP, PAR, PRE, PTLp, VIS, POST); AT is replaced
with a 7-area set (PERI, ORB, ECT, TEa, AI, PFC, OLF) that drops GU, VISC,
SS, ACA, and MO outright and folds PL and ILA into the single, native `PFC`
column (`PFC_length == PL_length + ILA_length`, verified to float
precision) rather than modeling them separately. This halves AT's area
count (13 to 7) and brings PM and AT to parity (7 and 7), at the cost of
dropping five areas from the analysis entirely and losing PL/ILA's
individual resolution.

Scripts 01-06 mirror the primary pipeline's model family exactly — same
eleven-model family, same cross-validation and posterior-predictive checks,
same broadcast-trait extension, same refined model, same spatial and
territory tests — against this 14-area definition (12 "other" areas,
excluding RSP/ORB, versus the primary pipeline's 18). Script 07 asks a
different, more granular question: not group-level statistics, but each
neuron's exact present/absent motif across a curated 8-area subset. Two
fit-specific constants needed re-deriving under the new definition: script
03's `KAPPA` (read off script 01's fitted graded-identity slope) and script
04's territory/propensity constants (via `scripts/fit_refined_params.py`).
This document summarizes what was found, focusing on where the alternate
definition replicates the primary pipeline's conclusions and where it
diverges. Per-figure numeric detail lives in `figure_legends.md`; see
`ANALYSIS.md` for the mechanical setup (data contract, area lists, shared
setup module).

## Results

### Fitting eleven candidate models under the 14-area definition (script 01)

The same eleven models were fit to smooth PM/AT probability curves and
simulated (**Figure `sim_models_fitted_grid.png`**). The graded identity
model fits kappa = 2.261 (used downstream in scripts 03 and 04), compared
to the primary definition's ~2.5-ish fit — a plausible shift given the
smaller, more balanced AT set. The graded-identity/choose-then-grow
degeneracy break again succeeds: real data give a small nonzero conditional
-length slope (coef = 0.073, SE = 0.010); graded identity reproduces a
comparable nonzero slope (0.409, SE 0.023) while choose-then-grow's slope is
near zero (−0.036, SE 0.026) — the same qualitative separation as the
primary pipeline, confirming the mechanism-level distinction survives the
area-set change.

### Cross-validated comparison replicates the two held-out failures (script 02)

Five-fold CV and posterior-predictive checks again show ML position is
real, generalizable information (**Figure `sim_models_comparison.png`**):
null (−1442.5) and area-intrinsic (−1429.1) trail every position-dependent
model, with wiring-cost/geometry (−1162.8), broadcast/expanded (−1165.1),
and two independent drives (−1167.5) leading. The two held-out failures
that motivate the rest of the pipeline reproduce as well: Spearman
rho(breadth, ML) is real = 0.017 here (much lower than the primary
definition's 0.20, consistent with a flatter, more symmetric breadth
profile at 14 areas), with most position-dependent models overshooting into
the upper percentiles (84th-98th); and the RSP/ORB "both minus single"
contrast (real = 2.93 over the other-12-area scale) is missed by every
model, landing at or near the 100th percentile — the same qualitative
pattern as the primary pipeline, on a rescaled axis.

### A position-dependent broadcast boost again recovers the group contrast (script 03)

The same three nested variants were built on the graded model (kappa =
2.261, re-derived from script 01) (**Figure `sim_broadcast_trait.png`**).
Variant B (uniform trait alone, sigma = 0.945) improves rho's percentile
somewhat (84.5th) but still fails the contrast (100th percentile). Variant
C (position-dependent boost: amp = 4.256, width = 0.301, plus a residual
uniform trait sigma = 1.015) is again the first to reproduce the group
contrast, landing at the 55.0th percentile, while its rho percentile
(99.0th) sits further into the tail than the primary definition's
equivalent (95.8th) — the position-dependent mechanism still does the
necessary work for the contrast, but the residual uniform-trait component
is doing proportionally more of variant C's overall fit here than under the
primary definition, where the positional boost dominated (sigma was only
0.069 there, versus 1.015 here).

### The refined model reproduces the signature, with a smaller residual gap (script 04)

The refined model (per-network ML territories plus a continuous per-neuron
propensity) was refit end-to-end for this definition (**Figure
`refined_model_fit.png`**). Fitted parameters: PM territory A_PM = 5.869,
k_PM = 4.000 (capped — the smallest cap within tolerance here, versus the
primary pipeline's cap of 6), c_PM = −0.009; AT territory A_AT = 8.316,
k_AT = 4.000 (capped), c_AT = −0.421 — the AT edge again sits medial of the
PM edge, producing the same central territory overlap; propensity SD sigma
= 1.513. The model recovers RSP-only (data 7.1, model 7.1) and ORB-only
(data 6.6, model 7.0) target counts closely, and closes most of the gap on
the both-group, though it still undershoots it (model 9.4 vs. real 9.8) —
proportionally smaller than the primary pipeline's undershoot (~5% there,
~4% here). The same pattern holds for axon length: real data (RSP-only 46.7
mm, ORB-only 73.6 mm, both 81.8 mm) is matched well on the both-group by the
refined model's median (82.7 mm) but undershot on the two single-network
groups (57.8 mm, 59.2 mm), while the flat ML-only reference (~66-68 mm
across all groups) captures none of the group structure. The
breadth-vs-ML panel again shows the refined model tracking the real
mid-axis peak that the ML-only reference misses.

### Is the breadth propensity itself spatially organized? (script 05)

The same delta_i spatial test was rerun at three gradient flexibilities
(**Figure `delta_spatial_test.png`**). Under a **linear** gradient, delta_i
again shows a nominally significant central concentration (corr(delta,
−t²) = +0.160, p = 0.043) — weaker than the primary definition's
linear-gradient signal (+0.306, p < 0.001) but still present. This
**vanishes** at higher flexibility: **quadratic** (+0.016, p = 0.843) and
**cubic** (+0.017, p = 0.835) both give a null result — the same
gradient-misfit-not-genuine-spatial-structure conclusion as the primary
pipeline. No monotonic trend appears at any flexibility (rho(delta, t) =
+0.018 to +0.020, all p > 0.79), and delta_i again correlates strongly with
breadth (r = 0.92-0.97), confirming it remains a valid propensity measure.

### Testing the territory asymmetry directly (script 06)

The per-network vs. reflection-symmetric territory test was rerun directly
against the 14-area data (**Figure `territory_test.png`**, plus standalone
**`territory_asymmetric_fit.png`**). The per-network model again fits
significantly better (chi2(3) = 20.03, p = 1.67×10⁻⁴; log-likelihood
−1048.63 vs. −1058.65), with AIC agreeing (dAIC = +14.0) and BIC essentially
a wash (dBIC = −3.1) — the same qualitative pattern as the primary pipeline
(dAIC = +18.0, dBIC = −0.2). The central co-occupancy zone is robust under
bootstrap (600 resamples): median width 0.667 ML units, 95% CI [0.537,
0.826], present in 100% of resamples — as robust as the primary pipeline's
result (0.607, CI [0.457, 0.766], also 100%). Held-out per-neuron
log-likelihood again shows only a negligible, non-significant advantage for
the per-network model (symmetric −6.723 vs. per-network −6.720; paired
Δ = +0.003 ± 0.044 SE) — the same modest, by-design gap as the primary
pipeline (+0.016 ± 0.060), consistent with the asymmetry's value lying in
the population-level breadth profile, not per-neuron prediction.

### The refined model struggles with one-network-only motifs (script 07)

The 8-area combinatorial motif-frequency figure (RSP, VIS, PTLp, PERI, AI,
TEa, ORB, PFC) (**Figure `motif_combination_frequency.png`**) surfaces a
mismatch distinct from the group-level statistics above: the refined
model badly **underpredicts the all-AT/no-PM motif** (PERI+AI+TEa+ORB+PFC
present, none of RSP/VIS/PTLp: real count 21, model median 10, 95% band
[5, 16] — the real value falls outside the band), while **reproducing the
all-8-present motif closely** (real 22, model median 25, band [17, 34]).
This is the same qualitative "struggles with one-network-only engagement"
pattern documented for the primary pipeline's own combinatorial
motif-frequency analysis, on the opposite network here (all-AT rather than
all-PM) — plausibly related to the same territory-amplitude effect that
causes the refined model to over-predict neurons reaching every PM area
simultaneously in the primary pipeline (a shared per-network territory
amplitude, applied uniformly to every area's intercept, cannot represent an
area that stays persistently under- or over-recruited regardless of
position). Not yet formally diagnosed for this area definition; flagged
here for visibility rather than smoothed over.

## Discussion

**Main conclusions.** Every headline qualitative conclusion from the
primary 20-area pipeline replicates under this smaller, more symmetric
14-area definition: ML position organizes network identity; no ML-only
model reproduces the broadcast breadth signature; a position-dependent
broadcast term (not a uniform trait) is specifically what recovers the
RSP/ORB group contrast; the refined model's per-network territories plus
continuous propensity jointly reproduce the target-number, axon-length, and
breadth-vs-ML signatures (with a residual, proportionally smaller
both-group undershoot); the propensity term carries no genuine spatial
signal once gradient misfit is ruled out; the territory asymmetry is a
real, bootstrap-robust, profile-justified feature rather than a
per-neuron-predictive one; and the refined model struggles with
one-network-only combinatorial motifs (here: all-AT), mirroring the primary
pipeline's own combinatorial mismatch on the opposite network. None of
these conclusions depends on the specific 20-area PM/AT membership choice —
they hold under a considerably different, smaller area set that also
changes PM/AT from a 7-vs-13 to a 7-vs-7 split.

**Where the two definitions diverge quantitatively.** The breadth-ML
Spearman rho is much lower here (real = 0.017 vs. 0.20), consistent with a
flatter breadth-vs-ML profile at 14 areas; script 03's fitted broadcast
variant C relies proportionally more on its residual uniform trait
(sigma = 1.015) than on the positional boost, compared to the primary
definition where the positional boost dominated and the uniform trait was
almost negligible (sigma = 0.069) — the mechanism split between "where" and
"how much" versus "who broadcasts more in general" shifts somewhat under
the smaller area set, though the position-dependent term remains necessary
either way. The refined model's fitted territory sharpness caps at a lower
boundary (KMAX = 4 vs. 6 in the primary pipeline), and its both-group
undershoot is proportionally smaller (~4% vs. ~5%) — plausibly because
there are fewer "other" areas to under- or over-predict across, so profile
mismatches average out differently.

**Limitations, carried over from the primary pipeline.** As with the
primary pipeline's refined model, this repository's script 04 parameters
are hardcoded from a separate offline fit (`scripts/fit_refined_params.py`,
PASS-confirmed against the pasted values in script 04) and will go stale if
the model or its inputs change without re-running that utility. The
territory asymmetry's statistical support (LRT, bootstrap) again outpaces
its held-out per-neuron predictive value, for the same by-design reason as
the primary pipeline. And, as throughout this pipeline, fixed random seeds
mean small run-to-run numeric drift is expected and not a bug. One
limitation specific to this analysis: dropping GU, VISC, SS, ACA, and MO,
and merging PL/ILA into a single PFC area, was not itself validated against
any alternative — it was taken as given, per the requested area list, and
this write-up does not attempt to judge whether that specific curation is
biologically preferable to the primary 20-area definition, only whether the
pipeline's conclusions are robust to making it.

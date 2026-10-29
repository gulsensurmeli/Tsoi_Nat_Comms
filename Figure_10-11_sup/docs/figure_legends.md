# Figure legends: alternate PM/AT area-group definition (scripts 01-07)

*Legends for the figures produced by this repository's scripts (see
`ANALYSIS.md`). This pipeline reruns the manuscript's primary projection-
modeling analysis (a separate pipeline, not included here) — same eleven-
model family, same cross-validation and posterior-predictive checks, same
broadcast-trait extension, same refined model, same spatial and territory
tests — against a smaller, symmetric curated PM(7)/AT(7) area grouping:
VISpor, RSP, PAR, PRE, PTLp, VIS, POST (PM) and PERI, ORB, ECT, TEa, AI, PFC,
OLF (AT), instead of the manuscript's primary PM(7)/AT(13) grouping. GU,
VISC, SS, ACA, and MO are excluded outright; PL and ILA are not excluded but
are folded into the single, native `PFC` column (`PFC_length == PL_length +
ILA_length`, verified to float precision). Fitted values and test statistics
are reported here as observed from a run of each script; per `ANALYSIS.md`'s
seed-dependence note, exact values can drift slightly run to run without
changing any conclusion.*

---

## `sim_models_fitted_grid.png` / `.pdf` (script 01)

**Eleven-model fit/simulate/plot grid under the 14-area alternate
definition.** Rows are real data followed by the eleven candidate models,
columns are (1) PM/AT bias vs ML position with Spearman rho annotated, (2)
breadth (out of 14 areas) binned by ML position, (3) mean per-area
projection probability by network, and (4-5) the RSP/ORB broadcast-group
contrast in number of other targets and total axon length (over the other
12 areas). Real data's bias-ML correlation is strong under this definition.
The graded-identity model fits kappa = 2.261 (used downstream as script
03/04's `KAPPA`/reference-model slope). The conditional-length degeneracy
breaker separates graded identity from choose-then-grow: real data coef =
0.073 (SE 0.010); graded identity reproduces a comparable nonzero slope
(0.409, SE 0.023) while choose-then-grow's slope is near zero and slightly
negative (−0.036, SE 0.026), as intended by construction. (Network-then-
area's slope is undefined — too few conditional observations survive its
degenerate fit under this smaller area set to estimate it.)

---

## `sim_models_comparison.png` / `.pdf` (script 02)

**Cross-validated comparison: ML position generalizes, but no position-only
model reproduces the broadcast contrast.** (A) Five-fold CV held-out
log-likelihood, sorted descending: null (−1442.5) and area-intrinsic
(−1429.1) trail every position-dependent model by 250-1000+ log-likelihood
units; wiring-cost/geometry (−1162.8), broadcast/expanded (−1165.1), and two
independent drives (−1167.5) form the top tier; trade-off allocation
(−1400.9) and network-then-area (−1605.8, worst) lag furthest behind. (B)
Posterior-predictive violins for Spearman rho(breadth, ML), real = 0.017
(much lower than the primary 20-area definition's ~0.20, reflecting the
smaller, more symmetric 14-area profile). Null/area-intrinsic sit
mid-distribution (57-59th percentile) while most position-dependent models
overshoot toward the upper tail (84-98th percentile) — trade-off allocation
is the exception, undershooting at the 0.0th percentile. (C)
Posterior-predictive violins for the RSP/ORB "both minus single" breadth
contrast, real = 2.93 (other-12-area scale). Every model sits at or near the
100th percentile — no ML-only model reproduces the real broadcast contrast,
the failure that motivates script 03's broadcast-trait extension.

---

## `sim_broadcast_trait.png` / `.pdf` (script 03)

**Adding a broadcast term recovers the group contrast that no ML-only model
could reach, though the position-dependent boost does less work here than
under the primary 20-area definition.** Three nested variants built on the
graded-identity model with kappa = 2.261 (re-derived from script 01):
**A. ML-only (graded)**; **B. + uniform trait**, fit sigma = 0.945; **C. +
position-dependent broadcast**, jointly fit amp = 4.256, width = 0.301,
sigma = 1.015 (a much larger residual uniform-trait sigma than variant C's
primary-definition fit, where the positional boost did more of the work).
Left panel: Spearman rho(breadth, ML), real = 0.017. Variant A overshoots at
the 95.2nd percentile; variant B is intermediate (84.5th); variant C
overshoots furthest (99.0th). Right panel: the RSP/ORB "both minus single"
contrast, real = 2.93. Variant A (100th percentile) and variant B (100th
percentile) both fail to reach it; variant C reaches the 55.0th percentile —
comfortably centered, showing that a position-dependent broadcast component
(not a uniform trait alone) is needed to reproduce the contrast, even though
the uniform-trait component does comparatively more of variant C's total
work under this smaller area set.

---

## `refined_model_fit.png` / `.pdf` (script 04)

**The refined model reproduces the group target-number, axon-length, and
breadth-vs-ML signatures under the alternate 14-area definition, with a
proportionally smaller undershoot than the primary 20-area definition.**
Real data (black points, mean ± 95% CI) vs. 600-replicate percentile bands
for the refined model (blue) and the ML-only reference model (grey line,
kappa = 2.261). (Left) Number of other targets (12) by RSP/ORB group. Real
data: RSP-only 7.1, ORB-only 6.6, both 9.8 (n = 57/48/55). The refined model
recovers RSP-only (7.1) and ORB-only (7.0) almost exactly and closes most of
the gap on the both-group (9.4) relative to the flat ML-only reference
(~7.6 across all three groups), though it still undershoots the real
both-group value (9.4 vs. 9.8), proportionally smaller (~4%) than the
primary definition's undershoot. (Middle) Total axon length over the other
12 areas by group, rescaled to mm. Real data: RSP-only 46.7 mm, ORB-only
73.6 mm, both 81.8 mm. The refined model's median simulated length is 57.8
mm (RSP-only), 59.2 mm (ORB-only), 82.7 mm (both) — it matches the
both-group closely but undershoots RSP-only/ORB-only somewhat, while the
flat ML-only reference (~66-68 mm across all groups) captures none of the
group structure. (Right) Breadth (all 14 areas) binned by ML position: the
refined model's blue curve tracks the real data's mid-axis peak, while the
ML-only reference (grey) stays comparatively flat and misses it. Fitted
parameters (hardcoded from `scripts/fit_refined_params.py`, PASS-checked
against these pasted values): PM territory A_PM = 5.869, k_PM = 4.000
(capped — the smallest cap within tolerance for this area set), c_PM =
−0.009; AT territory A_AT = 8.316, k_AT = 4.000 (capped), c_AT = −0.421 —
the AT edge sits medial of the PM edge, producing central territory
overlap; propensity delta_i ~ Normal(0, sigma²) with sigma = 1.513.

---

## `delta_spatial_test.png` / `.pdf` (script 05)

**The breadth propensity's apparent central concentration is a
gradient-misfit artifact that vanishes with a more flexible gradient.**
Three-panel layout (linear, quadratic, cubic gradient flexibility). Under a
**linear** gradient, delta_i shows a significant central concentration
(corr(delta, −t²) = +0.160, p = 0.043) — weaker than the primary
definition's linear-gradient signal (+0.306, p < 0.001), but still nominally
significant. This **vanishes** at higher flexibility: **quadratic** (+0.016,
p = 0.843) and **cubic** (+0.017, p = 0.835) both give a null result. No
monotonic trend is present at any flexibility (rho(delta, t) = +0.018 to
+0.020, all p > 0.79). Delta_i correlates strongly with total breadth at
every flexibility (r = 0.92-0.97), confirming it remains a valid propensity
measure under this area definition.

---

## `territory_test.png` / `.pdf` (script 06)

**Direct test of PM/AT territory asymmetry favors the per-network model
over the reflection-symmetric one, with a "profile-justified, not
per-neuron-predictive" caveat.** Three-panel layout. (A) Fitted per-network
territories (solid) vs. the nested reflection-symmetric alternative
(dashed) against binned real projection probability, with the central
co-occupancy zone (both networks' fitted mean probability > 0.5) shaded.
The per-network model fits significantly better (LRT: chi2(3) = 20.03, p =
1.67×10⁻⁴; log-likelihood −1048.63 vs. −1058.65), and AIC agrees (dAIC =
+14.0, favoring per-network) while BIC is a wash (dBIC = −3.1). Fitted
logistic midpoints: c_PM = +0.189, c_AT = −1.518 (AT's midpoint extrapolates
outside the sampled ML range, hence the co-occupancy-zone measure);
co-occupancy width 0.667 ML units, peak joint probability 0.76. (B)
Bootstrap distribution (600 neuron-resamples): width median 0.667, 95% CI
[0.537, 0.826], with a central zone present in 100% of resamples (P =
1.000). (C) Five-fold held-out per-neuron log-likelihood, symmetric
(−6.723) vs. per-network (−6.720): a paired difference of only +0.003 ±
0.044 SE, not significant — consistent with the asymmetry's value lying in
the population-level breadth-vs-ML profile and both-group signature (script
04), not per-neuron prediction.

---

## `territory_asymmetric_fit.png` / `.pdf` (script 06)

**Standalone version of panel A above**, showing only the fitted
per-network territories against binned real data and the shaded
co-occupancy zone (no symmetric-model comparison curves), for use where a
single, uncluttered panel is preferred.

---

## `motif_combination_frequency.png` / `.pdf` (script 07)

**UpSet-style combinatorial motif-frequency figure over an 8-area subset,
surfacing a distinct model mismatch from the group-level statistics
above.** Unlike the group-level target-number/axon-length/breadth
statistics in scripts 02-06, this figure asks about each neuron's *exact*
present/absent projection pattern across a curated 8-area subset drawn from
this pipeline's PM(7)/AT(7) definition: RSP, VIS, PTLp (PM) and PERI, AI,
TEa, ORB, PFC (AT). PFC is one of this definition's own native modeled
areas (not a derived PL+ILA combination), so no OR-merge special-casing is
needed for either real or simulated data — all 8 areas are read the same
generic way. The figure plots real data's most frequent observed motifs
against the refined model's Monte Carlo predicted count (600 replicates,
median + 95% band) for those same motifs, with a shared combination-matrix
panel showing which of the 8 areas each motif involves. The model notably
**underpredicts the all-AT/no-PM motif** (PERI+AI+TEa+ORB+PFC present, none
of RSP/VIS/PTLp: real count 21, model median 10, 95% band [5, 16] — the
real value falls outside the band), while **reproducing the all-8-present
motif closely** (real 22, model median 25, band [17, 34]). This mismatch is
not yet diagnosed; it is flagged here for visibility rather than smoothed
over — see `docs/write_up.md`'s discussion section.

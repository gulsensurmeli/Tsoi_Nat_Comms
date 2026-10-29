# ECL5a projection modeling — alternate PM/AT area-group definition

Context for reproducing and evaluating this analysis. Read this before
editing any script.

## What this is

This is one analysis from a larger study of ECL5a (entorhinal cortex layer
5a) projection neurons. The scientific question: how does a neuron's
mediolateral (ML) soma position drive its projection identity across two
cortical networks — PM (posteromedial) and AT (anterior-temporal) — and how
do we account for neurons that project broadly across *both* networks (the
"broadcast" phenotype)? The pipeline fits a family of candidate generative
models to single-neuron projection data and compares them.

The manuscript's primary analysis (a separate pipeline, not included in this
repository) answers this question using a PM(7)/AT(13) area grouping. This
repository reruns the exact same model family and analyses — same eleven
candidate models, same cross-validation and posterior-predictive checks,
same broadcast-trait extension, same refined model, same spatial and
territory tests, same combinatorial motif-frequency figure — against a
different, smaller, more symmetric PM(7)/AT(7) grouping, as a robustness
check: do the primary analysis's conclusions depend on the specific 20-area
membership choice, or do they hold under a considerably different area set?
See `docs/write_up.md` for the answer, area by area and figure by figure.

## Data contract (load-bearing — every script assumes this)

- File: `data/master_detailed_comment_updated__4_.csv`.
- Analysis uses the **160-neuron PMAT subset**: rows where the `comment`
  column is null (`df[df["comment"].isna()]`). Every script asserts n == 160.
- Per-area columns are `<AREA>_length` (axon length) and `<AREA>_endpoint`.
- **PM (7 areas, unchanged from the manuscript's primary grouping):**
  VISpor, RSP, PAR, PRE, PTLp, VIS, POST.
- **AT (7 areas, this analysis's alternate grouping):** PERI, ORB, ECT, TEa,
  AI, PFC, OLF. This drops GU, VISC, SS, ACA, and MO outright (no
  substitute — genuinely excluded from this analysis) relative to the
  primary pipeline's 13-area AT set. PL and ILA are *not* dropped — their
  signal is folded into a single `PFC` column instead: `PFC` here is a
  **genuine raw-CSV column** (`PFC_length`/`PFC_endpoint`), verified
  `PFC_length == PL_length + ILA_length` to float precision for all 160
  subset rows, loaded through the exact same generic per-area machinery as
  every other area. It is fit as its own ordinary area with its own
  logistic curve, not a derived combination the model builds post hoc.
- 14 areas total are modeled (12 "other" areas once RSP and ORB are
  excluded — see below).
- "Projects to area j" is defined as `<AREA>_length > 0`.
- ML coordinate is the `z` column, standardized (z-score) then divided by
  its max abs value so `t ∈ [−1, 1]`. Convention: medial → PM, lateral → AT,
  giving a NEGATIVE rho(t, PM/AT bias). Do not flip this sign.
- The two group-level statistics (target number, total axon length) are
  computed over the **other 12 areas only** — i.e. excluding RSP and ORB
  themselves, since those two define the grouping. RSP/ORB groups:
  RSP-only, ORB-only, both (in this subset every neuron projects to RSP or
  ORB, so there is no "neither" group).

## Pipeline (run in order; each script is standalone)

Paths inside each script resolve relative to the script location, so run
from anywhere. All outputs go to `outputs/`.

1. `scripts/01_sim_models_fitted.py` — 11 models fit to data, 6-panel grid.
   Produces `sim_models_fitted_grid.{png,pdf}`. Its printed "4. Graded
   identity (1 axis)" fitted kappa is the source for script 03's hardcoded
   `KAPPA` constant.
2. `scripts/02_sim_models_comparison.py` — CV log-likelihood + posterior-
   predictive checks. Produces `sim_models_comparison.{png,pdf}`.
3. `scripts/03_sim_broadcast_trait.py` — broadcast-variant comparison on
   held-out statistics. Produces `sim_broadcast_trait.{png,pdf}`.
4. `scripts/04_refined_model_fit.py` — refined model: asymmetric per-network
   ML territories + continuous per-neuron propensity. Produces
   `refined_model_fit.{png,pdf}`.
5. `scripts/05_delta_spatial_test.py` — per-neuron delta spatial test with
   gradient-robustness guardrail. Produces `delta_spatial_test.{png,pdf}`.
6. `scripts/06_territory_test.py` — direct test of the PM/AT ML-territory
   asymmetry the refined model uses (per-network vs reflection-symmetric
   LRT; bootstrap of the central co-occupancy zone; held-out comparison).
   Produces `territory_test.{png,pdf}` and
   `territory_asymmetric_fit.{png,pdf}` (a standalone version of the fitted-
   territories panel, with no symmetric-model comparison curves).
7. `scripts/07_motif_combination_frequency.py` — UpSet-style combination-
   frequency figure over an 8-area subset: RSP, VIS, PTLp (PM) and PERI, AI,
   TEa, ORB, PFC (AT) — real data's most frequent motifs vs. the refined
   model's Monte Carlo predicted count for those same motifs, plus a shared
   combination-matrix panel. Produces `motif_combination_frequency.{png,pdf}`.

`scripts/fit_refined_params.py` is a derivation utility (not a numbered
pipeline step): it re-fits the refined model's territory/propensity
constants end to end and prints a copy-paste block plus a PASS/DIFFERS
check against the constants currently pasted into `04_refined_model_fit.py`
(and duplicated in `07_motif_combination_frequency.py`). Run it only if you
change the model or its inputs — see "Hardcoded fitted parameters" below.

## Shared setup

Shared setup lives in `scripts/common.py`, imported by every numbered
script and by `fit_refined_params.py`. It contains: path constants
(`DATA_FILE`, `OUT_DIR`), data loading + the n==160 assert, area lists
(`PM`, `AT`, `AREAS`, `ispm`, `nPM`, `nAT`), `sig`/`logit` helpers, `Y`/`L`/
`t`/`a` arrays, RSP/ORB group indices, and `LMAX`/`KSHAPE`. If you make
substantive edits (e.g. change the "projects" definition, the area lists,
or the ML scaling), change `common.py` — it propagates to every script
automatically.

## Key methodological decisions (and why)

- **Fit to breadth statistics, not the monotonic PM/AT probability curves.**
  The aggregated PM/AT-vs-ML curves are monotonic by construction (logit-
  linear GLM target), so they cannot represent the mid-axis broadcast bump.
  Script 01 fits to those curves and therefore *suppresses* the broadcast
  mechanism (the broadcast model's boost fits to ~0 there); this is
  expected, not a bug. Scripts 03–05 fit to breadth-based statistics
  instead.
- **Saturated benchmark is a refit, not a model.** In script 01 the
  "saturated benchmark" row is a per-area logistic refit of the data — a
  ceiling reference, not a mechanistic hypothesis. Do not rank it against
  the models as a peer.
- **Graded vs choose-then-grow degeneracy** is broken by conditional axon
  length (length given the neuron projects, regressed on sign-corrected ML
  with area fixed effects). Graded → nonzero slope (same sign as real);
  choose-then-grow → ~0. This is the 6th column in script 01's grid.
- **Developmental gradient** is amplitude-matched to the graded fit (single
  free steepness parameter) because scale/steepness are jointly non-
  identifiable; its steep→0 limit is the graded model.
- **delta_i is marginalized / partially pooled**, never read off a neuron's
  own projections, to avoid circularity. In script 05, delta_i is a shrunk
  random intercept; the spatial test (delta vs t) is not self-predictive
  because t only enters the shared gradient.
- **Guardrail in script 05:** any delta–ML relationship must be stable
  across linear/quadratic/cubic gradient flexibility. An apparent central
  concentration that appears only at linear flexibility and vanishes at
  quadratic/cubic is gradient misfit, not spatial broadcasting — see
  `docs/write_up.md` for this analysis's actual result.
- **The per-network ML asymmetry is real but profile-justified (script
  06).** A per-network territory can beat a reflection-symmetric one on a
  likelihood-ratio test while adding almost nothing to held-out per-neuron
  likelihood — its value is reproducing the breadth-vs-ML profile / both-
  group (a conditional-mean feature), not per-neuron prediction.

## Known limitations (carry these forward; don't "fix" them silently)

- The refined model (script 04) reproduces single-groups, breadth mean,
  variance, rho, AND the breadth-vs-ML profile, but tends to still
  undershoot the both-group target number somewhat. See `docs/write_up.md`
  for the exact numbers from this area definition. Keep this visible in any
  figure legend; don't overstate the fit as exact.
- **Broad projectors are the upper tail of a continuous propensity, not a
  discrete cell type.** The refined model uses delta_i ~ Normal(0, sigma^2)
  (continuous), which matches the breadth variance while keeping the
  simulated breadth distribution unimodal. Describe broad projectors as "a
  subset of neurons that project more widely," not as a discrete
  subpopulation or with a fixed percentage.
- **delta_i is a scalar breadth propensity.** Script 05 establishes breadth
  is not spatial; it does NOT by itself address cross-network
  directionality. That is handled by the per-network ML territories in
  script 04: the AT territory edge sits medial of the PM edge, so they
  overlap centrally, and a neuron with a spatially-uniform propensity lands
  in an intermediate "both" group only near the centre (RSP and ORB sit at
  opposite poles).
- **Hardcoded fitted parameters.** Script 03 contains a fit constant
  (`KAPPA`, read off script 01's stdout) and script 04 contains fit
  constants (the per-network territory params `A_PM/K_PM/C_PM`,
  `A_AT/K_AT/C_AT` and the propensity `SIGMA`) pasted from a separate fit;
  script 07 duplicates script 04's constants rather than importing it (it
  has no `__main__` guard, so importing it would re-run the whole fit). In
  script 04, `K_PM`/`K_AT` are capped (see the script for the current cap)
  to avoid an overfit near-step boundary. If you change the model or
  re-fit, these can silently drift out of sync with the fit they came from.
  Re-derive and verify script 04/07's constants by running
  `scripts/fit_refined_params.py` — it re-runs the fit end-to-end, prints a
  copy-paste block, and checks against the pasted values (PASS/DIFFERS).
  Script 03's `KAPPA` is not scripted; re-derive it by hand from script 01's
  stdout.
- **Seed-dependence.** Scripts set fixed seeds. Small run-to-run
  differences in fraction/undershoot magnitude/percentile values are
  expected and are not bugs.

## Environment

`pip install -r requirements.txt`. statsmodels is load-bearing for script
05 (`BinomialBayesMixedGLM`). All scripts use a non-interactive matplotlib
backend and write PNG+PDF to `outputs/`.

# ECL5a projection modeling — alternate PM/AT area-group definition

This repository reruns a projection-modeling analysis of ECL5a (entorhinal
cortex layer 5a) neurons under an alternate, smaller, more symmetric PM/AT
area grouping. It was extracted from a larger analysis codebase so that this
specific analysis can be reviewed and reproduced independently. See
`ANALYSIS.md` for the full scientific context — data contract, area
definitions, methodological decisions, and known limitations. That file is
the primary reference; read it before editing any script.

## Layout

```
ecl5a-altpmat-pipeline/
  README.md            # this file — setup and run instructions
  ANALYSIS.md           # scientific context (read this first)
  requirements.txt      # pinned Python dependencies
  data/                 # master_detailed_comment_updated__4_.csv
  scripts/              # 01..07, numbered in pipeline order, plus common.py
                         # and fit_refined_params.py
  docs/
    figure_legends.md    # per-figure legend, statistic-by-statistic
    write_up.md           # narrative results write-up
  outputs/               # figures land here (PNG + PDF), created on run
```

## Setup

```
pip install -r requirements.txt
```

## Run

Each script is standalone and resolves its own paths relative to its own
location, so you can run from anywhere:

```
python3 scripts/01_sim_models_fitted.py
python3 scripts/02_sim_models_comparison.py
python3 scripts/03_sim_broadcast_trait.py
python3 scripts/04_refined_model_fit.py
python3 scripts/05_delta_spatial_test.py
python3 scripts/06_territory_test.py
python3 scripts/07_motif_combination_frequency.py
```

or all at once, in order:

```
for f in scripts/0*.py; do echo "=== $f ==="; python3 "$f"; done
```

Each script prints its key fitted numbers/statistics to stdout and writes
PNG + PDF into `outputs/`.

### Re-deriving hardcoded fit constants (only needed if you change the model)

Scripts 03, 04, and 07 contain a few fit constants pasted from separate,
one-off fits (see `ANALYSIS.md`'s "Hardcoded fitted parameters" section for
why). These do not need to be re-run for a normal reproduction of the
figures — they're already set correctly for the current model. If you
change the model definition or the underlying data, re-derive them with:

```
python3 scripts/01_sim_models_fitted.py     # source of script 03's KAPPA
python3 scripts/fit_refined_params.py       # source of script 04/07's territory + propensity constants
```

`fit_refined_params.py` prints a copy-paste block and a PASS/DIFFERS check
against the constants currently pasted in `04_refined_model_fit.py`.

## Script → figure map

| Script | Produces (`outputs/`) |
|---|---|
| `01_sim_models_fitted.py` | `sim_models_fitted_grid.{png,pdf}` |
| `02_sim_models_comparison.py` | `sim_models_comparison.{png,pdf}` |
| `03_sim_broadcast_trait.py` | `sim_broadcast_trait.{png,pdf}` |
| `04_refined_model_fit.py` | `refined_model_fit.{png,pdf}` |
| `05_delta_spatial_test.py` | `delta_spatial_test.{png,pdf}` |
| `06_territory_test.py` | `territory_test.{png,pdf}`, `territory_asymmetric_fit.{png,pdf}` |
| `07_motif_combination_frequency.py` | `motif_combination_frequency.{png,pdf}` |

Per-figure legends (what each panel shows, and the key statistics annotated
on it) are in `docs/figure_legends.md`. A narrative results write-up,
including where this alternate area definition's conclusions agree or
diverge from the manuscript's primary analysis, is in `docs/write_up.md`.

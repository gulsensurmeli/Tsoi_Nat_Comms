# Figure 7

This folder contains all the files needed to generate plots and analysis results related to Figure 7.

To set up the environment with the appropriate packages to run these scripts, set up the working directory in this folder and run:

```
conda env create -f set_up.yml
conda activate ec_analysis_fig_7
```

The files need to be run in the following order, as they rely on the datasets generated in earlier stages:

1. **`bulk_transcriptomics_fig7b-c.ipynb`** - generates **`bulk_gene_classes.csv`** with bulk sequencing results assigned to gene expression pattern groups. It also generates Figure 7 panels b-c.
2. **`merfish_gene_trend_analysis_fig7e.ipynb`** - generates **`merfish_final_gene_categories.csv`** with MERFISH data grouped into expression patterns across the entorhinal cortex mediolateral axis. Additionally, it produces `merfish_linear_gradient_detail.csv` with a linear pattern GLM, `merfish_quadratic_curvature_detail.csv` with a curved (quadratic) pattern GLM, and `merfish_bulk_concordance_summary.csv` summarising the agreement between categories of genes present in both MERFISH and bulk sequencing.
3. **`marker_gene_expression_fig7f.py`** - generates the renderings of all genes with soma coloured based on the gene's expression levels. Generates all categorised genes; meanwhile, the figure uses a few examples.
4. **`MERFISH_subtype_plotting_fig7h.py`** - generates the rendering of the cells in the MERFISH dataset, coloured by their super-type.
5. **`discrete_vs_continuous_model.ipynb`** - provides analysis of continuous model selection over discrete. It generates files `continuous_vs_discrete.csv` (gene categorisation by discrete, using EC compartments, versus continuous, using position on mediolateral (z) axis models), `continuous_vs_discrete_cross-validation.csv` (cross-validates the gene categories by training the model on the subset of cells and testing on the remaining ones), `continuous_vs_discrete_dorsoventral_covariate.csv` (inspects if the model assigns different categories to genes if the dorsoventral axis position is used as a covariate).

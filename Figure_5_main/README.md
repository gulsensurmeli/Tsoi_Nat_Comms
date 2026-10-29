# Figure 5 plots

This folder contains all the files needed for data processing as well as the output files and plots for Figure 5.

To set up the environment with appropriate packages to run these scripts, set up the working directory in this folder and run:

```
conda env create -f set_up.yml
conda activate ec_analysis_fig_5
```

## Single neuron reconstruction (SNR) files

The raw registered neuronal tracing files (swc format) are in a folder named **raw_snr_data**. These are the starting files for data preprocessing. These files are used to extract information about each neuron, such as axon length and axonal endpoint count.

**Ntype in files identifies the structural identity of that point within the neurite:**

| Ntype | Structure |
|---|---|
| 0 | undefined |
| 1 | soma |
| 2 | axon |
| 3 | (basal) dendrite |
| 4 | apical dendrite |
| 5 | axon fork point |
| 6 | axon end point |
| 7 | (basal) dendrite end point |
| 8 | apical dendrite end point |
| 11 | (basal) dendrite fork point |
| 12 | apical dendrite fork point |

## Data processing and plotting pipeline

The data is first processed using the **`preprocessing_of_files.ipynb`** script. It normalises the neurons to the same hemisphere and generates the following Excel files:

1. **`mouseNR_neuronID_neurites.xlsx`** -> has all the points for neurites, their branching data, and endpoints. It uses the atlas volume (`annotation_25.nrrd`) and atlas structure tree (`structure_tree_safe.csv`) to assign area names and their paths to each node of the reconstruction.
2. **`mouseNR_neuronID_length.xlsx`** -> length of axons in all areas as registered in the atlas (in this and previous files ID/id stands for the registered area ID). The files with "ipsi" in their names indicate that the neuron projects to both hemispheres, but the length counts are restricted to the ipsilateral hemisphere.

The processing workflow also pools data from all neurons into the master file (**`master_detailed_comment.csv`**), which contains each neuron ID, axonal length and endpoint information in selected areas, soma coordinates, projection group, and comments if the neuron has been removed from analysis.

The outputs are saved into `preprocessed_data`.


The plotting goes as follows:

**`Soma_plots.py`**, **`Density_maps.py`**, and **`Quantitative_plots_fig5e-h.ipynb`** utilise the processed data and can be run in any order.
   - **`Soma_plots_fig5b-c.py`** generates soma maps labelled based on their projection group with and without group centroids.
   - **`Density_maps_fig5d.py`** uses the Brainrender density function to create axon density maps based on axonal points.
   - **`Quantitative_plots_fig5e-h.ipynb`** generates the quantitative analysis plots, such as heatmaps and violin plots for axon length and endpoints for each group in the isocortex and other areas.
   
Figure 5a plots and the generation script can be found in the Sup Figure 6 folder.


`MI_bootstrap_analysis.Rmd` contains the analysis run for soma positioning analysis.

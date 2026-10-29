
##---------- Load all required packages --------------------

suppressPackageStartupMessages({

  library(mvtnorm)
  library(extraDistr)
  library(SciViews)
  library(doParallel)
  library(foreach)
  library(ggplot2)
  library(tidyverse)
  library(ggpubr)
  library(fields)
  library(Matrix)
  library(Rtsne)
  library(ggpubr)
  library(ggrepel)
  library(parallel)
  library(progress)
  library(pheatmap)
  library(RColorBrewer)
  library(readxl)
  library(truncnorm)
  library(clevr)
  library(pbapply)
  library(grid)
  library(gridExtra)
  library(coda)
  library(mcclust)
  library(ClusterR)
  library(latex2exp)
  library(pracma)
  library(cascsim)
})

`%notin%` <- Negate(`%in%`)


# MCMC
source('~/old_hbmap/R/mcmc_steps.R')
# source('R/mcmc_run_all.R')
source("~/old_hbmap/R/mcmc_run_all_prior2alpha.R")
source('~/old_hbmap/R/mcmc_run_post.R')
source('~/old_hbmap/R/mcmc_run_omega_JM.R')
source('~/old_hbmap/R/cluster_label.R')

# Posterior similarity matrix and plot
source('~/old_hbmap/R/similarity_matrix.R')
source('~/old_hbmap/R/posterior_similarity_plot.R')

# Optimal clustering
source('~/old_hbmap/R/optimal_clustering.R')
source('~/old_hbmap/R/clustering_estimate.R')
source('~/old_hbmap/R/clustering_estimate_plot.R')
source('~/old_hbmap/R/difference_in_omega_jm.R')
source('~/old_hbmap/R/difference_in_omega_jm_plot.R')

# Cluster labeling
source('~/old_hbmap/R/cluster_label.R')
source('~/old_hbmap/R/projection_by_EC.R')

# Projection probability by cluster
source('~/old_hbmap/R/pp_by_z.R')

# Posterior predictive checks
source('~/old_hbmap/R/ppc.R')
source('~/old_hbmap/R/ppc_single.R')

# k means and k medoids model
#source('~/old_hbmap/R/kmeans_model.R')

#source('~/old_hbamp/R/binom_cluster_reorder.R')
#source('~/old_hbmap/R/find_large_cluster.R')

#source('~/old_hbmap/R/projection_by_EVI.R')
#source('~/old_hbamp/R/binomial_model.R')

#source('jinlu/added_noise.R')
#source('jinlu/gel_plot.R')
#source('jinlu/pp_by_z_2.R')
#source('jinlu/projection_probability_code.R')
#source('R/weights_analysis.R')
#source("R/projection_prob_functions.R")
#source("jinlu/projection_probability_code.R")

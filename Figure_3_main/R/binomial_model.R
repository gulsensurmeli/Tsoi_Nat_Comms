

#' Binomial model 1 - assume the probability of projection is the same for all regions.
#'
#' @param data The input dataset of length M, where M is the total number of datasets. Each item has R rows and C columns, where R is the total number of variables and C is the total
#' number of observations.
#' @param threshold Only projections greater or equal to the defined threshold are considered as true projections. The method is used to remove potential noise
#' during the data collection process.
#'
#' @return The function returns a list containing: the allocation of neurons (allocation),
#' the projection regions of each cluster (cluster_label),
#' and the summary of each cluster, including the p-value of binomial test, the expected and observed counts of each cluster, and classify clusters into over-represented and under-represented clusters,
#' and significant and insignificant clusters (cluster_summary).
#' @export

binomial_model <- function(data,
                           threshold = 1){


  # If the data is supplied in a list
  if(is.list(data)){

    C <- sapply(1:length(data), function(m) ncol(data[[m]]))
    M <- length(data)
    data <- do.call(cbind, data)

  }else{

    C <- ncol(data)
    M <- 1
  }

  R <- nrow(data)



  # Calculate Nf, N_0 and p_ea
  Nf <- ncol(data)

  si <- sapply(1:nrow(data),
               function(r) length(which(data[r,] >= threshold)))


  N_0 <- round(fzero(function(x) Nf - x*(1-prod(1-si/x)), Nf)$x, 0)

  p_e <- fzero(function(x) Nf - N_0*(1-(1-x)^nrow(data)), 0.5)$x

  # Compute the expected count for different number of projecting motifs
  expected_number <- sapply(1:nrow(data),
                            function(r){

                              N_0*p_e^(r)*(1-p_e)^(R-r)
                            })



  # Find the set of unique clusters
  cluster_neuron <- sapply(1:ncol(data),
                           function(c) paste(rownames(data)[data[,c] >= threshold],
                                             collapse = " "))

  # Unique set of cluster
  cluster_label <- unique(cluster_neuron)

  # Region names
  region_names <- rownames(data)

  # Neuron allocation
  allocation <- sapply(1:ncol(data),
                       function(c) which(cluster_label == cluster_neuron[c]))


  # Separate allocations
  C_cumsum <- c(0, cumsum(C))

  allocation <- lapply(1:M,
                       function(m) allocation[(C_cumsum[m]+1):(C_cumsum[m+1])])


  cluster_summary <- data.frame(cluster_name = cluster_label,
                                N_projecting_regions = sapply(1:length(cluster_label),
                                                              function(i) length(str_split(cluster_label[[i]], pattern = ' ')[[1]]))
  )

  # Expected number of projections
  cluster_summary$size_expect <- expected_number[cluster_summary$N_projecting_regions]

  # Observed number of projections
  cluster_summary$size_observed <- sapply(1:nrow(cluster_summary),
                                          function(j) length(which(cluster_neuron == cluster_label[j])))

  cluster_summary <- cluster_summary %>%
    mutate(cluster.type = ifelse(size_expect < size_observed, 'over-represented', 'under-represented'))

  # p-value
  cluster_summary$p_value <- sapply(1:nrow(cluster_summary),
                                    function(j) binom.test(x = cluster_summary$size_observed[j],
                                                           n = ncol(data),
                                                           p = (p_e)^cluster_summary$N_projecting_regions[j]*(1-p_e)^(R-cluster_summary$N_projecting_regions[j]),
                                                           alternative = 'two.sided')$p.value)


  # Threshold
  threshold_p_value <- -log(0.05/nrow(cluster_summary), base = 10)

  cluster_summary$significant <- ifelse(-log(cluster_summary$p_value, base = 10) < threshold_p_value,
                                        'insignificant',
                                        'significant')



  return(list('allocation' = allocation,
              'cluster_label' = cluster_label,
              'cluster_summary' = cluster_summary))
}



#' @export
binomial_model_clustering_only <- function(data,
                                           threshold = 1){
  
  
  # If the data is supplied in a list
  if(is.list(data)){
    
    C <- sapply(1:length(data), function(m) ncol(data[[m]]))
    M <- length(data)
    data <- do.call(cbind, data)
    
  }else{
    
    C <- ncol(data)
    M <- 1
  }
  
  R <- nrow(data)
  
  
  # Find the set of unique clusters
  cluster_neuron <- sapply(1:ncol(data),
                           function(c) paste(rownames(data)[data[,c] >= threshold],
                                             collapse = " "))
  
  # Unique set of cluster
  cluster_label <- unique(cluster_neuron)
  
  # Region names
  region_names <- rownames(data)
  
  # Neuron allocation
  allocation <- sapply(1:ncol(data),
                       function(c) which(cluster_label == cluster_neuron[c]))
  
  
  # Separate allocations
  C_cumsum <- c(0, cumsum(C))
  
  allocation <- lapply(1:M,
                       function(m) allocation[(C_cumsum[m]+1):(C_cumsum[m+1])])
  
  
  cluster_summary <- data.frame(cluster_name = cluster_label,
                                N_projecting_regions = sapply(1:length(cluster_label),
                                                              function(i) length(str_split(cluster_label[[i]], pattern = ' ')[[1]]))
  )
  
  
  # Observed number of projections
  cluster_summary$size_observed <- sapply(1:nrow(cluster_summary),
                                          function(j) length(which(cluster_neuron == cluster_label[j])))
  
  
  return(list('allocation' = allocation,
              'cluster_label' = cluster_label,
              'cluster_summary' = cluster_summary))
}


#' Binomial model 2 - assume the probability of projection of each region is independent.
#'
#' @param data The input dataset of length M, where M is the total number of datasets. Each item has R rows and C columns, where R is the total number of variables and C is the total
#' number of observations.
#' @param threshold Only projections greater or equal to the defined threshold are considered as true projections. The method is used to remove potential noise
#' during the data collection process.
#'
#' @return The function returns a list containing: the allocation of neurons (allocation),
#' the projection regions of each cluster (cluster_label),
#' and the summary of each cluster, including the p-value of binomial test, the expected and observed counts of each cluster, and classify clusters into over-represented and under-represented clusters,
#' and significant and insignificant clusters (cluster_summary),
#' and the plot to show relationship between p-values and log ratio between the observed and expected counts (plot_output).
#' @export

binomial_model2 <- function(data,
                            threshold = 1){


  # If the data is supplied in a list
  if(is.list(data)){

    M <- length(data)
    C <- sapply(1:length(data), function(m) ncol(data[[m]]))
    data <- do.call(cbind, data)

  }else{

    M <- 1
    C <- ncol(data)
  }

  R <- nrow(data)



  # Calculate Nf, N_0 and p_ea
  Nf <- ncol(data)

  si <- sapply(1:nrow(data),
               function(r) length(which(data[r,] >= threshold)))

  N_0 <- round(fzero(function(x) Nf - x*(1-prod(1-si/x)), Nf)$x, 0)


  p_i <- si/N_0



  # Find the set of unique clusters
  cluster_neuron <- sapply(1:ncol(data),
                           function(c) paste(rownames(data)[data[,c] >= threshold],
                                             collapse = " "))

  # Unique set of cluster
  cluster_label <- unique(cluster_neuron)

  # Region names
  region_names <- rownames(data)

  # Neuron allocation
  allocation <- sapply(1:ncol(data),
                       function(c) which(cluster_label == cluster_neuron[c]))


  # Separate allocations
  C_cumsum <- c(0, cumsum(C))

  allocation <- lapply(1:M,
                       function(m) allocation[(C_cumsum[m]+1):(C_cumsum[m+1])])


  cluster_summary <- data.frame(cluster_name = cluster_label,
                                N_projecting_regions = sapply(1:length(cluster_label),
                                                              function(i) length(str_split(cluster_label[[i]], pattern = ' ')[[1]]))
  )



  # Expected number of projections
  cluster_summary$size_expect <- sapply(1:nrow(cluster_summary),
                                        function(j){

                                          prod(c(p_i[which(rownames(data) %in%  str_split(cluster_summary$cluster_name[[j]], pattern = ' ')[[1]])],
                                                 1-p_i[which(rownames(data) %notin%  str_split(cluster_summary$cluster_name[[j]], pattern = ' ')[[1]])]))*N_0
                                        }
  )


  # Observed number of projections
  cluster_summary$size_observed <- sapply(1:nrow(cluster_summary),
                                          function(j) length(which(cluster_neuron == cluster_label[j])))

  cluster_summary <- cluster_summary %>%
    mutate(cluster.type = ifelse(size_expect < size_observed, 'over-represented', 'under-represented'))


  # p-value
  cluster_summary$p_value <- sapply(1:nrow(cluster_summary),
                                    function(j) binom.test(x = cluster_summary$size_observed[j],
                                                           n = N_0,
                                                           p = prod(c(p_i[which(rownames(data) %in%  str_split(cluster_summary$cluster_name[[j]], pattern = ' ')[[1]])],
                                                                      1-p_i[which(rownames(data) %notin%  str_split(cluster_summary$cluster_name[[j]], pattern = ' ')[[1]])])),
                                                           alternative = 'two.sided')$p.value)


  threshold_p_value <- -log(0.05/(M*length(which(cluster_summary$N_projecting_regions > 1))), base = 10)

  plot_output <- cluster_summary %>%
    filter(N_projecting_regions > 1) %>%
    ggplot(mapping = aes(x = log(size_observed/size_expect, base = 2),
                         y = -log(p_value, base = 10)))+
    geom_point()+
    geom_hline(yintercept = threshold_p_value, linetype = 'dashed', color = 'grey')+
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey')+
    theme_bw()


  # Significance
  cluster_summary$significant <- ifelse(-log(cluster_summary$p_value, base = 10) < threshold_p_value,
                                        'insignificant',
                                        'significant')


  return(list('allocation' = allocation,
              'cluster_label' = cluster_label,
              'cluster_summary' = cluster_summary,
              'plot_output' = plot_output))
}




#' Reordering of binomial clusters depending on the region of strongest projection of each cluster.
#' @param Y The input dataset of length M, where M is the total number of datasets. Each item has R rows and C columns, where R is the total number of variables and C is the total
#' number of observations.
#' @param binomial_output Direct output from function binomial_model() or binomial_model2().
#' @return A list containing the following items: neuron allocations after reordering (Z_new),
#' cluster labels after reordering (cluster_label),
#' and clustering summary after reordering (cluster_summary).
#' @export

binom_cluster_reorder <- function(Y,
                                  binomial_output){

  Z <- binomial_output$allocation

  Z_unlist <- unlist(Z)
  Y_join <- do.call(cbind, Y)

  M <- length(Y)
  C <- sapply(1:M,
              function(m) ncol(Y[[m]]))

  # Calculate proportions
  Y_proportion <- Y_join
  N_data <- colSums(Y_join)

  for(i in 1:ncol(Y_join)){

    Y_proportion[,i] <- Y_proportion[,i]/N_data[i]
  }

  # Calculate the average projection strength of neurons in each cluster
  average_q <- lapply(sort(unique(Z_unlist)),
                      function(j){

                        data_j <- matrix(Y_proportion[,which(Z_unlist == j)],

                                         ncol = length(which(Z_unlist==j)))

                        matrix(rowMeans(data_j),
                               ncol = 1)
                      })

  average_q <- do.call(cbind, average_q)

  # Reorder clusters
  strong.proj <- sapply(1:ncol(average_q),
                        function(j) which(average_q[,j] == max(average_q[,j]))[1])



  df0 <- lapply(1:nrow(average_q),
                function(i){

                  # Clusters with the strongest projection to region i
                  strong.proj.i <- which(strong.proj==i)

                  if(length(strong.proj.i) != 0){

                    qj_estimate_i <- matrix(average_q[,strong.proj.i],
                                            nrow = length(strong.proj.i))

                    output <- strong.proj.i[order(-qj_estimate_i[,i])]
                  }else{

                    output <- NULL
                  }

                }
  )

  old_ordering <- unlist(df0)

  Z_new <- lapply(1:M,
                  function(d){

                    sapply(1:C[d],
                           function(c) which(old_ordering == Z[[d]][c]))
                  })


  return(list('allocation' = Z_new,
              'cluster_label' = binomial_output$cluster_label[old_ordering],
              'cluster_summary' = binomial_output$cluster_summary[old_ordering,]))
}



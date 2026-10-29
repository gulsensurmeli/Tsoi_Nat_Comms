
find_large_cluster <- function(mcmc_run_omega_output,
                               tau = 0.01,
                               top.n){
  
  # Trace of omega
  omega.trace <- mcmc_run_omega_output$omega_J_M_output
  
  # Number of clusters
  J <- nrow(omega.trace[[1]])
  
  # Posterior probability of having w_j > tau
  posterior.probability <- sapply(1:J,
                                  function(j){
                                    
                                    omega.trace.j <- sapply(1:length(omega.trace),
                                                            function(t) omega.trace[[t]][j])
                                    
                                    mean(omega.trace.j > tau)
                                  })
  
  return(which(rank(desc(posterior.probability)) <= top.n))
}


opt.clustering.frequency1_large <- function(clustering,
                                            main='',
                                            EC_label,
                                            large_cluster_index){
  
  
  cluster.names <- paste('cluster', 
                         1:length(unique(unlist(clustering))))
  
  
  # Table of LEC and MEC
  LEC_df <- data.frame(cluster = factor(cluster.names,
                                        levels = cluster.names),
                       EC_label = 'LEC',
                       counts = as.vector(table(factor(unlist(clustering)[which(unlist(EC_label) == 'LEC')],
                                                       levels = 1:length(cluster.names)))))
  
  MEC_df <- data.frame(cluster = factor(cluster.names,
                                        levels = cluster.names),
                       EC_label = 'MEC',
                       counts = as.vector(table(factor(unlist(clustering)[which(unlist(EC_label) == 'MEC')],
                                                       levels = 1:length(cluster.names)))))
  
  df <- rbind(LEC_df, MEC_df) %>%
    filter(cluster %in% paste('cluster', large_cluster_index))
  
  # Plot to show proportion
  
  ggplot(df, mapping = aes(x = cluster, y = counts, fill = EC_label))+
    geom_bar(stat="identity")+
    theme_bw()+
    ylab('Number of LEC and MEC neurons')+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
}


opt.clustering.frequency_large <- function(clustering,
                                           main = '',
                                           large_cluster_index){
  
  cluster.names <- paste('cluster', 1:length(unique(unlist(clustering))))
  
  loop.result <- lapply(1:length(clustering), function(m){
    
    data.frame(cluster = factor(cluster.names,
                                levels = cluster.names),
               mouse = paste('mouse', m),
               counts = as.vector(table(factor(clustering[[m]],
                                               levels = 1:length(cluster.names)))))
  })
  
  # Frequency table
  z.frequency <- do.call(rbind, loop.result) %>%
    filter(cluster %in% paste('cluster', large_cluster_index))
  
  
  # Bar chart
  z.frequency %>%
    ggplot(mapping = aes(x = cluster, y = counts, fill = mouse))+
    geom_bar(stat="identity")+
    theme_bw()+
    ylab('Number of neurons')+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(main)
}


opt.clustering.frequency2_large <- function(clustering,
                                            main='',
                                            EC_label,
                                            large_cluster_index){
                                      
  
  
  cluster.names <- paste('cluster', 
                         1:length(unique(unlist(clustering))))
  
  
  # Table of LEC and MEC
  LEC_df <- data.frame(cluster = factor(cluster.names,
                                        levels = cluster.names),
                       EC_label = 'LEC',
                       counts = as.vector(table(factor(unlist(clustering)[which(unlist(EC_label) == 'LEC')],
                                                       levels = 1:length(cluster.names)))))
  
  MEC_df <- data.frame(cluster = factor(cluster.names,
                                        levels = cluster.names),
                       EC_label = 'MEC',
                       counts = as.vector(table(factor(unlist(clustering)[which(unlist(EC_label) == 'MEC')],
                                                       levels = 1:length(cluster.names)))))
  
  df <- rbind(LEC_df, MEC_df) %>%
    filter(cluster %in% paste('cluster', large_cluster_index))
  
  df2 <- df %>%
    group_by(cluster) %>%
    summarise(sum_counts = sum(counts)) 
  
  
  
  # Plot to show proportion
  ggplot(df, 
         aes(x = cluster, y = counts, fill = EC_label)) +
    geom_col(position = "fill")+
    theme_bw()+
    ylab('Proportion of LEC and MEC neurons')+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    annotate(geom = 'text', 
             label = df2$sum_counts,
             x = (1:nrow(df2)), y = 0.85, hjust = 0, vjust = 1, size = 5,
             angle = 90)+
    geom_hline(yintercept = mean(unlist(EC_label) == 'MEC'), col = "red")
  
}

##-- Estimated projection strengths of large clusters
estimated_projection_strength_large <- function(mcmc_run_unique_output,
                                                large_cluster_index,
                                                region.name){
  
  estimated.projection.df <- mcmc_run_unique_output$estimated.projection.df %>%
    filter(cluster %in% paste('cluster', large_cluster_index))
  
  ggplot2::ggplot(estimated.projection.df,
                  mapping = aes(x = factor(region, levels = region.name),
                                y = projection.med,
                                group = cluster,
                                color = class)) +
    geom_line()+
    geom_point()+
    geom_errorbar(aes(ymin = projection.lower,
                      ymax = projection.upper),
                  width = 0.1)+
    theme_bw()+
    ylim(c(0,1))+
    ylab('projection strength')+
    xlab('region')+
    facet_wrap(~cluster)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=10),
          plot.title = element_text(size=10))
}


q_tilde_plot_large <- function(mcmc_run_unique_output,
                               large_cluster_index,
                               region.name){
  
  q_tilde_001 <- mcmc_run_unique_output$q_tilde_001
  J <- nrow(q_tilde_001)
  
  
  q_tilde_001_table <- data.frame(region = rep(region.name,
                                               each = J),
                                  cluster = rep(paste('cluster', 1:J),
                                                R),
                                  probability = as.vector(q_tilde_001))
  
  q_tilde_001_table$region <- factor(q_tilde_001_table$region,
                                     levels = region.name)
  
  q_tilde_001_table$cluster <- factor(q_tilde_001_table$cluster,
                                      levels = paste('cluster', 1:J))
  
  q_tilde_001_table %>%
    filter(cluster %in% paste('cluster', large_cluster_index)) %>%
    ggplot(mapping = aes(x = cluster,
                         y = region,
                         fill = probability))+
    geom_tile()+
    theme_bw()+
    scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                         values = scales::rescale(c(0, 0.25, 1)),
                         limits = c(0,1))+
    xlab('cluster')+
    ylab('regions')+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
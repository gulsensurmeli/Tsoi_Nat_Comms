# Plot of line graph to show projection strengths of neurons in each cluster, color-coded by the injection site

projection_by_EC <- function(Y,
                             EC_label,
                             Z,
                             region_name){
  
  
  M <- length(Y)
  
  # Projection strength of each neuron
  Y_prop <- lapply(1:M,
                   function(m){
                     
                     matrix(Y[[m]]/rep(colSums(Y[[m]]), each = nrow(Y[[m]])),
                            
                            nrow = nrow(Y[[m]]),
                            ncol = ncol(Y[[m]]))
                   })
  
  # Data frame
  list0 <- lapply(1:M,
                  function(m){
                    
                    data.frame(projection_strength = as.vector(Y_prop[[m]]),
                               region_name = rep(region_name, ncol(Y[[m]])),
                               cluster = rep(Z[[m]], each = nrow(Y[[m]])),
                               EC_label = rep(EC_label[[m]], each = nrow(Y[[m]])))
                  })
  
  df <- do.call(rbind, list0)
  
  # Change neuron into factor
  df$neuron <- factor(rep(1:length(unlist(Z)), each = nrow(Y[[1]])),
                      levels = rep(1:length(unlist(Z))))
  
  # Change cluster into factors
  df$cluster <- factor(paste('cluster', df$cluster),
                       levels = paste('cluster', 1:max(unlist(Z))))
  
  # Change region names into factors
  df$region_name <- factor(df$region_name,
                           levels = region_name)
  
  # Line graph
  ggplot(df)+
    geom_line(mapping = aes(x = region_name,
                            y = projection_strength,
                            colour = EC_label,
                            group = interaction(neuron, EC_label)))+
    facet_wrap(~cluster)+
    theme_bw()+
    xlab('region')+
    ylab('projection strengths')+ 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
}
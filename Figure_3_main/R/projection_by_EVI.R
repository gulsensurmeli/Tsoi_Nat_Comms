# Create a function that computes the contribution of each data point to the
# expected VI

evi.contribution = function(Zmat, Zhat){
  n = length(Zhat)
  S = dim(Zmat)[1]
  evi.i = function(i,Zmat, Zhat){
    c = 1/n*log2(sum(Zhat==Zhat[i])/n) + sum(apply(Zmat,1,function(x){1/n*log2(sum(x==x[i])/n)}))/S -2/n*sum(apply(Zmat,1,function(x,y){log2(sum((x==x[i])&(y==y[i]))/n)},y=Zhat))/S
    return(c)
  }
  contrib = unlist(lapply(c(1:n), evi.i, Zmat=Zmat, Zhat=Zhat))
  return(contrib)
}


projection_by_evi <- function(Y,
                              evi,
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
                               evi_label = rep(evi[[m]], each = nrow(Y[[m]])))
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
                            colour = evi_label,
                            group = interaction(neuron, evi_label)))+
    facet_wrap(~cluster)+
    scale_color_gradient2(mid='blue',high='red') +
    theme_bw()+
    xlab('region')+
    ylab('projection strengths') +
    labs(colour = "EVI")
  
}


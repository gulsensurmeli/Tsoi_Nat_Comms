#' Computation of the probability that the projection strength of neurons
#' in each cluster is greater than a small threshold for each brain region

q_tilde <- function(neuron_projection_df,
                    epsilon = 0.01){

  # Trace of q
  q_trace <- neuron_projection_df$q_star_1_J_output
  gamma_trace <- neuron_projection_df$gamma_star_1_J_output

  # Estimated allocation
  Z <- neuron_projection_df$Z

  # Trace length
  trace.length <- length(q_trace)

  df0 <- lapply(1:nrow(q_trace[[1]]),
                function(j){

                  df0_j <- lapply(1:trace.length,
                                  function(i){

                                    vec0 <- 1-pbeta(q = epsilon,
                                                    shape1 = q_trace[[i]][j,]*gamma_trace[[i]][j],
                                                    shape2 = (1-q_trace[[i]][j,])*gamma_trace[[i]][j],
                                                    lower.tail = TRUE)

                                    matrix(vec0,
                                           nrow = 1)
                                  })

                  df0_j <- do.call(rbind,
                                   df0_j)

                  matrix(round(colMeans(df0_j),3),
                         nrow = 1)

                })

  # Dimension: J x R
  df0 <- do.call(rbind,
                 df0)

  return(df0)
}

q_tilde2 <- function(neuron_projection_df, 
                     n = 100){
  
  # Trace of q
  q_trace <- neuron_projection_df$q_star_1_J_output
  gamma_trace <- neuron_projection_df$gamma_star_1_J_output
  
  # Estimated allocation
  Z <- neuron_projection_df$Z
  R <- ncol(q_trace[[1]])
  
  # Trace length
  trace.length <- length(q_trace)
  
  df0 <- lapply(1:nrow(q_trace[[1]]),
                function(j){
                  
                  df0_j <- lapply(1:trace.length,
                                  function(i){
                                    alphaj = q_trace[[i]][j,]*gamma_trace[[i]][j]
                                    vec0 <- sapply(c(1:R), function(r){1- ddirmultinomial(c(0,n),n,c(alphaj[r],sum(alphaj)-alphaj[r]))})
                                    
                                    matrix(vec0,
                                           nrow = 1)
                                  })
                  
                  df0_j <- do.call(rbind,
                                   df0_j)
                  
                  matrix(colMeans(df0_j),
                         nrow = 1)
                  
                })
  
  # Dimension: J x R
  df0 <- do.call(rbind,
                 df0)
  
  return(df0)
}


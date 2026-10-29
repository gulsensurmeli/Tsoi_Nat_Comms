#' Synthetic data generation

synthetic_data_generation <- function(mcmc_run_all_output,
                      Y,
                      regions.name = NULL,
                      M,
                      C,
                      noise.levels = NULL){
  
  
  
  length_per_chain <- length(mcmc_run_all_output$Z_output)
  
  C_data <- mcmc_run_all_output$C
  R <- mcmc_run_all_output$R
  M_data <- mcmc_run_all_output$M
  
  if(is.null(regions.name)){
    
    regions.name <- paste('region', 1:R)
  }
  
  
  # Sum of counts for each neuron cell
  N_CM <- lapply(1:M,
                 function(m) colSums(Y[[m]]))
  
  # theta
  target.index <- sample(length_per_chain,
                         size = 1,
                         replace = FALSE)
  
  
  theta <- list(q_j_star = mcmc_run_all_output$q_star_1_J_output[[target.index]],
                gamma_j_star = mcmc_run_all_output$gamma_star_1_J_output[[target.index]],
                omega = mcmc_run_all_output$omega_output[[target.index]],
                alpha = mcmc_run_all_output$alpha_output[target.index],
                Z = mcmc_run_all_output$Z[[target.index]])
  
  # Distribution for counts
  mu = max(mean(unlist(N_CM)-5),1)
  disp = mu2^2/max(var(unlist(N_CM)-5)-mu2,1)
  # Cluster-specific paramters
  mu_j_star = sapply(1:J, function(j){
    if(sum(unlist(theta$Z)==j)<1){
      mu_j = mu
    }else{
      mu_j = max(mean(unlist(N_CM)[unlist(theta$Z)==j]-5),1)
    }
    mu_j
  })
  disp_j_star = sapply(1:J, function(j){
    if(sum(unlist(theta$Z)==j)<5){
      disp_j = disp
    }else{
      disp_j = mu_j_star[j]^2/max(var(unlist(N_CM)[unlist(theta$Z)==j]-5)-mu_j_star[j], 1)
    }
  })
  
  # Simulate synthetic data
  Y_rep = list()
  Z_rep = list()
  omega_MJ_rep = list()
  
  for(m in 1:M){
    omega_m = rdirichlet(1, theta$omega*theta$alpha)
    omega_MJ_rep[[m]] = omega_m
    
    Z_rep_m <- apply(rmultinom(C[m],1,omega_m),2,function(x){which(x==1)})
    Z_rep[[m]] = Z_rep_m
    
    N_rep_m = rnbinom(C[m], mu = mu_j_star[Z_rep_m], size = disp_j_star[Z_rep_m]) +5
    
    Y_rep_m <- do.call(cbind, lapply(1:C[m],
                                     function(c){
                                       
                                       
                                       prob.c <- rdirichlet(n = 1,
                                                            theta$q_j_star[Z_rep_m[c],]*theta$gamma_j_star[Z_rep_m[c]])
                                       
                                       rmultinom(n = 1, size = N_rep_m[c], prob = prob.c)
                                     }))
    
    rownames(Y_rep_m) <- regions.name
    Y_rep[[m]] = Y_rep_m
  }
  
  # Add some noise to the Y_rep
  noisy_data_list = NULL
  if(length(noise.levels)>0){
    noisy_data_list <- lapply(noise.levels,
                              function(noise.level){
                                
                                lapply(1:M,
                                       function(m){
                                         
                                         noisy_data_m <- lapply(1:R,
                                                                function(r){
                                                                  
                                                                  noisy_data_m0 <- lapply(1:ncol(Y_rep[[m]]),
                                                                                          function(i){
                                                                                            
                                                                                            prob_mi <- rep(0, R)
                                                                                            prob_mi[r] <- (1-noise.level)
                                                                                            prob_mi[-r] <- noise.level/(R-1)
                                                                                            
                                                                                            matrix(rmultinom(n = 1, size = Y_rep[[m]][r,i], prob = prob_mi),
                                                                                                   ncol = 1)
                                                                                          })
                                                                  
                                                                  do.call(cbind, noisy_data_m0)
                                                                })
                                         
                                         # Add up and round to integer
                                         noisy_data_m <- Reduce('+', noisy_data_m)
                                         rownames(noisy_data_m) <- regions.name
                                         
                                         noisy_data_m
                                       })
                              })
  }
  
  
  return(list('synthetic_data' = Y_rep,
              'noisy_synthetic_data' = noisy_data_list,
              'Z_synthetic' = Z_rep,
              'theta' = theta))
  
  
  
}

# Compare the variance of the brain-specific weights across mice with a null model,
# where the null model assumes equal weight across mice

# Save to server hard drive
# save(mcmc_all_hans,
#      mcmc_unique_hans,
#      omega_JM_mcmc,
#      file = "C:/Users/s1847159/OneDrive - University of Edinburgh/BrainConnectivity - WADE Sara's files/Results/local_weights_input.RData")

# local_weights_analysis_vc <- local_weights_analysis(N = 100,
#                                                     mcmc_run_all_output = mcmc_all_hans,
#                                                     mcmc_run_post_output = mcmc_unique_hans,
#                                                     mcmc_run_omega_output = omega_JM_mcmc)

local_weights_analysis <- function(N = 100,
                                   mcmc_run_all_output,
                                   mcmc_run_post_output,
                                   mcmc_run_omega_output,
                                   mouse_indices=NULL){
  
  J <- max(unlist(mcmc_run_post_output$Z))
  # Proportion of neurons in each cluster across all mice
  if(is.null(mouse_indices)){mouse_indices = c(1:length(mcmc_run_post_output$Z))}
  Z <- mcmc_run_post_output$Z[mouse_indices]
  M <- length(mouse_indices)
  
  # w <- lapply(1:M, function(m){
  #   unlist(lapply(1:J, function(j){
  #     sum(Z[[m]]==j)/length(Z[[m]])
  #   }))
  # })
  # w = matrix(unlist(w), J, M)
  # w = rowMeans(w)
  
  # w <- sapply(1:J, function(j) length(which(unlist(Z)==j)))/length(unlist(Z))
  
  w =  mcmc_run_omega_output$mean_omega_J
  
  # Empirical variance
  w_jm_empirical <- mcmc_run_omega_output$omega_J_M_output
  
  w_jm_variance_empirical <- lapply(1:length(w_jm_empirical),
                                    function(t) matrix(apply(w_jm_empirical[[t]][,mouse_indices], 1, var),
                                                       nrow = 1))
  
  
  w_jm_variance_empirical <- do.call(rbind, w_jm_variance_empirical)
  
  w_jm_variance_empirical <- colMeans(w_jm_variance_empirical)
  
  
  # Simulate a number of allocations under null
  Z_null <- lapply(1:N,
                   function(n){
                     
                     # Allocation
                     #J_n <- 0
                     
                     #while(J_n < J){
                       
                      Z_n <- lapply(1:M,
                                     function(m) rcat(n = length(Z[[m]]),
                                                      prob = w))
                       
                       J_n <- length(unique(unlist(Z_n)))
                       
                       empty.ind = sapply(1:J, function(j) length(which(unlist(Z_n)==j)))==0
                       
                     #}
                     
                     
                     mcmc_run_post_n <- mcmc_run_post_output
                     mcmc_run_post_n$Z <- Z_n
                     
                     # Obtain posterior samples of w_jm
                     w_sampling <- mcmc_run_omega_JM(mcmc_run_all_output = mcmc_run_all_output,
                                                     mcmc_run_post_output = mcmc_run_post_n,
                                                     J = J,
                                                     thinning = 5,
                                                     burn_in = 5000,
                                                     number_iter = 15000)
                     
                     w_jm_samples <- w_sampling$omega_J_M_output
                     
                     # Compute the variance of each iteration
                     w_jm_variance <- lapply(1:length(w_jm_samples),
                                             function(t) matrix(apply(w_jm_samples[[t]], 1, var),
                                                                nrow = 1))
                     
                     w_jm_variance <- do.call(rbind, w_jm_variance)
                     
                     
                     # Return posterior mean
                     matrix(colMeans(w_jm_variance), nrow = 1)
                     
                   })
  
  # Combine
  Z_null <- do.call(rbind, Z_null)
  
  # Proportion of times the variance is larger than the null
  proportion_large <- sapply(1:J, 
                             function(j) mean(w_jm_variance_empirical[j] > Z_null[,j]))
  
  return(data.frame('cluster' = 1:J,
                    'probability' = proportion_large))
  
}
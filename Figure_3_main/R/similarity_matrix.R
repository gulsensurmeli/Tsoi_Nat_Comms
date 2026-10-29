
#' Posterior similarity matrix

similarity_matrix <- function(mcmc_run_all_output){
  
  # Input from MCMC
  M <- mcmc_run_all_output$M
  C <- mcmc_run_all_output$C
  
  C_cumsum = c(0,cumsum(C))
  
  # Put z in a matrix
  Zmat = matrix(unlist(mcmc_run_all_output$Z_output), length(mcmc_run_all_output$Z_output), sum(C),byrow = TRUE)
  
  # Posterior similarity between all observations
  psm <- comp.psm(Zmat)
  
  # Posterior similarity between observations across 2 datasets
  psm_within <- list(M)
  
  for(m1 in c(1:M)){
    
    psm_within[[m1]] = list(M)
    
    for(m2 in c(1:M)){
      
      psm_within[[m1]][[m2]] = psm[(C_cumsum[m1]+1):C_cumsum[m1+1],(C_cumsum[m2]+1):C_cumsum[m2+1]]
      
    }
  }
  
  ##-------------------------------- Return both -----------------------------------------
  
  return(list('psm.within' = psm_within,
              'psm.combined' = psm))
  
}

similarity_matrix_mini <- function(allocation_trace){
  
  loop.result <- lapply(1:length(allocation_trace),
                        function(iter){
                          
                          
                          
                          psm.empty <- matrix(0, 
                                              nrow = length(allocation_trace[[1]]),
                                              ncol = length(allocation_trace[[1]]))
                          
                          for(i in 1:length(allocation_trace[[1]])){
                            
                            psm.empty[i,] <- psm.empty[i,] + ifelse(allocation_trace[[iter]] == allocation_trace[[iter]][i],
                                                                    1,
                                                                    0)
                          }
                          
                          psm.empty
                          
                        })
  
  
  # Sum
  loop.result <- Reduce('+', loop.result)
  psm <- loop.result/length(allocation_trace)
  
  return(psm)
}
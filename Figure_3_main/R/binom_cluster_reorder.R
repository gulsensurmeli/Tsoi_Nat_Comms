
# Function to reorder clusters
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





# Function to reorder clusters
k_means_reorder <- function(Y,
                            Z){
  
  
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
  
  
  return(Z_new)
}


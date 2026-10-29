
# Euclidean distance with log transformation
k_means_model <- function(Y,
                          k){

  Y_cbind <- do.call(cbind, Y)

  # Projection probabilities
  df <- t(apply(t(Y_cbind), 1, function(x){return(x/sum(x))}))

  # Number of neurons in each cluster
  C <- sapply(1:M, function(m) ncol(Y[[m]]))

  # Cumulative number of neurons in a vector form
  C_cumsum <- c(0, cumsum(C))

  # Number of mouse
  M <- length(Y)
  
  # Neuron allocation
  allocation <- kmeans(log(df+1), k, nstart = 50, iter.max = 100)

  # Neuron allocations
  Z <- lapply(1:M,
              function(m) allocation$cluster[(C_cumsum[m]+1):C_cumsum[m+1]])

  return(Z)

}

cosine_normalize <- function(mat) {

  # Normalize each row by its magnitude
  return(mat / sqrt(rowSums(mat^2)))
}

# Cosine distance
k_medoids_model <- function(Y,
                            k){

  # Length of input data
  M <- length(Y)
  
  C <- sapply(1:M, function(m) ncol(Y[[m]]))

  # Normalize data
  data_normalized <- lapply(1:M,
                            function(m){

                                normalized_mat <- apply(Y[[m]], 2, function(col) col / max(col))
                                normalized_mat
                            })

  # Combine normalized data by columns
  data_normalized_bind <- t(do.call(cbind, data_normalized))

  # Cosine distance
  data_normalized_bind <- cosine_normalize(mat = data_normalized_bind)

  # Allocations
  kmeans_result <- kmeans(data_normalized_bind, centers = k, iter.max = 100, nstart = 50)

  # Cumulative number of neurons in a vector form
  C_cumsum <- c(0, cumsum(C))

  # Neuron allocations
  Z <- lapply(1:M,
              function(m) kmeans_result$cluster[(C_cumsum[m]+1):C_cumsum[m+1]])

  return(Z)

}



#' Function to reorder clusters found by k-means.
#'
#' @param Y A list of input dataset.
#' @param Z A list of neuron allocations.
#' @return A list of neuron allocations after reordering.
#' @export

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



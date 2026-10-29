
##---------------- MCMC Simulation for Dirichlet Multinomial --------------

ddirmultinomial <- function(x, n, alpha, log = FALSE){

  # ind = alpha>0
  # 
  # # zero density if x_r>0 and alpha_r=0
  # if(sum(x[!ind])>0){ log.prob = -Inf}
  # else{
  #   log.prob <- lgamma(sum(alpha[ind])) + lgamma(n+1) - lgamma(n+sum(alpha[ind])) + sum(lgamma(x[ind] + alpha[ind]) - lgamma(alpha[ind]) - lgamma(x[ind]+1))
  # }  
  log.prob <- lgamma(sum(alpha)) + lgamma(n+1) - lgamma(n+sum(alpha)) + sum(lgamma(x + alpha) - lgamma(alpha) - lgamma(x+1))
  
  if(isFALSE(log)){

    return(exp(log.prob))
    
  }else{

    return(log.prob)
  }
}


##--------------- Simulation of allocations ----------------
allocation_variables_dirmult_mcmc <- function(omega_J_M,
                                              q_star_1_J,
                                              gamma_1_J_star,
                                              Y){

  # Create an empty list
  M <- length(Y)
  C <- unlist(lapply(Y, ncol))
  J <- nrow(omega_J_M)
  R <- ncol(q_star_1_J)

  ## Initialize Z
  Z <- NULL
  allocation.prob <- NULL

  for(m in 1:M){

    result.rbind <- lapply(1:C[m],
                           function(c){
                             
                             
                             ##-- log_p_tilde
                             log.p.tilde <- sapply(1:J, function(j){
                               
                               ddirmultinomial(Y[[m]][,c],
                                               n = sum(Y[[m]][,c]),
                                               alpha = q_star_1_J[j,]*gamma_1_J_star[j],
                                               log = TRUE)+
                                 log(omega_J_M[j,m])
                             })
                             
                             ##-- logK
                             logK <- -max(log.p.tilde)
                             
                             ##-- probability
                             allo.prob <- exp(log.p.tilde + logK)/sum(exp(log.p.tilde + logK))
                             
                             ##-- allocation
                             Z_cd <- extraDistr::rcat(n = 1,
                                                      prob = allo.prob)
                             
                             Z_cd_prob <- allo.prob[Z_cd]
                             
                             ##-- prepare for output
                             c(Z_cd, Z_cd_prob)
                             
                           })
    
    result.rbind <- do.call(rbind, result.rbind)


    ##-- Store Z
    Z[[m]] <- result.rbind[,1]

    ##-- Store allocation probability
    allocation.prob[[m]] <- result.rbind[,2]

  }

  ##-- Return Z
  return(list('Z' = Z,
              'allocation.prob' = allocation.prob))
}

##------------ Simulation of Data-Specific component probabilities --------------------
dataset_specific_mcmc <- function(Z,
                                  omega,
                                  alpha){

  J <- length(omega)
  M <- length(Z)

  loop.result <- lapply(1:M, function(m) {

    parameter <- as.vector(table(factor(Z[[m]], levels = c(1:J)))) + alpha*omega
    
    omega_J_M <- extraDistr::rdirichlet(n = 1,
                                             alpha = parameter)
    

    # shifted the zero values
    omega_J_M <- ifelse(omega_J_M == 0, .Machine$double.eps, omega_J_M)
    omega_J_M <- omega_J_M/rowSums(omega_J_M)

    return(omega_J_M)
  })

  omega_J_M <- matrix(unlist(loop.result), nrow=J, ncol=M)

  # Returning omega_J_M
  return(omega_J_M)
}



##----------- Simulation of Component probabilities ----------------
component_log_prob <- function(omega,
                               omega_J_M,
                               alpha_zero,
                               alpha){

  J <- nrow(omega_J_M)
  M <- ncol(omega_J_M)

  lprod <- sum((alpha_zero/J-1)*log(omega)) +

    sum(unlist(lapply(1:M, function(m){

      sum(alpha*omega*log(omega_J_M[,m])-lgamma(alpha*omega))
    })))

  # outputs
  return(lprod)
}

component_probabilities_mcmc <- function(omega,
                                         omega_J_M,
                                         alpha_zero,
                                         alpha,
                                         covariance,
                                         mean_x,
                                         tilde_s,
                                         iter_num,
                                         adaptive_prop = 0.01){

  J <- nrow(omega_J_M)
  M <- ncol(omega_J_M)

  # Define the inputs
  omega_old <- omega
  covariance_old <- covariance
  mean_X_d_old <- mean_x
  tilde_s_old <- tilde_s

  X_old <- log(omega_old[1:J-1]/omega_old[J]) # Length = J-1

  # Specify the iteration number
  n <- iter_num

  # Adaptive step
  if(n <= 100){

    X_new <- rmvnorm(n = 1,
                     mean = X_old,
                     sigma = diag(x=0.01, nrow=J-1, ncol=J-1))

  }else{

    X_new <- rmvnorm(n = 1,
                     mean = X_old,
                     sigma = 2.4^2/(J-1)*(covariance_old + adaptive_prop*diag(1, nrow = J-1, ncol = J-1)))

  }

  # Compute omega_new (Length = J) from X_new
  omega_new <- c(exp(X_new)/(1+sum(exp(X_new))),1/(1+sum(exp(X_new))))
  
  # shifted the zero values
  omega_new <- ifelse(omega_new == 0, .Machine$double.eps, omega_new)
  omega_new <- omega_new/sum(omega_new)

  # Compute acceptance probability
  log_acceptance <- component_log_prob(omega_new,
                                       omega_J_M,
                                       alpha_zero,
                                       alpha) -
    component_log_prob(omega_old,
                       omega_J_M,
                       alpha_zero,
                       alpha) +

    sum(log(omega_new)) - sum(log(omega_old))

    if(!is.finite(log_acceptance)){
    print('issues with inf/nan')
  }
  
  # Random Bernoulli
  outcome <- rbinom(n = 1,
                    size = 1,
                    prob=min(1, exp(log_acceptance)))

  if(outcome == 0){

    X_new <- X_old
    omega_new <- omega_old
  }

  # Update Covariance Structure
  tilde_s_new <- tilde_s_old + matrix(X_new, ncol = 1)%*%matrix(X_new, nrow = 1)
  mean_x_new <- mean_X_d_old*(1-1/n) + 1/n*matrix(X_new, nrow = 1)
  covariance_new <- 1/(n-1)*tilde_s_new - n/(n-1)*t(mean_x_new)%*%mean_x_new

  # Output
  return(list('omega_new' = omega_new,
              'tilde_s_new' = tilde_s_new,
              'mean_x_new' = mean_x_new,
              'covariance_new' = covariance_new,
              'accept' = outcome))
}


## ----------------------- Simulation of Concentration parameter alpha -----------------------
alpha_log_prob <- function(omega_J_M,
                           omega,
                           alpha,
                           a_alpha,
                           b_alpha){

  M <- ncol(omega_J_M)

  # lprod <- -alpha + M*lgamma(alpha) +
  # 
  #   sum(unlist(lapply(1:M, function(m) {
  #     sum(alpha*omega*log(omega_J_M[,m])-lgamma(alpha*omega))
  #   })))
   lprod <- -alpha*b_alpha + M*lgamma(alpha) + (a_alpha-1)*log(alpha) +
   
     sum(unlist(lapply(1:M, function(m) {
       sum(alpha*omega*log(omega_J_M[,m])-lgamma(alpha*omega))
     })))

  # output
  return(lprod)
}


alpha_mcmc <- function(omega_J_M,
                       omega,
                       alpha,
                       a_alpha,
                       b_alpha,
                       X_mean,
                       M_2,
                       variance,
                       iter_num,
                       adaptive_prop){

  M <- ncol(omega_J_M)

  # Defining the inputs
  alpha_old <- alpha
  X_old <- log(alpha_old)
  X_mean_old <- X_mean
  M_2_old <- M_2
  variance_old <- variance

  # Defining the dimensions
  n <- iter_num

  if(n <= 100){

    X_new <- rnorm(n = 1,
                   mean = X_old,
                   sd = 0.01)
  }else{

    X_new <- rnorm(n = 1,
                   mean = X_old,
                   sd = sqrt(2.4^2*(variance_old + adaptive_prop)))
  }

  # alpha_new
  alpha_new <- exp(X_new)

  # Compute log acceptance probability
  log_acceptance <- alpha_log_prob(omega_J_M = omega_J_M,
                                   omega = omega,
                                   alpha = alpha_new,
                                   a_alpha =  a_alpha,
                                   b_alpha =  b_alpha) -

    alpha_log_prob(omega_J_M = omega_J_M,
                   omega = omega,
                   alpha = alpha_old,
                   a_alpha =  a_alpha,
                   b_alpha =  b_alpha) +

    log(alpha_new) - log(alpha_old)
  
  if(!is.finite(log_acceptance)){
    print('issues with inf/nan')
  }

  # Random Bernoulli
  outcome <- rbinom(n = 1,
                    size = 1,
                    prob = min(1, exp(log_acceptance)))

  
  if(outcome == 0){

    X_new <- X_old
    alpha_new <- alpha_old
  }

  # Update Covariance Structure
  X_mean_new <- (1-1/n)*X_mean_old + 1/n*X_new
  # M_2_new <- M_2_old + (X_new-X_mean_old)*(X_new-X_mean_new)
  # variance_new <- 1/(n-1)*M_2_new
  M_2_new <- M_2_old +X_new^2
  variance_new <- 1/(n-1)*M_2_new - (1-1/n)*X_mean_new^2  

  # Output
  return(list('alpha_new' = alpha_new,
              'X_mean_new' = X_mean_new,
              'M_2_new' = M_2_new,
              'variance_new' = variance_new,
              'accept' = outcome))
}

##---------------------- Simulation of concentration parameter alpha_zero -----------------
alpha_zero_log_prob <- function(omega,
                                alpha_zero,
                                a_alpha0,
                                b_alpha0){

  # dimension
  J <- length(omega)

  # log probability
  lprob <- -b_alpha0*alpha_zero + lgamma(alpha_zero) - J*lgamma(alpha_zero/J) + sum(alpha_zero/J*log(omega)) + (a_alpha0-1)*log(alpha_zero)

  # output
  return(lprob)
}

alpha_zero_mcmc <- function(omega,
                            alpha_zero,
                            a_alpha0,
                            b_alpha0,
                            X_mean,
                            M_2,
                            variance,
                            iter_num,
                            adaptive_prop){

  # dimension
  J <- length(omega)
  n <- iter_num

  # previous values
  alpha_zero_old <- alpha_zero
  X_old <- log(alpha_zero_old)
  variance_old <- variance
  M_2_old <- M_2
  X_mean_old <- X_mean

  # Simulation
  if(n <= 100){

    X_new <- rnorm(n = 1,
                   mean = X_old,
                   sd = 0.01)
  }else{

    X_new <- rnorm(n = 1,
                   mean = X_old,
                   sd = sqrt(2.4^2*(variance_old + adaptive_prop)))
  }

  # alpha_zero_new
  alpha_zero_new <- exp(X_new)

  # log acceptance probability
  log_acceptance <- alpha_zero_log_prob(omega = omega,
                                        alpha_zero = alpha_zero_new, 
                                        a_alpha0 = a_alpha0,
                                        b_alpha0 = b_alpha0) -

    alpha_zero_log_prob(omega = omega,
                        alpha_zero = alpha_zero_old,
                        a_alpha0 = a_alpha0,
                        b_alpha0 = b_alpha0) +

    log(alpha_zero_new) - log(alpha_zero_old)
  
  if(!is.finite(log_acceptance)){
    print('issues with inf/nan')
  }

  # Random Bernoulli
  outcome <- rbinom(n = 1,
                    size = 1,
                    prob = min(1, exp(log_acceptance)))


  
  if(outcome == 0){

    X_new <- X_old
    alpha_zero_new <- alpha_zero_old
  }

  # Update Covariance Structure
  X_mean_new <- (1-1/n)*X_mean_old + 1/n*X_new
  # M_2_new <- M_2_old + (X_new-X_mean_old)*(X_new-X_mean_new)
  # variance_new <- 1/(n-1)*M_2_new
  M_2_new <- M_2_old +X_new^2
  variance_new <- 1/(n-1)*M_2_new - (1-1/n)*X_mean_new^2

  # Output
  return(list('alpha_zero_new' = alpha_zero_new,
              'X_mean_new' = X_mean_new,
              'M_2_new' = M_2_new,
              'variance_new' = variance_new,
              'accept' = outcome))
}

##------------------------------ Gamma ----------------------------------

gamma_logprob <- function(Y,
                          Z,
                          q_j_star,
                          gamma_j_star,
                          a_gamma,
                          b_gamma,
                          j){
  
  # Number of regions
  R <- nrow(Y[[1]])

  # Bind Y by columns
  Y_bind <- do.call(cbind,Y)

  # Initialize the log probability - gamma prior
  log_prob <- (a_gamma-1)*log(gamma_j_star) - b_gamma*gamma_j_star

  # Find the cell index of the ones with allocation equal to j
  cell_index_j <- which(unlist(Z)==j)

  # Sum of neuron counts across all regions for these cells
  N_CM_j <- colSums(matrix(Y_bind[,cell_index_j],
                           nrow = R,
                           ncol = length(cell_index_j)))


  # log - probability
  log_prob <- log_prob + sum(lgamma(gamma_j_star) + lgamma(N_CM_j + 1) - lgamma(N_CM_j + gamma_j_star) +

                               colSums(lgamma(matrix(Y_bind[,cell_index_j],
                                                     nrow = R,
                                                     ncol= length(cell_index_j)) + matrix(rep(q_j_star*gamma_j_star,
                                                                                            length(cell_index_j)),
                                                                                        nrow = R)
                               )) - sum(lgamma(q_j_star*gamma_j_star)) - colSums(lgamma(matrix(Y_bind[,cell_index_j] + 1,
                                                                                               nrow = R,
                                                                                               ncol = length(cell_index_j)))
                                                                                 )
  )

  return(log_prob)

}

gamma_mcmc <- function(Y,
                       Z,
                       q_star_1_J,
                       gamma_1_J_star,
                       a_gamma,
                       b_gamma,
                       X_mean,
                       M_2,
                       variance,
                       iter_num,
                       adaptive_prop){

  J <- nrow(q_star_1_J)
  n <- iter_num
  R <- ncol(q_star_1_J)

  # previous values
  gamma_1_J_star_old <- gamma_1_J_star
  X_old <- log(gamma_1_J_star_old)
  variance_old <- variance
  M_2_old <- M_2
  X_mean_old <- X_mean

  # Update
  X_mean_new <- NULL
  M_2_new <- NULL
  variance_new <- NULL
  gamma_count <- 0
  gamma_1_J_star_new <- NULL
  
  # For each component j
  for(j in 1:J){
    
    if(length(which(unlist(Z)==j)) != 0){
      
      if(n < 100){
        
        X_j_new <- rnorm(n = 1,
                         mean = X_old[j],
                         sd = 0.01)
        
        
      }else{
        
        X_j_new <- rnorm(n = 1,
                         mean = X_old[j],
                         sd = sqrt(2.4^2*(variance_old[j] + adaptive_prop)))
        
      }
      
      if(X_old[j]<log(2)){
        print(paste('old =', exp(X_old[j]), ', new =', exp(X_j_new)))
      }
      
      # gamma_j_new
      gamma_j_new <- exp(X_j_new)
      
      # Acceptance probability
      log.accept.prob <- gamma_logprob(Y = Y,
                                       Z = Z,
                                       q_j_star = q_star_1_J[j,],
                                       gamma_j_star = gamma_j_new,
                                       a_gamma = a_gamma,
                                       b_gamma = b_gamma,
                                       j = j)-
        
        gamma_logprob(Y = Y,
                      Z = Z,
                      q_j_star = q_star_1_J[j,],
                      gamma_j_star = gamma_1_J_star_old[j],
                      a_gamma = a_gamma,
                      b_gamma = b_gamma,
                      j = j)-
        
        log(gamma_1_J_star_old[j]) + log(gamma_j_new)
      
      accept.prob <- exp(log.accept.prob)
      
    }else{
      
      # If component j is empty
      
      # Simulate from the prior
      gamma_j_new <- rgamma(n = 1,
                            shape = a_gamma,
                            rate = b_gamma)
      
      # X_new
      X_j_new <- log(gamma_j_new)
      
      # Always accept
      accept.prob <- 1
      
    }
    
    if(!is.finite(accept.prob)){
      print('issues with inf/nan')
    }
 
    # Random Bernoulli
    outcome <- rbinom(n = 1,
                      size = 1,
                      prob = min(1,accept.prob))

    # If outcome is to reject
    if(outcome == 0){
      
      gamma_j_new <- gamma_1_J_star_old[j]
      X_j_new <- X_old[j]
    }
    
    X_mean_new[j] <- (1-1/n)*X_mean_old[j] + 1/n*X_j_new
    #M_2_new[j] <- M_2_old[j] + (X_j_new-X_mean_old[j])*(X_j_new-X_mean_new[j])
    #variance_new[j] <- 1/(n-1)*M_2_new[j]
    M_2_new[j] <- M_2_old[j] +X_j_new^2 
    variance_new[j] <- 1/(n-1)*M_2_new[j] - (1-1/n)*X_mean_new[j]^2
    gamma_count <- gamma_count + outcome
    gamma_1_J_star_new[j] <- gamma_j_new
    
  }
  


  # Return
  return(list('gamma_1_J_star_new' = gamma_1_J_star_new,
              'variance_gamma_new' = variance_new,
              'M_2_gamma_new' = M_2_new,
              'X_mean_gamma_new' = X_mean_new,
              'gamma_count' = gamma_count))

}

##-------------------------- q_star ------------------------------


q_star_logprob <- function(Y,
                           Z,
                           q_j_star,
                           gamma_j_star,
                           alpha_h,
                           j){


  # Bind Y by columns
  Y_bind <- do.call(cbind,Y)
  
  
  # Find the cell index of the ones with allocation equal to j
  cell_index_j <- which(unlist(Z)==j)
  
  # Sum of neuron counts across all regions for these cells
  N_CM_j <- colSums(as.matrix(Y_bind[,cell_index_j]))
  
  # Number of regions
  R <- nrow(Y[[1]])

  # Initialize the log probability - Dirichlet prior
  log_prob <- sum((alpha_h-1)*log(q_j_star))

  # log probability - distribution of Y
  log_prob <- log_prob + sum(colSums(lgamma(as.matrix(Y_bind[,cell_index_j]) + matrix(rep(q_j_star*gamma_j_star,
                                                                                          length(cell_index_j)),
                                                                                      nrow = R)
  ))) - sum(lgamma(q_j_star*gamma_j_star))*length(cell_index_j) -
    
    sum(colSums(lgamma(as.matrix(Y_bind[,cell_index_j] + 1))))

                 
  return(log_prob)

}

q_star_mcmc <- function(Y,
                        Z,
                        q_star_1_J,
                        gamma_1_J_star,
                        alpha_h,
                        covariance,
                        mean_x,
                        tilde_s,
                        iter_num,
                        adaptive_prop){

  J <- nrow(q_star_1_J)
  R <- ncol(q_star_1_J)
  n <- iter_num

  # Define the inputs
  q_star_1_J_old <- q_star_1_J
  covariance_old <- covariance
  mean_X_old <- mean_x
  tilde_s_old <- tilde_s
  X_old <- matrix(0, nrow = J, ncol = R-1) # Dimension = J x R-1
  for(j in 1:J){

    X_old[j,] <- log(q_star_1_J_old[j,1:(R-1)]/q_star_1_J_old[j,R])
  }

  q_star_1_J_new <- matrix(0, J, R)
  tilde_s_new <- NULL
  mean_x_new <- NULL
  covariance_new <- NULL
  q_star_count <- 0
  
  # For each cluster
  for(j in 1:J){
    
    
    if(length(which(unlist(Z)==j)) != 0){
      
      if(n < 100){
        
        # Length of R - 1
        X_j_new <- as.vector(rmvnorm(n = 1,
                                     mean = X_old[j,],
                                     sigma = diag(0.01,nrow = R-1, ncol = R-1)))
        
      }else{
        
        # Length of R - 1
        X_j_new <- as.vector(rmvnorm(n = 1,
                                     mean = X_old[j,],
                                     sigma = (2.4^2)/(R-1)*(covariance_old[[j]] + adaptive_prop*diag(1,nrow = R-1, ncol = R-1))))
        
      }
      
      # q_star_new
      q_star_j_new <- c(exp(X_j_new)/(1+sum(exp(X_j_new))),1/(1+sum(exp(X_j_new))))
      # shifted the zero values
      q_star_j_new <- ifelse( q_star_j_new == 0, .Machine$double.eps,  q_star_j_new)
      q_star_j_new <-  q_star_j_new/sum(q_star_j_new)
      
      # Acceptance probability
      log.accept.prob <- q_star_logprob(Y = Y,
                                        Z = Z,
                                        q_j_star = q_star_j_new,
                                        gamma_j_star = gamma_1_J_star[j],
                                        alpha_h = alpha_h,
                                        j = j)-
        
        q_star_logprob(Y = Y,
                       Z = Z,
                       q_j_star = q_star_1_J_old[j,],
                       gamma_j_star = gamma_1_J_star[j],
                       alpha_h = alpha_h,
                       j = j)+
        
        sum(log(q_star_j_new)) - sum(log(q_star_1_J_old[j,]))
      
      
      accept.prob <- exp(log.accept.prob)
      
    }else{
      
      # If component j is empty
      
      # Simulate from the prior
      q_star_j_new <- as.vector(extraDistr::rdirichlet(n = 1,
                                                       alpha = alpha_h))
      # shifted the zero values
      q_star_j_new <- ifelse( q_star_j_new == 0, .Machine$double.eps,  q_star_j_new)
      q_star_j_new <-  q_star_j_new/sum( q_star_j_new)
      
      # X_new
      X_j_new <- log(q_star_j_new[1:(R-1)]/q_star_j_new[R])
      
      # Always accept
      accept.prob <- 1
      
    }
    
    if(!is.finite(accept.prob)){
      print('issues with inf/nan')
    }
    
    
    # Random Bernoulli
    outcome <- rbinom(n = 1,
                      size = 1,
                      prob = min(1,accept.prob))

    # If outcome is to reject
    if(outcome == 0){
      
      q_star_j_new <- q_star_1_J_old[j,]
      X_j_new <- X_old[j,]
      
    }
    
    # Update
    tilde_s_new[[j]] <- tilde_s_old[[j]] + matrix(X_j_new, ncol = 1)%*%matrix(X_j_new, nrow = 1)
    mean_x_new[[j]] <- mean_X_old[[j]]*(1-1/n) + 1/n*matrix(X_j_new, nrow = 1)
    covariance_new[[j]] <- 1/(n-1)*tilde_s_new[[j]] - n/(n-1)*t(mean_x_new[[j]])%*%mean_x_new[[j]]
    q_star_count <- q_star_count + outcome
    q_star_1_J_new[j,] <- q_star_j_new
  }
  
  
  
  # Return
  return(list('tilde_s_new' = tilde_s_new,
              'mean_x_new' = mean_x_new,
              'covariance_new' = covariance_new,
              'q_star_count' = q_star_count,
              'q_star_1_J_new' = q_star_1_J_new))


}

##---------------------------- alpha_h -------------------------------------


alpha_h_logprob <- function(alpha_h,
                            q_star_1_J,
                            a_a,
                            b_a){

  R <- length(alpha_h)
  J <- nrow(q_star_1_J)

  # value <- -sum(alpha_h)
  value <- sum((a_a-1)*log(alpha_h)-b_a*alpha_h)

  for(j in 1:J){

    value <- value + lgamma(sum(alpha_h))+sum((alpha_h-1)*log(q_star_1_J[j,])-lgamma(alpha_h))
  }

  return(value)
}

alpha_h_mcmc <- function(alpha_h,
                         mean_X_alpha_h,
                         tilde_s_alpha_h,
                         q_star_1_J,
                         iter_num,
                         a_a,
                         b_a,
                         covariance,
                         adaptive_prop){

  ##-- Dimensions
  R <- length(alpha_h)
  J <- nrow(q_star_1_J)
  n <- iter_num

  ##-- Previous Covariance structure
  mean_X_alpha_h_old <- mean_X_alpha_h
  tilde_s_alpha_h_old <- tilde_s_alpha_h
  covariance_old <- covariance

  ##-- Previous alpha_h
  alpha_h_old <- alpha_h
  X_alpha_h_old <- log(alpha_h_old)

  ##-- Simulation
  if(n <= 100){

    X_alpha_h_new <- rmvnorm(n = 1,
                             mean = X_alpha_h_old,
                             sigma = diag(0.01,nrow = R, ncol = R))
  }else{

    X_alpha_h_new <- rmvnorm(n = 1,
                             mean = X_alpha_h_old,
                             sigma = (2.4^2)/R*(covariance_old + adaptive_prop*diag(1,nrow = R, ncol = R)))
  }

  ##-- alpha_h_new
  alpha_h_new <- exp(X_alpha_h_new)
  # shifted the zero values
  alpha_h_new <- ifelse( alpha_h_new == 0, .Machine$double.eps,  alpha_h_new)
  

  ##-- Acceptance probability
  accep.log.prob <- alpha_h_logprob(alpha_h = as.vector(alpha_h_new),
                                    q_star_1_J = q_star_1_J,
                                    a_a = a_a,
                                    b_a = b_a)-

    alpha_h_logprob(alpha_h = as.vector(alpha_h_old),
                    q_star_1_J = q_star_1_J,
                    a_a = a_a,
                    b_a = b_a)-

    sum(log(alpha_h_old)) + sum(log(alpha_h_new))
  
  if(!is.finite(accep.log.prob)){
    print('issues with inf/nan')
  }
  
  # Random Bernoulli
  outcome <- rbinom(n = 1,
                    size = 1,
                    prob = min(1,exp(accep.log.prob)))

  if(outcome == 0){

    X_alpha_h_new <- X_alpha_h_old
    alpha_h_new <- alpha_h_old
  }

  ##-- Update Covariance Structure
  tilde_s_alpha_h_new <- tilde_s_alpha_h_old + matrix(X_alpha_h_new, ncol=1) %*%
    matrix(X_alpha_h_new, nrow=1)

  mean_X_alpha_h_new <- mean_X_alpha_h_old*(1-1/n) + 1/n*matrix(X_alpha_h_new, nrow=1)

  covariance_new <- 1/(n-1)*tilde_s_alpha_h_new - n/(n-1)*t(mean_X_alpha_h_new)%*%mean_X_alpha_h_new

  ##-- Output
  return(list('alpha_h_new' = alpha_h_new,
              'accept' = outcome,
              'tilde_s_alpha_h_new' = tilde_s_alpha_h_new,
              'mean_X_alpha_h_new' = mean_X_alpha_h_new,
              'covariance_alpha_h_new' = covariance_new))
}

##---------------------------- a_h -------------------------------------


a_h_logprob <- function(a_h,
                        q_star_1_J,
                        a_a,
                        b_a,
                        mx){
  
  R <- length(mx)
  J <- nrow(q_star_1_J)
  
  value <- (a_a-1)*log(a_h)-b_a*a_h + 
    J*lgamma(a_h) -J*sum(lgamma(a_h*mx)) +
    a_h*sum(mx*colSums(log(q_star_1_J)))
  
  return(value)
}

a_h_mcmc <- function(a_h,
                     mean_X_a_h,
                     M_2_a_h,
                     q_star_1_J,
                     iter_num,
                     a_a,
                     b_a,
                     mx,
                     variance,
                     adaptive_prop){
  
  ##-- Dimensions
  R <- length(mx)
  J <- nrow(q_star_1_J)
  n <- iter_num
  
  ##-- Previous variance structure
  mean_X_a_h_old <- mean_X_a_h
  M_2_a_h_old <-  M_2_a_h
  variance_old <- variance
  
  ##-- Previous a_h
  a_h_old <- a_h
  X_a_h_old <- log(a_h_old)
  
  ##-- Simulation
  if(n <= 100){
    
    X_a_h_new <- rnorm(n = 1,
                       mean = X_a_h_old,
                       sd = 0.01)
  }else{
    
    X_a_h_new <- rnorm(n = 1,
                         mean = X_a_h_old,
                         sd = sqrt(2.4^2*(variance_old + adaptive_prop)))
  }
  
  ##-- a_h_new
  a_h_new <- exp(X_a_h_new)
  # shifted the zero values
  a_h_new <- ifelse( a_h_new == 0, .Machine$double.eps,  a_h_new)
  
  
  ##-- Acceptance probability
  accep.log.prob <- a_h_logprob(a_h = a_h_new,
                                q_star_1_J = q_star_1_J,
                                a_a = a_a,
                                b_a = b_a,
                                mx = mx)-
    
    a_h_logprob(a_h = a_h_old,
                q_star_1_J = q_star_1_J,
                a_a = a_a,
                b_a = b_a,
                mx = mx)-
    
    log(a_h_old) + log(a_h_new)
  
  if(!is.finite(accep.log.prob)){
    print('issues with inf/nan')
  }
  
  # Random Bernoulli
  outcome <- rbinom(n = 1,
                    size = 1,
                    prob = min(1,exp(accep.log.prob)))
  
  if(outcome == 0){
    
    X_a_h_new <- X_a_h_old
    a_h_new <- a_h_old
  }
  
  ##-- Update Covariance Structure
  M_2_a_h_new <- M_2_a_h_old + X_a_h_new^2
  
  mean_X_a_h_new <- mean_X_a_h_old*(1-1/n) + 1/n*X_a_h_new
  
  variance_new <- 1/(n-1)*M_2_a_h_new - n/(n-1)*mean_X_a_h_new^2
  
  ##-- Output
  return(list('a_h_new' = a_h_new,
              'accept' = outcome,
              'M_2_a_h_new' = M_2_a_h_new,
              'mean_X_a_h_new' = mean_X_a_h_new,
              'variance_a_h_new' = variance_new))
}

##---------------------------- split-and-merge -------------------------------------

split_merge_mcmc <- function(omega_J_M,
                             q_1_J_star,
                             gamma_1_J_star,
                             Y,
                             Z,
                             alpha_h,
                             a_gamma,
                             b_gamma,
                             alpha,
                             omega,
                             alpha_zero,
                             split_accept,
                             merge_accept){
  
  M <- length(Y)
  C <- unlist(lapply(Y, ncol))
  J <- nrow(omega_J_M)
  R <- ncol(q_1_J_star)
  
  # Initialize output
  gamma_1_J_star_new <- gamma_1_J_star 
  q_1_J_star_new <- q_1_J_star
  omega_J_M_new <- omega_J_M
  omega_new <- omega
  Z_new <- Z
  
  if(runif(1)<0.5){
    # Merge move
    
    # Compute the distance between the q star 
    # L1 but change to the symmetric cross-entropy?
    # TODO: restrict to only occupied components (or choose j1 and j2 based on likelihood if not too expensive)
    nj = unlist(lapply(1:J, function(j){sum(unlist(Z)==j)}))
    ind = which(nj>0)
    q_star_dist = as.matrix(dist(q_1_J_star[ind,], method = "manhattan"))
    temp = 10
    
    # Choose j1
    q_star_dist_min = unlist(lapply(1:length(ind), function(j){min(q_star_dist[j,-j])}))
    j1_ind = which(rmultinom(1,1, exp(-temp*q_star_dist_min))==1)
    j1 = ind[j1_ind]
    
    # Choose j2
    j2_ind = which(rmultinom(1,1, exp(-temp*q_star_dist[j1_ind,-j1_ind]))==1)
    j2_ind = j2_ind + (j2_ind>=j1_ind)
    j2 = ind[j2_ind]
    lproposal_j12= -temp*q_star_dist_min[j1_ind] - log(sum(exp(-temp*q_star_dist_min)))
    lproposal_j12 = lproposal_j12 + -temp*q_star_dist[j1_ind,j2_ind]- log(sum(exp(-temp*q_star_dist[j1_ind,-j1_ind])))
    lproposal_j21 = -temp*q_star_dist_min[j2_ind] - log(sum(exp(-temp*q_star_dist_min)))
    lproposal_j21 = lproposal_j21 + -temp*q_star_dist[j2_ind,j1_ind]- log(sum(exp(-temp*q_star_dist[j2_ind,-j2_ind])))
    lproposal_j = log(exp( lproposal_j12) + exp( lproposal_j21))
    
    # Propose merged values
    n_j1 = nj[j1]
    n_j2 = nj[j2]
    
    # Merged Z
    j_keep = j1*(n_j1>=n_j2) + j2*(n_j1<n_j2) 
    j_change = j1*(n_j1<=n_j2) + j2*(n_j1>n_j2)
    Z_merge = lapply(Z, function(z){z[z==j_change]=j_keep
                                     z})    
    # Z_merge = lapply(Z, function(z){z[z==j2]=j1
    #                                 z})
    q_1_J_star_merge = q_1_J_star
    gamma_1_J_star_merge = gamma_1_J_star
    omega_J_M_merge = omega_J_M
    
    # Merged q
    # # Option 1
    # q_j1_hat = q_1_J_star[j1,]*n_j1/(n_j1+ n_j2) + q_1_J_star[j2,]*n_j2/(n_j1+ n_j2)
    # q_1_J_star_merge[j1,] = extraDistr::rdirichlet(1, (n_j1+ n_j2)*q_j1_hat)
    # q_1_J_star_merge[j1,] = ifelse( q_1_J_star_merge[j1,] == 0, .Machine$double.eps,  q_1_J_star_merge[j1,])
    # q_1_J_star_merge[j1,] = q_1_J_star_merge[j1,]/sum(q_1_J_star_merge[j1,])
    # lproposal_qj1 = extraDistr::ddirichlet(q_1_J_star_merge[j1,], alpha = (n_j1+ n_j2)*q_j1_hat, log = TRUE)
    # lproposal_qj1_r = extraDistr::ddirichlet(q_1_J_star_merge[j1,], alpha = alpha_h, log = TRUE)
    # # propose from prior
    # q_1_J_star_merge[j2,] = extraDistr::rdirichlet(1, alpha = alpha_h)
    # q_1_J_star_merge[j2,] = ifelse( q_1_J_star_merge[j2,] == 0, .Machine$double.eps,  q_1_J_star_merge[j2,])
    # q_1_J_star_merge[j2,] = q_1_J_star_merge[j2,]/sum(q_1_J_star_merge[j2,])
    # lproposal_qj2 = extraDistr::ddirichlet(q_1_J_star_merge[j2,], alpha = alpha_h, log = TRUE)
    # lproposal_qj2_r = extraDistr::ddirichlet(q_1_J_star_merge[j2,], alpha = (n_j1+ n_j2)*q_j1_hat, log = TRUE)
    # Option 2
    q_1_J_star_merge[j_keep,] = q_1_J_star[j_keep,]
    q_1_J_star_merge[j_change,] = extraDistr::rdirichlet(1, alpha = alpha_h)
    q_1_J_star_merge[j_change,] = ifelse( q_1_J_star_merge[j_change,] == 0, .Machine$double.eps,  q_1_J_star_merge[j_change,])
    q_1_J_star_merge[j_change,] = q_1_J_star_merge[j_change,]/sum(q_1_J_star_merge[j_change,])
    lproposal_q = extraDistr::ddirichlet(q_1_J_star_merge[j_change,], alpha = alpha_h, log = TRUE)
    
    # Merged gamma
    # # Option 1
    # gamma_j1_hat = min(gamma_1_J_star[j1], gamma_1_J_star[j2])
    # gamma_1_J_star_merge[j1] = rgamma(1, shape = 10*gamma_j1_hat, rate = 10)
    # gamma_1_J_star_merge[j2] = rgamma(1,shape = a_gamma, rate = b_gamma)
    # lproposal_gammaj1 = dgamma(gamma_1_J_star_merge[j1], shape = 10*gamma_j1_hat, rate = 10,log = TRUE)
    # lproposal_gammaj1_r = dgamma(gamma_1_J_star_merge[j1], shape = a_gamma, rate = b_gamma,log = TRUE)
    # lproposal_gammaj2 = dgamma(gamma_1_J_star_merge[j2], shape = a_gamma, rate = b_gamma,log = TRUE)
    # lproposal_gammaj2_r = dgamma(gamma_1_J_star_merge[j2], shape = 10*gamma_j1_hat, rate = 10,log = TRUE)
    # Option 2
    gamma_1_J_star_merge[j_keep] = gamma_1_J_star[j_keep]
    gamma_1_J_star_merge[j_change] = rgamma(1,shape = a_gamma, rate = b_gamma)
    lproposal_gamma = dgamma(gamma_1_J_star_merge[j_change], shape = a_gamma, rate = b_gamma,log = TRUE)
    
    # Merged omega mouse specific
    omega_J_M_merge = matrix(0, J, M)
    lproposal_omega1 = 0
    # lproposal_omega1_r = 0
    for(m in 1:M){
      omega_temp = rbeta(1, shape1 = 1 + sum(Z[[m]]==j1)+sum(Z[[m]]==j2), shape2 = 1)
      omega_temp = c(omega_temp, 1-omega_temp)
      # omega_temp = extraDistr::rdirichlet(1, alpha = alpha*omega[c(j1,j2)] + c(sum(Z[[m]]==j1)+sum(Z[[m]]==j2),0))
      omega_temp = ifelse( omega_temp == 0, .Machine$double.eps,  omega_temp)
      omega_temp = ifelse( omega_temp >(1-10^(-16)), 1-10^(-16),  omega_temp)
      lproposal_omega1 = lproposal_omega1 + dbeta(omega_temp[1], shape1 = 1 + sum(Z[[m]]==j1)+sum(Z[[m]]==j2), shape2 = 1, log = TRUE)
      # lproposal_omega1 = lproposal_omega1 + dbeta(omega_temp[1], shape1 = 1 + sum(Z[[m]]==j1)+sum(Z[[m]]==j2), shape2 = 1, log = TRUE)
      # lproposal_omega1_r = lproposal_omega1_r + dbeta(omega_temp[1], shape2 = 1 + sum(Z[[m]]==j1)+sum(Z[[m]]==j2), shape1 = 1, log = TRUE)
      # # lproposal_omega = lproposal_omega + extraDistr::ddirichlet(omega_temp, alpha = alpha*omega[c(j1,j2)] + c(sum(Z[[m]]==j1)+sum(Z[[m]]==j2),0), log = TRUE)
      omega_temp = omega_temp*sum(omega_J_M[c(j1,j2),m])
      omega_out = omega_J_M[,m]
      # omega_out[c(j1,j2)] = omega_temp
      omega_out[c(j_keep,j_change)] = omega_temp
      omega_J_M_merge[,m] = omega_out/sum(omega_out)
    }
    
    # Merged omega 
    omega_merge = omega
    s = M/mean(omega_J_M_merge[j1,] + omega_J_M_merge[j2,])
    omega_temp = rbeta(1, shape1 = s*mean(omega_J_M_merge[j1,]), shape2 =  s*mean(omega_J_M_merge[j2,]))
    omega_temp = c(omega_temp, 1-omega_temp)
    omega_temp = ifelse( omega_temp == 0, .Machine$double.eps,  omega_temp)
    omega_temp = ifelse( omega_temp >(1-10^(-16)), 1-10^(-16),  omega_temp)
    lproposal_omega2 = dbeta(omega_temp[1], shape1 = s*mean(omega_J_M_merge[j1,]), shape2 = s*mean(omega_J_M_merge[j2,]), log = TRUE)
    omega_temp = omega_temp*sum(omega[c(j1,j2)])
    omega_merge[c(j1,j2)] = omega_temp
    omega_merge = omega_merge/sum(omega_merge)
    
    # Compute acceptance probability
    
    # Log of merge proposal
    # lproposal_merge_ordering = lproposal_qj1 +  lproposal_qj2 + lproposal_gammaj1 + lproposal_gammaj2 + lproposal_omega1 
    # lproposal_merge_ordering_reverse = lproposal_qj1_r +  lproposal_qj2_r + lproposal_gammaj1_r + lproposal_gammaj2_r + lproposal_omega1_r 
    # lproposal_merge = lproposal_j +  lproposal_merge_ordering - log(2) + log(1+0.5*exp(lproposal_merge_ordering_reverse- lproposal_merge_ordering)) + lproposal_omega2
    lproposal_merge = lproposal_j + lproposal_gamma + lproposal_q + lproposal_omega1 + lproposal_omega2
    
    # Log of split proposal
    temp = 1
    nj_merge = unlist(lapply(1:J, function(j){sum(unlist(Z_merge)==j)}))
    ind_merge = which(nj_merge>1)
    # lproposal_s_j1 = -temp*gamma_1_J_star_merge[j1] - log(sum(exp(-temp*gamma_1_J_star_merge[ind_merge])))
    lproposal_s_j1 = -temp*gamma_1_J_star_merge[j_keep] - log(sum(exp(-temp*gamma_1_J_star_merge[ind_merge])))
    
    # for q
    # Option 1 (more expensive)
    # data.cbind <- matrix(do.call(cbind, Y), nrow = R)
    # ind_merge_j1 = which(unlist(Z_merge)==j1)
    # data.cbind  = data.cbind[,ind_merge_j1]
    # pdata.j1 <- apply(data.cbind, 2, function(x){x/sum(x)})
    # pdata.j1.dist <- as.matrix(dist(t(pdata.j1), method = "manhattan"))
    # temp = 10
    # lprob.unnorm = temp*pdata.j1.dist
    # proposal_s_q = lapply(ind_merge_j1, function(c1){
    #   lprob_c1 = unlist(lapply(ind_merge_j1, function(c2){
    #     if(c1!=c2){
    #       lprob = -log(length(ind_merge_j1)) + lprob.unnorm[which(ind_merge_j1==c1),which(ind_merge_j1==c2)]- log(sum(exp(lprob.unnorm[which(ind_merge_j1==c1),-which(ind_merge_j1==c1)])))
    #       m = sum(c1>cumsum(C))+1
    #       c1m = c1 
    #       if(m>1){c1m = c1m- cumsum(C)[m-1]}
    #       lprob = lprob + extraDistr::ddirichlet(q_1_J_star[j1,], alpha = Y[[m]][,c1m] +1, log = TRUE)
    #       m = sum(c2>cumsum(C))+1
    #       c2m = c2 
    #       if(m>1){c2m = c2m- cumsum(C)[m-1]}
    #       lprob = lprob + extraDistr::ddirichlet(q_1_J_star[j2,], alpha = Y[[m]][,c2m] +1, log = TRUE)
    #       lprob
    #     }
    #   }))
    #   aux  = max(lprob_c1)
    #   sum(exp(lprob_c1-aux))*exp(aux)
    # })
    # lproposal_s_q = log(sum(unlist(proposal_s_q)))
    # # Option 2: cheaper but not as good proposal
    # s = 1
    # lproposal_s_q1  = extraDistr::ddirichlet(q_1_J_star[j1,], s*q_1_J_star_merge[j1,], log = TRUE)
    # lproposal_s_q2  = extraDistr::ddirichlet(q_1_J_star[j2,], s*q_1_J_star_merge[j1,], log = TRUE)
    # lproposal_s_q = lproposal_s_q1  +lproposal_s_q2
    # # Option 3 (in the middle between option 1 and 2)
    # ind_merge_j1 = which(unlist(Z_merge)==j1)
    # data.cbind <- matrix(do.call(cbind, Y), nrow = R)
    # data.cbind  = data.cbind[,ind_merge_j1]
    # pdata.j1 <- apply(data.cbind, 2, function(x){x/sum(x)})
    # pdata.j1.dist <- apply(t(pdata.j1), 1, function(p){ 
    #   as.matrix(dist(rbind(p, q_1_J_star_merge[j1,]), method = "manhattan"))[1,2]})
    # temp = 10
    # prob = exp(temp*pdata.j1.dist)/sum(exp(temp*pdata.j1.dist))
    # s = sum(unlist(Z_merge)==j1)/2
    # lproposal_s_q1  = extraDistr::ddirichlet(q_1_J_star[j1,], s*q_1_J_star_merge[j1,], log = TRUE)
    # lproposal_s_q2_vec = log(prob)+ apply(data.cbind, 2, function(y){
    #   extraDistr::ddirichlet(q_1_J_star[j2,], alpha = y +1, log = TRUE)})
    # lproposal_s_q2_max = max(lproposal_s_q2_vec)
    # lproposal_s_q2 = lproposal_s_q2_max + log(sum(exp(lproposal_s_q2_vec-lproposal_s_q2_max)))
    # lproposal_s_q2_r  = extraDistr::ddirichlet(q_1_J_star[j2,], s*q_1_J_star_merge[j1,], log = TRUE)
    # lproposal_s_q1_vec_r = log(prob)+ apply(data.cbind, 2, function(y){
    #   extraDistr::ddirichlet(q_1_J_star[j1,], alpha = y +1, log = TRUE)})
    # lproposal_s_q1_max_r = max(lproposal_s_q1_vec_r)
    # lproposal_s_q1_r = lproposal_s_q1_max_r + log(sum(exp(lproposal_s_q1_vec_r-lproposal_s_q1_max_r)))
    # ls_max = max(lproposal_s_q1 + lproposal_s_q2,lproposal_s_q1_r + lproposal_s_q2_r)
    # ls_min = min(lproposal_s_q1 + lproposal_s_q2,lproposal_s_q1_r + lproposal_s_q2_r)
    # lproposal_s_q = ls_max - log(2) - exp(1 + 0.5*exp(ls_min - ls_max))
    #  Option 4
    ind_merge_j1 = which(unlist(Z_merge)==j_keep)
    data.cbind <- matrix(do.call(cbind, Y), nrow = R)
    data.cbind  = as.matrix(data.cbind[,ind_merge_j1],nrow = R,ncol = length(ind_merge_j1))
    pdata.j1 <- apply(data.cbind, 2, function(x){x/sum(x)})
    pdata.j1.dist <- apply(t(pdata.j1), 1, function(p){
      as.matrix(dist(rbind(p, q_1_J_star_merge[j_keep,]), method = "manhattan"))[1,2]})
    temp = 10
    prob = exp(temp*pdata.j1.dist)/sum(exp(temp*pdata.j1.dist))
    lproposal_s_q_vec = log(prob)+ apply(data.cbind, 2, function(y){
      extraDistr::ddirichlet(q_1_J_star[j_change,], alpha = y +1, log = TRUE)})
    lproposal_s_q_max = max(lproposal_s_q_vec)
    lproposal_s_q = lproposal_s_q_max + log(sum(exp(lproposal_s_q_vec-lproposal_s_q_max)))
    
    # for gamma
    # s = 10
    # gamma_hat = gamma_1_J_star_merge[j1]* nj_merge[j1]/(nj_merge[j1]+ 2) + a_gamma/b_gamma*2/( nj_merge[j1]+ 2)
    # lproposal_s_gammaj1 = dgamma(gamma_1_J_star[j1], shape = s*gamma_hat, rate = s, log = TRUE)
    # lproposal_s_gammaj2 = dgamma(gamma_1_J_star[j2], shape = s*gamma_hat, rate = s, log = TRUE)
    # Option 2
    lproposal_s_gamma = dgamma(gamma_1_J_star[j_change], shape = a_gamma, rate = b_gamma, log = TRUE)
    
    # for omega mouse specific
    # Proposed values for omega
    lproposal_s_omega1 = 0 
    for(m in 1:M){
      omega_temp = omega_J_M[c(j1,j2),m]/sum(omega_J_M[c(j1,j2),m])
      omega_temp = ifelse( omega_temp == 0, .Machine$double.eps,  omega_temp)
      omega_temp = ifelse( omega_temp >(1-10^(-16)), 1-10^(-16),  omega_temp)
      lproposal_s_omega1 = lproposal_s_omega1 +  dbeta(omega_temp[1], shape1 = 1+ sum(Z[[m]]==j1), shape2 = 1+ sum(Z[[m]]==j2), log = TRUE)
      # # lproposal_s_omega = lproposal_s_omega + extraDistr::ddirichlet(omega_temp, alpha = alpha*omega[c(j1,j2)]+ c(sum(Z[[m]]==j1),sum(Z[[m]]==j2)), log = TRUE)
    }
    
    # for omega
    omega_temp = omega[c(j1,j2)]/sum(omega[c(j1,j2)])
    omega_temp = ifelse( omega_temp == 0, .Machine$double.eps,  omega_temp)
    omega_temp = ifelse( omega_temp >(1-10^(-16)), 1-10^(-16),  omega_temp)
    s = M/mean(omega_J_M[j1,] + omega_J_M[j2,])
    lproposal_s_omega2 =  dbeta(omega_temp[1], shape1 = s*mean(omega_J_M[j1,]), shape2 = s*mean(omega_J_M[j2,]), log = TRUE)
     
    lproposal_s = lproposal_s_j1 + lproposal_s_q + lproposal_s_gamma + lproposal_s_omega1 + lproposal_s_omega2
    
    # Ratio of log posterior
    # Part for y plus proposal of split for z 
    lpost_y_prop_split_z = 0 
    for(m in 1:M){
      
      ind_j1j2 = c(which(Z[[m]]==j1), which(Z[[m]]==j2))
      
      if(length(ind_j1j2)>0){
        result.rbind <- lapply(ind_j1j2,
                               function(c){
                                 
                                 
                                 ##-- log_p_tilde numerator
                                 # log.p.tilde.num <- ddirmultinomial(Y[[m]][,c],
                                 #                                    n = sum(Y[[m]][,c]),
                                 #                                    alpha = q_1_J_star_merge[j1,]*gamma_1_J_star_merge[j1],
                                 #                                    log = TRUE)
                                 log.p.tilde.num <- ddirmultinomial(Y[[m]][,c],
                                                                    n = sum(Y[[m]][,c]),
                                                                    alpha = q_1_J_star_merge[j_keep,]*gamma_1_J_star_merge[j_keep],
                                                                    log = TRUE)
                                 
                                 ##-- log_p_tilde demoninator
                                 log.p.tilde.denom <- ddirmultinomial(Y[[m]][,c],
                                                                      n = sum(Y[[m]][,c]),
                                                                      alpha = q_1_J_star[Z[[m]][c],]*gamma_1_J_star[Z[[m]][c]],
                                                                      log = TRUE)
                                 
                                 ##-- log_p_tilde for z proposal
                                 log.p.tilde.split <- sapply(c(j1,j2), function(j){
                                   
                                   ddirmultinomial(Y[[m]][,c],
                                                   n = sum(Y[[m]][,c]),
                                                   alpha = q_1_J_star[j,]*gamma_1_J_star[j],
                                                   log = TRUE)
                                 })
                                 
                                 ##-- logK
                                 logK <- -max(log.p.tilde.split)
                                 
                                 ##-- probability
                                 allo.prob <- exp(log.p.tilde.split + logK)/sum(exp(log.p.tilde.split + logK))
                                 
                                 
                                 log.p.tilde.num - log.p.tilde.denom + log(allo.prob[which(Z[[m]][c] == c(j1,j2))])
                               })
        
        result.rbind <- do.call(rbind, result.rbind)
        
        lpost_y_prop_split_z =  lpost_y_prop_split_z + sum(result.rbind)
      }
    }
      
      # Part for Z and omega
      lpost_zomega_num = log(omega_J_M_merge[j1,])*(unlist(lapply(Z_merge,function(z){sum(z==j1)})) + alpha*omega_merge[j1]- 1) + 
        log(omega_J_M_merge[j2,])*(unlist(lapply(Z_merge,function(z){sum(z==j2)})) + alpha*omega_merge[j2]- 1) -
        M*lgamma(omega_merge[j1]) -M*lgamma(omega_merge[j2]) 
      lpost_zomega_denom = log(omega_J_M[j1,])*(unlist(lapply(Z,function(z){sum(z==j1)})) + alpha*omega[j1]- 1) + 
        log(omega_J_M[j2,])*(unlist(lapply(Z,function(z){sum(z==j2)})) + alpha*omega[j2]- 1) -
        M*lgamma(omega[j1]) -M*lgamma(omega[j2]) 
      lpost_zomega = sum(lpost_zomega_num - lpost_zomega_denom)
      lpost_zomega = lpost_zomega + log(omega_merge[j1])*(alpha_zero/J-1) + log(omega_merge[j2])*(alpha_zero/J-1) -
        log(omega[j1])*(alpha_zero/J-1) - log(omega[j2])*(alpha_zero/J-1)
      
      
      # Part for q
      lpost_q_num = extraDistr::ddirichlet(q_1_J_star_merge[j1,], alpha = alpha_h, log = TRUE) + extraDistr::ddirichlet(q_1_J_star_merge[j2,], alpha = alpha_h, log = TRUE)
      lpost_q_denom = extraDistr::ddirichlet(q_1_J_star[j1,], alpha = alpha_h, log = TRUE) + extraDistr::ddirichlet(q_1_J_star[j2,], alpha = alpha_h, log = TRUE)
      lpost_q = lpost_q_num - lpost_q_denom
      
      # Part for gamma
      lpost_gamma_num = dgamma(gamma_1_J_star_merge[j1], shape = a_gamma, rate = b_gamma,log = TRUE) + dgamma(gamma_1_J_star_merge[j2],shape = a_gamma, rate = b_gamma,log = TRUE)
      lpost_gamma_denom = dgamma(gamma_1_J_star[j1], shape = a_gamma, rate = b_gamma,log = TRUE) + dgamma(gamma_1_J_star[j2], shape = a_gamma, rate =b_gamma,log = TRUE)
      lpost_gamma = lpost_gamma_num - lpost_gamma_denom
      
      lpost = lpost_y_prop_split_z + lpost_zomega + lpost_q + lpost_gamma
      
      accep.log.prob = lpost + lproposal_s - lproposal_merge
      
      print(paste('Merge accept prob:', round(min(1,exp(accep.log.prob)),4)))
      
      if(!is.finite(accep.log.prob)){
        print('issues with inf/nan')
      }
        
      # Random Bernoulli
      outcome <- rbinom(n = 1,
                        size = 1,
                        prob = min(1,exp(accep.log.prob)))
      
      if(outcome == 1){
        
        gamma_1_J_star_new <- gamma_1_J_star_merge 
        q_1_J_star_new <- q_1_J_star_merge
        omega_J_M_new <- omega_J_M_merge
        omega_new <- omega_merge
        Z_new <- Z_merge
        merge_accept = merge_accept +1
        
        print('Merge move accepted')
      }
    
  }else{
    
    # Split move
    # Check if empty components, if not can't propose split
    nj = unlist(lapply(1:J, function(j){sum(unlist(Z)==j)}))
    ind = which(nj>1)
    if (length(ind)<J){
      
      # Choose j
      temp = 1
      j1_ind = which(rmultinom(1,1, exp(-temp*gamma_1_J_star[ind]))==1)
      j1 = ind[j1_ind]
      j2 = which(nj==0)[1]
      lproposal_j1 = -temp*gamma_1_J_star[j1] - log(sum(exp(-temp*gamma_1_J_star[ind])))
      
      # Proposed values for q 
      q_1_J_star_split = q_1_J_star
      # # Option 1 (expensive)
      # # Pick randomly a neuron allocated to j1, and propose j1 based on the normalized counts
      # ind_j1 = which(unlist(Z)==j1)
      # c1 = sample(ind_j1,1)
      # m = sum(c1>cumsum(C))+1
      # c1m = c1 
      # if(m>1){c1m = c1m- cumsum(C)[m-1]}
      # q_1_J_star_split[j1,] = extraDistr::rdirichlet(1, alpha = Y[[m]][,c1m] +1)
      # q_1_J_star_split[j1,] = ifelse( q_1_J_star_split[j1,] == 0, .Machine$double.eps,  q_1_J_star_split[j1,])
      # q_1_J_star_split[j1,] = q_1_J_star_split[j1,]/sum(q_1_J_star_split[j1,])
      # # Pick a second neuron allocated to j1 with probability proportiontal to distance
      # # And propose j2 based on the normalized counts
      # data.cbind <- matrix(do.call(cbind, Y), nrow = R)
      # data.cbind  = data.cbind[,ind_j1]
      # pdata.j1 <- apply(data.cbind, 2, function(x){x/sum(x)})
      # pdata.j1.dist <- as.matrix(dist(t(pdata.j1), method = "manhattan"))
      # temp = 10
      # prob = exp(temp*pdata.j1.dist[which(ind_j1==c1),-which(ind_j1==c1)])/sum(exp(temp*pdata.j1.dist[which(ind_j1==c1),-which(ind_j1==c1)]))
      # c2 = sample(ind_j1[-which(ind_j1==c1)],1, prob = prob)
      # m = sum(c2>cumsum(C))+1
      # c2m = c2 
      # if(m>1){c2m = c2m- cumsum(C)[m-1]}
      # q_1_J_star_split[j2,] =  extraDistr::rdirichlet(1, alpha = Y[[m]][,c2m] +1)
      # q_1_J_star_split[j2,] = ifelse( q_1_J_star_split[j2,] == 0, .Machine$double.eps,  q_1_J_star_split[j2,])
      # q_1_J_star_split[j2,] = q_1_J_star_split[j2,]/sum(q_1_J_star_split[j2,])
      # 
      # # log proposal
      # lprob.unnorm = temp*pdata.j1.dist
      # proposal_q = lapply(ind_j1, function(c1){
      #   lprob_c1 = unlist(lapply(ind_j1, function(c2){
      #     if(c1!=c2){
      #       lprob = -log(length(ind_j1)) + lprob.unnorm[which(ind_j1==c1),which(ind_j1==c2)]- log(sum(exp(lprob.unnorm[which(ind_j1==c1),-which(ind_j1==c1)])))
      #       m = sum(c1>cumsum(C))+1
      #       c1m = c1 
      #       if(m>1){c1m = c1m- cumsum(C)[m-1]}
      #       lprob = lprob + extraDistr::ddirichlet(q_1_J_star_split[j1,], alpha = Y[[m]][,c1m] +1, log = TRUE)
      #       m = sum(c2>cumsum(C))+1
      #       c2m = c2 
      #       if(m>1){c2m = c2m- cumsum(C)[m-1]}
      #       lprob = lprob + extraDistr::ddirichlet(q_1_J_star_split[j2,], alpha = Y[[m]][,c2m] +1, log = TRUE)
      #       lprob
      #     }
      #   }))
      #   aux  = max(lprob_c1)
      #   sum(exp(lprob_c1-aux))*exp(aux)
      # })
      # lproposal_q = log(sum(unlist(proposal_q)))
      
      # # Option 2 (not as good proposal but cheaper)
      # s = 1
      # q_1_J_star_split[j1,] = extraDistr::rdirichlet(1, s*q_1_J_star[j1,])
      # q_1_J_star_split[j1,] = ifelse( q_1_J_star_split[j1,] == 0, .Machine$double.eps,  q_1_J_star_split[j1,])
      # q_1_J_star_split[j1,] = q_1_J_star_split[j1,]/sum(q_1_J_star_split[j1,])
      # q_1_J_star_split[j2,] = extraDistr::rdirichlet(1, s*q_1_J_star[j1,])
      # q_1_J_star_split[j2,] = ifelse( q_1_J_star_split[j2,] == 0, .Machine$double.eps,  q_1_J_star_split[j2,])
      # q_1_J_star_split[j2,] = q_1_J_star_split[j2,]/sum(q_1_J_star_split[j2,])
      # lproposal_q1  = extraDistr::ddirichlet(q_1_J_star_split[j1,], s*q_1_J_star[j1,], log = TRUE)
      # lproposal_q2  = extraDistr::ddirichlet(q_1_J_star_split[j2,], s*q_1_J_star[j1,], log = TRUE)
      # lproposal_q = lproposal_q1  +lproposal_q2
      
      # # Option 3 (in the middle between option 1 and 2)
      # # Pick a neuron allocated to j1, and based on its distance to qj1
      # ind_j1 = which(unlist(Z)==j1)
      # data.cbind <- matrix(do.call(cbind, Y), nrow = R)
      # data.cbind  = data.cbind[,ind_j1]
      # pdata.j1 <- apply(data.cbind, 2, function(x){x/sum(x)})
      # pdata.j1.dist <- apply(t(pdata.j1), 1, function(p){ 
      #   as.matrix(dist(rbind(p, q_1_J_star[j1,]), method = "manhattan"))[1,2]})
      # temp = 10
      # prob = exp(temp*pdata.j1.dist)/sum(exp(temp*pdata.j1.dist))
      # c1 = sample(ind_j1,1, prob = prob)
      # m = sum(c1>cumsum(C))+1
      # c1m = c1 
      # if(m>1){c1m = c1m- cumsum(C)[m-1]}
      # s = nj[j1]/2
      # q_1_J_star_split[j1,] = extraDistr::rdirichlet(1, s*q_1_J_star[j1,])
      # q_1_J_star_split[j1,] = ifelse( q_1_J_star_split[j1,] == 0, .Machine$double.eps,  q_1_J_star_split[j1,])
      # q_1_J_star_split[j1,] = q_1_J_star_split[j1,]/sum(q_1_J_star_split[j1,])
      # q_1_J_star_split[j2,] =  extraDistr::rdirichlet(1, alpha = Y[[m]][,c1m] +1)
      # q_1_J_star_split[j2,] = ifelse( q_1_J_star_split[j2,] == 0, .Machine$double.eps,  q_1_J_star_split[j2,])
      # q_1_J_star_split[j2,] = q_1_J_star_split[j2,]/sum(q_1_J_star_split[j2,])
      # # # log proposal
      # lproposal_q1  = extraDistr::ddirichlet(q_1_J_star_split[j1,], s*q_1_J_star[j1,], log = TRUE)
      # lproposal_q2_vec = log(prob)+ apply(data.cbind, 2, function(y){
      #   extraDistr::ddirichlet(q_1_J_star_split[j2,], alpha = y +1, log = TRUE)})
      # lproposal_q2_max = max(lproposal_q2_vec)
      # lproposal_q2 = lproposal_q2_max + log(sum(exp(lproposal_q2_vec-lproposal_q2_max)))
      # lproposal_q2_r  = extraDistr::ddirichlet(q_1_J_star_split[j2,], s*q_1_J_star[j1,], log = TRUE)
      # lproposal_q1_vec_r = log(prob)+ apply(data.cbind, 2, function(y){
      #   extraDistr::ddirichlet(q_1_J_star_split[j1,], alpha = y +1, log = TRUE)})
      # lproposal_q1_max_r = max(lproposal_q1_vec_r)
      # lproposal_q1_r = lproposal_q1_max_r + log(sum(exp(lproposal_q1_vec_r-lproposal_q1_max_r)))
      # lproposal_q = lproposal_q1 + lproposal_q2 - log(2) - exp(1 + 0.5*exp(lproposal_q1_r + lproposal_q2_r - lproposal_q1 - lproposal_q2))
      # 

      # Option 4 (in the middle between option 1 and 2)
      # Pick a neuron allocated to j1, and based on its distance to qj1
      ind_j1 = which(unlist(Z)==j1)
      data.cbind <- matrix(do.call(cbind, Y), nrow = R)
      data.cbind  = as.matrix(data.cbind[,ind_j1],nrow = R,ncol = length(ind_j1))
      pdata.j1 <- apply(data.cbind, 2, function(x){x/sum(x)})
      pdata.j1.dist <- apply(t(pdata.j1), 1, function(p){
        as.matrix(dist(rbind(p, q_1_J_star[j1,]), method = "manhattan"))[1,2]})
      temp = 10
      prob = exp(temp*pdata.j1.dist)/sum(exp(temp*pdata.j1.dist))
      c1 = sample(ind_j1,1, prob = prob)
      m = sum(c1>cumsum(C))+1
      c1m = c1
      if(m>1){c1m = c1m- cumsum(C)[m-1]}
      q_1_J_star_split[j1,] = q_1_J_star[j1,]
      q_1_J_star_split[j2,] =  extraDistr::rdirichlet(1, alpha = Y[[m]][,c1m] +1)
      q_1_J_star_split[j2,] = ifelse( q_1_J_star_split[j2,] == 0, .Machine$double.eps,  q_1_J_star_split[j2,])
      q_1_J_star_split[j2,] = q_1_J_star_split[j2,]/sum(q_1_J_star_split[j2,])

      # # log proposal
      lproposal_q_vec = log(prob)+ apply(data.cbind, 2, function(y){
        extraDistr::ddirichlet(q_1_J_star_split[j2,], alpha = y +1, log = TRUE)})
      lproposal_q_max = max(lproposal_q_vec)
      lproposal_q = lproposal_q_max + log(sum(exp(lproposal_q_vec-lproposal_q_max)))
      
      # Proposed values for gamma
      # gamma_1_J_star_split = gamma_1_J_star
      # s = 10
      # gamma_hat = gamma_1_J_star[j1]*nj[j1]/(nj[j1]+ 2) + a_gamma/b_gamma*2/(nj[j1]+ 2)
      # gamma_1_J_star_split[j1] = rgamma(1, shape = s*gamma_hat, rate = s)
      # gamma_1_J_star_split[j2] = rgamma(1, shape = s*gamma_hat, rate = s)
      # lproposal_gammaj1 = dgamma(gamma_1_J_star_split[j1], shape = s*gamma_hat, rate = s, log = TRUE)
      # lproposal_gammaj2 = dgamma(gamma_1_J_star_split[j2], shape = s*gamma_hat, rate = s, log = TRUE)
      # Option 2
      gamma_1_J_star_split = gamma_1_J_star
      gamma_1_J_star_split[j1] = gamma_1_J_star[j1]
      gamma_1_J_star_split[j2] = rgamma(1, shape = a_gamma, rate = b_gamma)
      lproposal_gamma = dgamma(gamma_1_J_star_split[j2], shape = a_gamma, rate = b_gamma, log = TRUE)
      
      # Proposed values for Z
      lproposal_z = 0 
      Z_split = Z
      for(m in 1:M){
        
        ind_j1 = which(Z[[m]]==j1)
        
        if(length(ind_j1)>0){
          result.rbind <- lapply(ind_j1,
                                 function(c){
                                   
                                   
                                   ##-- log_p_tilde
                                   log.p.tilde <- sapply(c(j1,j2), function(j){
                                     
                                     ddirmultinomial(Y[[m]][,c],
                                                     n = sum(Y[[m]][,c]),
                                                     alpha = q_1_J_star_split[j,]*gamma_1_J_star_split[j],
                                                     log = TRUE)
                                   })
                                   
                                   ##-- logK
                                   logK <- -max(log.p.tilde)
                                   
                                   ##-- probability
                                   allo.prob <- exp(log.p.tilde + logK)/sum(exp(log.p.tilde + logK))
                                   
                                   ##-- allocation
                                   Z_cd <- extraDistr::rcat(n = 1,
                                                            prob = allo.prob)
                                   
                                   Z_cd_prob <- allo.prob[Z_cd]
                                   
                                   Z_cd <-c(j1,j2)[Z_cd]
                                   
                                   ##-- prepare for output
                                   c(Z_cd, Z_cd_prob)
                                   
                                 })
          
          result.rbind <- do.call(rbind, result.rbind)
          
          lproposal_z =  lproposal_z + sum(log(result.rbind[,2]))
          
          ##-- Store Z
          Z_split[[m]][ind_j1] <- result.rbind[,1]
        }
      }
      
      # Check if both are occupied
      nj1_s = sum(unlist(Z_split)==j1)
      nj2_s = sum(unlist(Z_split)==j2)
      if(nj1_s ==0 |nj2_s==0){
        print('bad split proposal')
      }else{
        # Proposed values for omega mouse specific
        omega_J_M_split = matrix(0, J, M)
        lproposal_omega1 = 0 
        for(m in 1:M){
          omega_temp = rbeta(1, shape1 = 1 + sum(Z_split[[m]]==j1), shape2 = 1 + sum(Z_split[[m]]==j2) )
          omega_temp = c(omega_temp, 1-omega_temp)
          # omega_temp = extraDistr::rdirichlet(1, alpha = alpha*omega[c(j1,j2)] + c(sum(Z_split[[m]]==j1),sum(Z_split[[m]]==j2)))
          omega_temp = ifelse( omega_temp == 0, .Machine$double.eps,  omega_temp)
          omega_temp = ifelse( omega_temp >(1-10^(-16)), 1-10^(-16),  omega_temp)
          lproposal_omega1 = lproposal_omega1 + dbeta(omega_temp[1], shape1 = 1 + sum(Z_split[[m]]==j1), shape2 = 1 + sum(Z_split[[m]]==j2), log = TRUE)
          # lproposal_omega = lproposal_omega + extraDistr::ddirichlet(omega_temp, alpha = alpha*omega[c(j1,j2)] + c(sum(Z_split[[m]]==j1),sum(Z_split[[m]]==j2)), log = TRUE)
          omega_temp = omega_temp*sum(omega_J_M[c(j1,j2),m])
          omega_out = omega_J_M[,m]
          omega_out[c(j1,j2)] = omega_temp
          omega_J_M_split[, m] = omega_out/sum(omega_out)
        }
        
        omega_split = omega
        s = M/mean(omega_J_M_split[j1,] + omega_J_M_split[j2,])
        omega_temp = rbeta(1, shape1 = s*mean(omega_J_M_split[j1,]), shape2 = s*mean(omega_J_M_split[j2,]))
        omega_temp = c(omega_temp, 1-omega_temp)
        omega_temp = ifelse( omega_temp == 0, .Machine$double.eps,  omega_temp)
        omega_temp = ifelse( omega_temp >(1-10^(-16)), 1-10^(-16),  omega_temp)
        lproposal_omega2 =  dbeta(omega_temp[1], shape1 = s*mean(omega_J_M_split[j1,]), shape2 = s*mean(omega_J_M_split[j2,]), log = TRUE)
        omega_temp = omega_temp*sum(omega[c(j1,j2)])
        omega_split[c(j1,j2)] = omega_temp
        omega_split = omega_split/sum(omega_split)
        
        # Compute acceptance probability
        
        # Log of split proposal
        # lproposal_split = lproposal_j1 +  lproposal_q + lproposal_gammaj1 + lproposal_gammaj2 + lproposal_z + lproposal_omega1 + lproposal_omega2
        lproposal_split = lproposal_j1 +  lproposal_q + lproposal_gamma + lproposal_z + lproposal_omega1 + lproposal_omega2
        
        # Log of merge proposal
        nj_s = unlist(lapply(1:J, function(j){sum(unlist(Z_split)==j)}))
        ind_s = which(nj_s>0)
        q_star_dist = as.matrix(dist(q_1_J_star_split[ind_s,], method = "manhattan"))
        temp = 10
        
        # Choose j1, j2
        q_star_dist_min = unlist(lapply(1:length(ind_s), function(j){min(q_star_dist[j,-j])}))
        lproposal_merge_j12= -temp*q_star_dist_min[which(ind_s ==j1)] - log(sum(exp(-temp*q_star_dist_min)))
        lproposal_merge_j12 = lproposal_merge_j12 + -temp*q_star_dist[which(ind_s ==j1),which(ind_s ==j2)]- log(sum(exp(-temp*q_star_dist[which(ind_s ==j1),-which(ind_s ==j1)])))
        lproposal_merge_j21= -temp*q_star_dist_min[which(ind_s ==j2)] - log(sum(exp(-temp*q_star_dist_min)))
        lproposal_merge_j21 = lproposal_merge_j21 + -temp*q_star_dist[which(ind_s ==j2),which(ind_s ==j1)]- log(sum(exp(-temp*q_star_dist[which(ind_s ==j2),-which(ind_s ==j2)])))
        lproposal_merge_j = log(exp(lproposal_merge_j12)+ exp(lproposal_merge_j12))
        
        # Merged q
        # q_j1_hat = q_1_J_star_split[j1,]*nj_s[j1]/(nj_s[j1]+ nj_s[j2]) + q_1_J_star_split[j2,]*nj_s[j2]/(nj_s[j1]+ nj_s[j2])
        # lproposal_merge_qj1 = extraDistr::ddirichlet(q_1_J_star[j1,], alpha = (nj_s[j1]+ nj_s[j2])*q_j1_hat, log = TRUE)
        # lproposal_merge_qj2 = extraDistr::ddirichlet(q_1_J_star[j2,], alpha = alpha_h, log = TRUE)
        # lproposal_merge_qj1_r = extraDistr::ddirichlet(q_1_J_star[j1,], alpha = alpha_h, log = TRUE)
        # lproposal_merge_qj2_r = extraDistr::ddirichlet(q_1_J_star[j2,], alpha = (nj_s[j1]+ nj_s[j2])*q_j1_hat, log = TRUE)
        j_m_keep = j1*(nj1_s>=nj2_s) + j2*(nj1_s<nj2_s) 
        j_m_change = j1*(nj1_s<=nj2_s) + j2*(nj1_s>nj2_s)
        if(j_m_keep != j1){lproposal_merge_q = -Inf}else{
          lproposal_merge_q = extraDistr::ddirichlet(q_1_J_star[j2,], alpha = alpha_h, log = TRUE)
        }
        
        # Merged gamma
        # gamma_j1_hat = min(gamma_1_J_star_split[j1], gamma_1_J_star_split[j2])
        # lproposal_merge_gammaj1 = dgamma(gamma_1_J_star[j1], shape = 10*gamma_j1_hat, rate = 10,log = TRUE)
        # lproposal_merge_gammaj2 = dgamma(gamma_1_J_star[j2],shape =a_gamma, rate = b_gamma,log = TRUE)
        # lproposal_merge_gammaj1_r = dgamma(gamma_1_J_star[j1], shape =a_gamma, rate = b_gamma,log = TRUE)
        # lproposal_merge_gammaj2_r = dgamma(gamma_1_J_star[j2],shape = 10*gamma_j1_hat, rate = 10,log = TRUE)
        # Option 2
        if(j_m_keep != j1){lproposal_merge_gamma = -Inf}else{
          lproposal_merge_gamma = dgamma(gamma_1_J_star[j2], shape = a_gamma, rate = b_gamma,log = TRUE)
        }

        # Merged omega mouse specific
        lproposal_merge_omega1 = 0
        for(m in 1:M){
          # omega_temp = omega_J_M[c(j1,j2),m]/sum(omega_J_M[c(j1,j2),m])
          omega_temp = omega_J_M[c(j_m_keep,j_m_change),m]/sum(omega_J_M[c(j1,j2),m])
          omega_temp = ifelse( omega_temp == 0, .Machine$double.eps,  omega_temp)
          omega_temp = ifelse( omega_temp >(1-10^(-16)), 1-10^(-16),  omega_temp)
          lproposal_merge_omega1 = lproposal_merge_omega1 + dbeta(omega_temp[1], shape1 = 1 + sum(Z_split[[m]]==j1)+sum(Z_split[[m]]==j2), shape2 = 1, log = TRUE)
          # lproposal_merge_omega1_r = lproposal_merge_omega1 + dbeta(omega_temp[1], shape2 = 1 + sum(Z_split[[m]]==j1)+sum(Z_split[[m]]==j2), shape1 = 1, log = TRUE)
          # # lproposal_merge_omega = lproposal_merge_omega + extraDistr::ddirichlet(omega_temp, alpha = alpha*omega[c(j1,j2)] + c(sum(Z_split[[m]]==j1)+sum(Z_split[[m]]==j2),0), log = TRUE)
        }
        
        # Merged omega 
        omega_temp = omega[c(j1,j2)]/sum(omega[c(j1,j2)])
        omega_temp = ifelse( omega_temp == 0, .Machine$double.eps,  omega_temp)
        omega_temp = ifelse( omega_temp >(1-10^(-16)), 1-10^(-16),  omega_temp)
        s =M/mean(omega_J_M[j1,] + omega_J_M[j2,])
        lproposal_merge_omega2 = dbeta(omega_temp[1], shape1 = s*mean(omega_J_M[j1,]), shape2 = s*mean(omega_J_M[j2,]), log = TRUE)
        
        # Compute acceptance probability
        
        # Log of merge proposal
        # lproposal_merge_o1 = lproposal_merge_qj1 +  lproposal_merge_qj2 + lproposal_merge_gammaj1 + lproposal_merge_gammaj2 + lproposal_merge_omega1 
        # lproposal_merge_o2 = lproposal_merge_qj1_r +  lproposal_merge_qj2_r + lproposal_merge_gammaj1_r + lproposal_merge_gammaj2_r + lproposal_merge_omega1_r 
        # lm_max = max(lproposal_merge_o1,lproposal_merge_o2)
        # lm_min = min(lproposal_merge_o1,lproposal_merge_o2)
        # lproposal_merge = lproposal_merge_j +  lm_max - log(2) + log(1+0.5*exp(lm_min- lm_max)) + lproposal_merge_omega2
        lproposal_merge = lproposal_merge_j + lproposal_merge_q + lproposal_merge_gamma + lproposal_merge_omega1 + lproposal_merge_omega2
          
        # Ratio of log posterior
        lpost_y = 0 
        for(m in 1:M){
          
          ind_j1 = which(Z[[m]]==j1)
          
          if(length(ind_j1)>0){
            result.rbind <- lapply(ind_j1,
                                   function(c){
                                     
                                     
                                     ##-- log_p_tilde numerator
                                     log.p.tilde.num <- ddirmultinomial(Y[[m]][,c],
                                                                        n = sum(Y[[m]][,c]),
                                                                        alpha = q_1_J_star_split[Z_split[[m]][c],]*gamma_1_J_star_split[Z_split[[m]][c]],
                                                                        log = TRUE)
                                     
                                     ##-- log_p_tilde demoninator
                                     log.p.tilde.denom <- ddirmultinomial(Y[[m]][,c],
                                                                          n = sum(Y[[m]][,c]),
                                                                          alpha = q_1_J_star[j1,]*gamma_1_J_star[j1],
                                                                          log = TRUE)
                                     
                                     log.p.tilde.num - log.p.tilde.denom
                                   })
            
            result.rbind <- do.call(rbind, result.rbind)
            
            lpost_y =  lpost_y + sum(result.rbind)
          }
        }
        
        # Part for Z and omega
        lpost_zomega_num = log(omega_J_M_split[j1,])*(unlist(lapply(Z_split,function(z){sum(z==j1)})) + alpha*omega_split[j1]- 1) + 
          log(omega_J_M_split[j2,])*(unlist(lapply(Z_split,function(z){sum(z==j2)})) + alpha*omega_split[j2]- 1) -
          M*lgamma(omega_split[j1]) -M*lgamma(omega_split[j2]) 
        lpost_zomega_denom = log(omega_J_M[j1,])*(unlist(lapply(Z,function(z){sum(z==j1)}))+ alpha*omega[j1]- 1) + 
          log(omega_J_M[j2,])*(alpha*omega[j2]- 1) -
          M*lgamma(omega[j1]) -M*lgamma(omega[j2]) 
        lpost_zomega = sum(lpost_zomega_num - lpost_zomega_denom)
        lpost_zomega = lpost_zomega + log(omega_split[j1])*(alpha_zero/J-1) + log(omega_split[j2])*(alpha_zero/J-1) -
          log(omega[j1])*(alpha_zero/J-1) - log(omega[j2])*(alpha_zero/J-1)
        
        # Part for q
        lpost_q_num = extraDistr::ddirichlet(q_1_J_star_split[j1,], alpha = alpha_h, log = TRUE) + extraDistr::ddirichlet(q_1_J_star_split[j2,], alpha = alpha_h, log = TRUE)
        lpost_q_denom = extraDistr::ddirichlet(q_1_J_star[j1,], alpha = alpha_h, log = TRUE) + extraDistr::ddirichlet(q_1_J_star[j2,], alpha = alpha_h, log = TRUE)
        lpost_q = lpost_q_num - lpost_q_denom
        
        # Part for gamma
        lpost_gamma_num = dgamma(gamma_1_J_star_split[j1], shape = a_gamma, rate = b_gamma,log = TRUE) + dgamma(gamma_1_J_star_split[j2],shape = a_gamma, rate = b_gamma,log = TRUE)
        lpost_gamma_denom = dgamma(gamma_1_J_star[j1], shape = a_gamma, rate = b_gamma,log = TRUE) + dgamma(gamma_1_J_star[j2], shape = a_gamma, rate = b_gamma,log = TRUE)
        lpost_gamma = lpost_gamma_num - lpost_gamma_denom
        
        lpost = lpost_y + lpost_zomega + lpost_q + lpost_gamma
        
        accep.log.prob = lpost +  lproposal_merge -  lproposal_split
        
        print(paste('Split accept prob:', round(min(1,exp(accep.log.prob)),4)))
        
        outcome = 0
        if(!is.finite(accep.log.prob)){
          print('issues with inf/nan')
        }else{
          # Random Bernoulli
          outcome <- rbinom(n = 1,
                            size = 1,
                            prob = min(1,exp(accep.log.prob)))
        }
        
        if(outcome == 1){
          
          gamma_1_J_star_new <- gamma_1_J_star_split 
          q_1_J_star_new <- q_1_J_star_split
          omega_J_M_new <- omega_J_M_split
          omega_new <- omega_split
          Z_new <- Z_split
          split_accept = split_accept +1
          
          print('Split move accepted')
        }
      }
      }
    }
  
  ##-- Return Z
  return(list('Z' = Z_new,
              'omega_J_M' = omega_J_M_new,
              'omega' = omega_new,
              'q_1_J_star' = q_1_J_star_new,
              'gamma_1_J_star' = gamma_1_J_star_new,
              'split_accept' = split_accept,
              'merge_accept' = merge_accept))
}
#' MCMC Simulation
#'

mcmc_run_all <- function(Y,
                         J,
                         number_iter = 6000,
                         thinning = 5,
                         burn_in = 3000,
                         adaptive_prop = 0.001,
                         print_Z = FALSE,
                         iter_update = 100,
                         a_gamma,
                         b_gamma,
                         a_a = 1,
                         b_a= 1,
                         a_alpha = 1,
                         b_alpha= 1,
                         a_alpha0 = 1,
                         b_alpha0= 1,
                         Z.init = NULL){




  ##--------------- dimensions --------------------

  M <- length(Y)

  ##-- Need to make sure that all datasets have the same number of regions
  stopifnot(length(unique(sapply(1:M, function(m) nrow(Y[[m]]))))==1)

  ##-- Number of regions
  R <- nrow(Y[[1]])

  ##-- Number of neuron cells - vector of length M
  C <- unlist(lapply(1:M, function(m) ncol(Y[[m]])))

  #------------------------ Prepare for outputs -----------------

  Z_output <- NULL
  omega_J_M_output <- NULL
  omega_output <- NULL
  alpha_output <- 0
  alpha_zero_output <- 0
  allocation_prob_output <- NULL
  q_star_1_J_output <- NULL
  gamma_star_1_J_output <- NULL
  # alpha_h_output <- NULL
  a_h_output <- NULL

  #--------------------------------------------------------------------------------

  # Acceptance probability
  # acceptance_prob_list <- data.frame(omega = rep(0, number_iter-1),
  #                                    alpha = rep(0, number_iter-1),
  #                                    alpha_zero = rep(0, number_iter-1),
  #                                    q_star = rep(0, number_iter-1),
  #                                    gamma_star = rep(0, number_iter-1),
  #                                    alpha_h = rep(0, number_iter-1))
  acceptance_prob_list <- data.frame(omega = rep(0, number_iter-1),
                                     alpha = rep(0, number_iter-1),
                                     alpha_zero = rep(0, number_iter-1),
                                     q_star = rep(0, number_iter-1),
                                     gamma_star = rep(0, number_iter-1),
                                     a_h = rep(0, number_iter-1))


  #----------------------- Step 2: Set starting values ----------------------------
  if (is.null(Z.init)){
    Z_new <- lapply(1:M,
                    function(m) sample(1:J, size = C[m], replace = TRUE))
  }else{
    Z_new <- Z.init
  }
  
  # Empirical values based on counts
  omega_hat = function(z=z,J=J){
    omega = rep(1,J)
    tbz = table(z)
    omega[as.numeric(names(tbz))] = omega[as.numeric(names(tbz))] + tbz
    return(omega/sum(omega))
  }
  
  omega_J_M_new <- matrix(unlist(lapply(Z_new, omega_hat, J=J)), J, M, byrow=FALSE)
  
  omega_new <- omega_hat(unlist(Z_new), J=J)
  
  # Empirical values based on DP approximation
  alpha_zero_new  <- mean(unlist(lapply(Z_new,function(x){length(unique(x))}))/log(C))
  alpha_new <-length(unique(unlist(Z_new)))/log(sum(unlist(lapply(Z_new,function(x){length(unique(x))}))))
  
  # Hyperparameters for q
  # Compute empirical mean of projection strength
  data.j.prior <- matrix(do.call(cbind, Y), nrow = R)
  pdata.j.prior <- apply(data.j.prior, 2, function(x){x/sum(x)})
  data.j.mean <- rowMeans(pdata.j.prior)
  mx <- matrix(data.j.mean/sum(data.j.mean), nrow = 1)
  a_h_new = a_a/b_a
  alpha_h_new <- a_h_new*mx
  
  # Set q empirically based on Z_new
  q_star_1_J_new <- lapply(1:J,
                          function(j){
                            
                            if(length(which(unlist(Z_new)==j)) != 0){
                              
                              data.j.prior <- matrix(do.call(cbind, Y)[,which(unlist(Z_new)==j)],
                                                     nrow = R)
                              
                              pdata.j.prior <- apply(data.j.prior, 2, function(x){x/sum(x)})
                              data.j.mean <- rowMeans(pdata.j.prior)
                              
                              mx2 <- matrix(data.j.mean/sum(data.j.mean),
                                           nrow = 1)
                              
                            }else{
                              
                              mx2 <- matrix(rdirichlet(1, alpha_h_new), nrow = 1, ncol = R)
                            }
                            
                            mx2
                            
                          })
  
  q_star_1_J_new <- do.call(rbind, q_star_1_J_new)
  q_star_1_J_new <- ifelse(q_star_1_J_new == 0, .Machine$double.eps,  q_star_1_J_new)
  q_star_1_J_new <-  q_star_1_J_new/rowSums(q_star_1_J_new)
  
  # Set gamma empirically based on Z_new
  gamma_1_J_star_new  <- lapply(1:J,
                           function(j){
                             
                             if(length(which(unlist(Z_new)==j)) > 1){
                               
                               data.j.prior <- matrix(do.call(cbind, Y)[,which(unlist(Z_new)==j)],
                                                      nrow = R)
                               
                               pdata.j.prior <- apply(data.j.prior, 2, function(x){x/sum(x)})
                               data.j.var <- apply(pdata.j.prior,1,var)
                               data.j.mean <- rowMeans(pdata.j.prior)
                               
                               mx2 <- mean(data.j.mean*(1-data.j.mean)/(data.j.var+ 0.0001)  - 1)
                               mx2 <- mx2*(mx2>0) + (mx2<=0)
                               
                             }else{
                               
                               mx2 <- rgamma(1,shape = a_gamma, rate = b_gamma)
                             }
                             
                             mx2
                             
                           })
  
  gamma_1_J_star_new <- do.call(rbind, gamma_1_J_star_new)
  
  # alpha_h_new <- rep(a_a/b_a,R)


  # Total count of acceptance
  omega_count <- 0
  alpha_count <- 0
  alpha_zero_count <- 0
  q_star_count <- 0
  gamma_star_count <- 0
  # alpha_h_count <- 0
  a_h_count <- 0
  split_accept = 0
  merge_accept = 0

  #----------------------- Step 3: Prepare for the covariance update -----------------------------

  # omega
  mean_X_component_new <- matrix(log(omega_new[1:J-1]/omega_new[J]),
                                 nrow = 1)
  tilde_s_component_new <- t(mean_X_component_new)%*%mean_X_component_new
  covariance_component_new <- matrix(0,
                                     nrow = J-1,
                                     ncol = J-1)

  # alpha
  mean_X_alpha_new <- log(alpha_new)
  M_2_alpha_new <- mean_X_alpha_new^2
  variance_alpha_new <- 0

  # alpha_zero
  mean_X_alpha_zero_new <- log(alpha_zero_new)
  M_2_alpha_zero_new <- mean_X_alpha_zero_new^2
  variance_alpha_zero_new <- 0

  # q star
  mean_x_q_new <- lapply(1:J, function(j) matrix(log(q_star_1_J_new[j,1:(R-1)]/q_star_1_J_new[j,R]),
                                                 nrow = 1))
  tilde_s_q_new <- lapply(1:J, function(j) t(mean_x_q_new[[j]])%*%mean_x_q_new[[j]])
  covariance_q_new <- lapply(1:J, function(j) matrix(0, nrow = R-1, ncol = R-1))

  # gamma star
  X_mean_gamma_new <- log(gamma_1_J_star_new)
  M_2_gamma_new <-  X_mean_gamma_new^2
  variance_gamma_new <- rep(0,J)

  # alpha_h
  # mean_X_alpha_h_new <- matrix(log(alpha_h_new),
  #                              nrow = 1)
  # tilde_s_alpha_h_new <- t(mean_X_alpha_h_new)%*%mean_X_alpha_h_new
  # covariance_alpha_h_new <- matrix(0, nrow = R, ncol = R)
  mean_X_a_h_new <- log(a_h_new)
  M_2_a_h_new <- mean_X_a_h_new^2
  variance_a_h_new <- 0


  #----------------------- Step 4: Updates -----------------------------

  ##-- Set initial index
  output_index <- 0

  # Iteration starts with number 2
  for(iter in 2:number_iter){

    if(iter >= burn_in & iter%%thinning == 0){

      output_index <- output_index + 1
      update <- TRUE
    }else{

      update <- FALSE
    }

    if(iter %% iter_update == 0) print(paste("iter:", iter))
    print(paste0('Interations = ', iter))

    ##----------------------------------- Allocation variables ------------------------------------

    allocation_output <- allocation_variables_dirmult_mcmc(omega_J_M = omega_J_M_new,
                                                           q_star_1_J = q_star_1_J_new,
                                                           gamma_1_J_star = gamma_1_J_star_new,
                                                           Y = Y)

    Z_new <- allocation_output$Z
    allocation.prob <- allocation_output$allocation.prob

    ##-- Print table(Z) if TRUE
    if(isTRUE(print_Z)){
      for(m in 1:M) print(table(Z_new[[m]]))
    }
    
    print(paste0('J' = length(unique(unlist(Z_new)))))
    

    ##------------------------ Dataset specific component probabilities ------------------------------

    dataset_specfic_output <- dataset_specific_mcmc(Z = Z_new,
                                                    omega = omega_new,
                                                    alpha = alpha_new)

    omega_J_M_new <- dataset_specfic_output

    
    ##--------------------------- q_star ----------------------------------

    q_output_sim <- q_star_mcmc(Y = Y,
                                Z = Z_new,
                                q_star_1_J = q_star_1_J_new,
                                gamma_1_J_star = gamma_1_J_star_new,
                                alpha_h = alpha_h_new,
                                covariance = covariance_q_new,
                                mean_x = mean_x_q_new,
                                tilde_s = tilde_s_q_new,
                                iter_num = iter,
                                adaptive_prop = adaptive_prop)

    covariance_q_new <- q_output_sim$covariance_new
    mean_x_q_new <- q_output_sim$mean_x_new
    tilde_s_q_new <- q_output_sim$tilde_s_new
    q_star_1_J_new <- q_output_sim$q_star_1_J_new
    q_star_count <- q_output_sim$q_star_count + q_star_count
    acceptance_prob_list$q_star[iter-1] <- q_star_count/((iter-1)*J)


    ##-------------------------- gamma ---------------------------------

    gamma_output_sim <- gamma_mcmc(Y = Y,
                                   Z = Z_new,
                                   q_star_1_J = q_star_1_J_new,
                                   gamma_1_J_star = gamma_1_J_star_new,
                                   a_gamma = a_gamma,
                                   b_gamma = b_gamma,
                                   X_mean = X_mean_gamma_new,
                                   M_2 = M_2_gamma_new,
                                   variance = variance_gamma_new,
                                   iter_num = iter,
                                   adaptive_prop = adaptive_prop)

    gamma_1_J_star_new <- gamma_output_sim$gamma_1_J_star_new
    variance_gamma_new <- gamma_output_sim$variance_gamma_new
    M_2_gamma_new <- gamma_output_sim$M_2_gamma_new
    X_mean_gamma_new <- gamma_output_sim$X_mean_gamma_new
    gamma_star_count <- gamma_star_count + gamma_output_sim$gamma_count
    acceptance_prob_list$gamma_star[iter-1] <- gamma_star_count/((iter-1)*J)
    
    print(paste0('max gamma = ', max(gamma_1_J_star_new), 'min gamma = ', min(gamma_1_J_star_new)))

    
    
    ##----------------------------------- Split and merge moves ------------------------------------
    
    split_merge_output = split_merge_mcmc(omega_J_M = omega_J_M_new,
                                 q_1_J_star = q_star_1_J_new,
                                 gamma_1_J_star = gamma_1_J_star_new,
                                 Y = Y,
                                 Z = Z_new,
                                 alpha_h = alpha_h_new,
                                 a_gamma = a_gamma,
                                 b_gamma = b_gamma,
                                 alpha = alpha_new,
                                 omega = omega_new,
                                 alpha_zero = alpha_zero_new,
                                 split_accept = split_accept,
                                 merge_accept = merge_accept)

    gamma_1_J_star_new <- split_merge_output$gamma_1_J_star
    q_star_1_J_new <- split_merge_output$q_1_J_star
    omega_J_M_new <- split_merge_output$omega_J_M
    omega_new <- split_merge_output$omega
    Z_new <- split_merge_output$Z
    split_accept = split_merge_output$split_accept
    merge_accept = split_merge_output$merge_accept
    
    ##------------------------ Component probabilities -----------------------------------------------
    
    component_output <- component_probabilities_mcmc(omega = omega_new,
                                                     omega_J_M = omega_J_M_new,
                                                     alpha_zero = alpha_zero_new,
                                                     alpha = alpha_new,
                                                     covariance = covariance_component_new,
                                                     mean_x = mean_X_component_new,
                                                     tilde_s = tilde_s_component_new,
                                                     iter_num = iter,
                                                     adaptive_prop = adaptive_prop)
    
    omega_new <- component_output$omega_new
    tilde_s_component_new <- component_output$tilde_s_new
    mean_X_component_new <- component_output$mean_x_new
    covariance_component_new <- component_output$covariance_new
    omega_count <- omega_count + component_output$accept
    acceptance_prob_list$omega[iter-1] <- omega_count/(iter-1)
    
    ##-------------------------- Update alpha ----------------------------------
    
    alpha_output_sim <- alpha_mcmc(omega_J_M = omega_J_M_new,
                                   omega = omega_new,
                                   alpha = alpha_new,
                                   a_alpha = a_alpha,
                                   b_alpha = b_alpha,
                                   X_mean = mean_X_alpha_new,
                                   M_2 = M_2_alpha_new,
                                   variance = variance_alpha_new,
                                   iter_num = iter,
                                   adaptive_prop = adaptive_prop)
    
    alpha_new <- alpha_output_sim$alpha_new
    mean_X_alpha_new <- alpha_output_sim$X_mean_new
    M_2_alpha_new <- alpha_output_sim$M_2_new
    variance_alpha_new <- alpha_output_sim$variance_new
    alpha_count <- alpha_count + alpha_output_sim$accept
    acceptance_prob_list$alpha[iter-1] <- alpha_count/(iter-1)
    print(paste('alpha acceptance prob:', alpha_count/(iter-1), 'alpha=', alpha_new))
    
    ##----------------------- Update alpha_zero --------------------------------
    
    alpha_zero_output_sim <- alpha_zero_mcmc(omega = omega_new,
                                             alpha_zero = alpha_zero_new,
                                             a_alpha0 = a_alpha0,
                                             b_alpha0 = b_alpha0,
                                             X_mean = mean_X_alpha_zero_new,
                                             M_2 = M_2_alpha_zero_new,
                                             variance = variance_alpha_zero_new,
                                             iter_num = iter,
                                             adaptive_prop = adaptive_prop)
    
    alpha_zero_new <- alpha_zero_output_sim$alpha_zero_new
    mean_X_alpha_zero_new <- alpha_zero_output_sim$X_mean_new
    M_2_alpha_zero_new <- alpha_zero_output_sim$M_2_new
    variance_alpha_zero_new <- alpha_zero_output_sim$variance_new
    alpha_zero_count <- alpha_zero_count + alpha_zero_output_sim$accept
    acceptance_prob_list$alpha_zero[iter-1] <- alpha_zero_count/(iter-1)
    print(paste('alpha_zero acceptance prob:', alpha_zero_count/(iter-1), 'alpha_0=', alpha_zero_new))
    
    ##--------------------------- a_h ---------------------------------
    
    a_h_output_sim <- a_h_mcmc(a_h = a_h_new,
                               mean_X_a_h =  mean_X_a_h_new,
                               M_2_a_h = M_2_a_h_new,
                               q_star_1_J = q_star_1_J_new,
                               iter_num = iter,
                               a_a = a_a,
                               b_a = b_a,
                               mx = mx,
                               variance = variance_a_h_new,
                               adaptive_prop = adaptive_prop)
    
    a_h_new <- a_h_output_sim$a_h_new
    M_2_a_h_new <- a_h_output_sim$M_2_a_h_new
    mean_X_a_h_new <- a_h_output_sim$mean_X_a_h_new
    variance_a_h_new <- a_h_output_sim$variance_a_h_new
    a_h_count <- a_h_count + a_h_output_sim$accept
    acceptance_prob_list$a_h[iter-1] <- a_h_count/(iter-1)
    print(paste('a_h acceptance prob:', a_h_count/(iter-1),'a=', a_h_new))
    alpha_h_new = a_h_new*mx
    
    # alpha_h_output_sim <- alpha_h_mcmc(alpha_h = alpha_h_new,
    #                                    mean_X_alpha_h = mean_X_alpha_h_new,
    #                                    tilde_s_alpha_h = tilde_s_alpha_h_new,
    #                                    q_star_1_J = q_star_1_J_new,
    #                                    iter_num = iter,
    #                                    a_a = a_a,
    #                                    b_a = b_a,
    #                                    covariance = covariance_alpha_h_new,
    #                                    adaptive_prop = adaptive_prop)
    # 
    # alpha_h_new <- alpha_h_output_sim$alpha_h_new
    # tilde_s_alpha_h_new <- alpha_h_output_sim$tilde_s_alpha_h_new
    # mean_X_alpha_h_new <- alpha_h_output_sim$mean_X_alpha_h_new
    # covariance_alpha_h_new <- alpha_h_output_sim$covariance_alpha_h_new
    # alpha_h_count <- alpha_h_count + alpha_h_output_sim$accept
    # acceptance_prob_list$alpha_h[iter-1] <- alpha_h_count/(iter-1)
    # print(paste('alpha_h acceptance prob:', alpha_h_count/(iter-1)))
    
    

    #-------------------------- Step 5: Update simulated values ------------------------
    if(update == TRUE){

      Z_output[[output_index]] <- Z_new
      omega_J_M_output[[output_index]] <- omega_J_M_new
      omega_output[[output_index]] <- omega_new
      alpha_output[output_index] <- alpha_new
      alpha_zero_output[output_index] <- alpha_zero_new
      q_star_1_J_output[[output_index]] <- q_star_1_J_new
      allocation_prob_output[[output_index]] <- allocation.prob
      # alpha_h_output[[output_index]] <- alpha_h_new
      a_h_output[[output_index]] <- a_h_new
      gamma_star_1_J_output[[output_index]] <- gamma_1_J_star_new
    }

  }


  my_list <- list('acceptance_prob' = acceptance_prob_list,
                  'M' = M,
                  'C' = C,
                  'R' = R,
                  'allocation_probability' = allocation_prob_output,
                  'Z_output' = Z_output,
                  'omega_J_M_output' = omega_J_M_output,
                  'omega_output' = omega_output,
                  'alpha_output' = alpha_output,
                  'alpha_zero_output' = alpha_zero_output,
                  # 'alpha_h_output' = alpha_h_output,
                  'a_h_output' = a_h_output,
                  'q_star_1_J_output' = q_star_1_J_output,
                  'gamma_star_1_J_output' = gamma_star_1_J_output




  )

  return(my_list)

}




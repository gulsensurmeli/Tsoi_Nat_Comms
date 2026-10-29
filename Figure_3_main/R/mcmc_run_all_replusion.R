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
                         lb_gamma = 1,
                         a = 2,
                         tau = 0.4,
                         nu = 1/20,
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

  #--------------------------------------------------------------------------------

  # Acceptance probability
  acceptance_prob_list <- data.frame(omega = rep(0, number_iter-1),
                                     alpha = rep(0, number_iter-1),
                                     alpha_zero = rep(0, number_iter-1),
                                     q_star = rep(0, number_iter-1),
                                     gamma_star = rep(0, number_iter-1),
                                     split_merge = rep(0, number_iter-1))


  #----------------------- Step 2: Set starting values ----------------------------
  if (is.null(Z.init)){
    # Z_new <- lapply(1:M,
    #                 function(m) sample(1:J, size = C[m], replace = TRUE))
    # Initialize similar to k-means++
    data.j.prior <- matrix(do.call(cbind, Y), nrow = R)
    pdata.j.prior <- apply(data.j.prior, 2, function(x){x/sum(x)})
    Ycbind = do.call(cbind, Y)
    q_temp = matrix(0,1,R)
    c = sample(1:sum(C), size=1)
    q_temp[1,] = rdirichlet(1,alpha = Ycbind[,c] + 0.1 )
    ps = 1000
    q_prior = rdirichlet(ps, alpha = rep(1/J,J))
    for(j in 2:J){
      q_dist = sapply(1:sum(C), function(c){
        min(sapply(1:dim(q_temp)[1], function(l){
          c(dist(rbind(pdata.j.prior[,c],q_temp[l,]),method = "manhattan"))
        }))
      })
      temp = 10
      if(max(q_dist)<tau){
        q_dist_prior = sapply(1:ps, function(c){
          min(sapply(1:dim(q_temp)[1], function(l){
            c(dist(rbind(q_prior[c,],q_temp[l,]),method = "manhattan"))
          }))
        })
        myprob = exp(temp*q_dist_prior)*(q_dist_prior>tau)
        myprob = myprob/sum(myprob)
        c = sample(1:ps, size =1 , prob = myprob)
        q_temp = rbind(q_temp, q_prior[c,])
      }else{
        myprob = exp(temp*q_dist)*(q_dist>tau)
        myprob = myprob/sum(myprob)
        c =sample(1:sum(C), size=1, prob = myprob)
        q_temp = rbind(q_temp, rdirichlet(1,alpha = Ycbind[,c] + 0.1 ))
      }
      mindist = min(dist(q_temp, method = "manhattan"))
      if(mindist<tau){
        print(paste('min dist between initial q is ', mindist))
      }
    }
    # Allocate z based on closeness to q
    Z_new = lapply(1:M, function(m){
       Z.new.m = sapply(1:C[m], function(c){
         which.min(sapply(1:dim(q_temp)[1], function(l){
           c(dist(rbind(pdata.j.prior[,c],q_temp[l,]),method = "manhattan"))
         }))
       })
    })
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
  alpha_h <- a*mx

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
                              
                              mx2 <- matrix(0, nrow = 1, ncol = R)
                              # mx <- matrix(rdirichlet(1, alpha_h), nrow = 1, ncol = R)
                            }
                            
                            mx2
                            
                          })
  
  q_star_1_J_new <- do.call(rbind, q_star_1_J_new)
  
  mindist = min(dist(q_star_1_J_new, method = "manhattan"))
  if(mindist<tau){
    print(paste('min dist between q is ', mindist))
  }
  
  ind = c(1:J)[sapply(1:J, function(j){length(which(unlist(Z_new)==j)) == 0})]
  ind2 = c(1:J)[sapply(1:J, function(j){length(which(unlist(Z_new)==j)) > 0})]
  q_temp = q_star_1_J_new[ind2,]
  # Initialize empty ones that are far from occupied ones
  Ycbind = do.call(cbind, Y)
  for(j in ind){
    # Find the data point further from the current q
    q_dist = sapply(1:sum(C), function(c){
      min(sapply(1:dim(q_temp)[1], function(l){
        c(dist(rbind(pdata.j.prior[,c],q_temp[l,]),method = "manhattan"))
      }))
    })
    c =which.max(q_dist)
    q_star_1_J_new[j,] = rdirichlet(1,alpha = Ycbind[,c] + 0.1 )
    q_temp = rbind(q_temp,q_star_1_J_new[j,])
  }
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
                               mx2 = max(mx2, lb_gamma +1 )
                               
                             }else{
                               
                               mx2 <- rtgamma(1,shape = a_gamma, scale = 1/b_gamma, min = lb_gamma)
                             }
                             
                             mx2
                             
                           })
  
  gamma_1_J_star_new <- do.call(rbind, gamma_1_J_star_new)


  # Total count of acceptance
  omega_count <- 0
  alpha_count <- 0
  alpha_zero_count <- 0
  q_star_count <- 0
  gamma_star_count <- 0
  split_accept = 0
  merge_accept = 0

  #----------------------- Step 3: Prepare for the covariance update -----------------------------

  # omega
  # mean_X_component_new <- matrix(log(omega_new[1:J-1]/omega_new[J]),
  #                                nrow = 1)
  # tilde_s_component_new <- t(mean_X_component_new)%*%mean_X_component_new
  # covariance_component_new <- matrix(0,
  #                                    nrow = J-1,
  #                                    ncol = J-1)
  mean_X_component_new <- matrix(0,
                                nrow = 1,
                                ncol = J-1)
  tilde_s_component_new <- matrix(0,
                                  nrow = J-1,
                                  ncol = J-1)
  covariance_component_new <- matrix(0,
                                     nrow = J-1,
                                     ncol = J-1)
  scaling_omega = 0.1
  
  # alpha
  # mean_X_alpha_new <- log(alpha_new)
  # M_2_alpha_new <- mean_X_alpha_new^2
  # variance_alpha_new <- 0
  # mean_X_alpha_new <- 0
  # M_2_alpha_new <- 0
  # variance_alpha_new <- 0
  variance_alpha_new <- 0.1
  
  # alpha_zero
  # mean_X_alpha_zero_new <- log(alpha_zero_new)
  # M_2_alpha_zero_new <- mean_X_alpha_zero_new^2
  # variance_alpha_zero_new <- 0
  # mean_X_alpha_zero_new <- 0
  # M_2_alpha_zero_new <- 0
  # variance_alpha_zero_new <- 0
  variance_alpha_zero_new <- 0.1
  
  # q star
  # mean_x_q_new <- lapply(1:J, function(j) matrix(log(q_star_1_J_new[j,1:(R-1)]/q_star_1_J_new[j,R]),
  #                                                nrow = 1))
  # tilde_s_q_new <- lapply(1:J, function(j) t(mean_x_q_new[[j]])%*%mean_x_q_new[[j]])
  # covariance_q_new <- lapply(1:J, function(j) matrix(0, nrow = R-1, ncol = R-1))
  mean_x_q_new <- lapply(1:J, function(j) matrix(0,
                                                 nrow = 1,
                                                 ncol = R-1))
  tilde_s_q_new <- lapply(1:J, function(j) matrix(0,
                                                  nrow = R-1,
                                                  ncol = R-1))
  covariance_q_new <- lapply(1:J, function(j) matrix(0, nrow = R-1, ncol = R-1))
  q_star_iter_occupied = rep(1,J)
  scaling_q = rep(0.1, J)
  
  # gamma star
  # X_mean_gamma_new <- log(gamma_1_J_star_new)
  # M_2_gamma_new <-  X_mean_gamma_new^2
  # variance_gamma_new <- rep(0,J)
  # X_mean_gamma_new <- rep(0,J)
  # M_2_gamma_new <-  rep(0,J)
  # variance_gamma_new <- rep(0,J)
  variance_gamma_new <- rep(0.1,J)
  gamma_star_iter_occupied = rep(1,J)

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
                                alpha_h = alpha_h,
                                tau = tau,
                                nu = nu,
                                covariance = covariance_q_new,
                                mean_x = mean_x_q_new,
                                tilde_s = tilde_s_q_new,
                                scaling = scaling_q,
                                iter_num = q_star_iter_occupied,
                                adaptive_prop = adaptive_prop)

    covariance_q_new <- q_output_sim$covariance_new
    mean_x_q_new <- q_output_sim$mean_x_new
    tilde_s_q_new <- q_output_sim$tilde_s_new
    scaling_q <- q_output_sim$scaling
    q_star_1_J_new <- q_output_sim$q_star_1_J_new
    q_star_count <- q_output_sim$q_star_count + q_star_count
    acceptance_prob_list$q_star[iter-1] <- q_star_count/((iter-1)*J)
    q_star_iter_occupied = q_output_sim$iter_occupied

    mindist = min(dist(q_star_1_J_new, method = "manhattan"))
    if(mindist<tau){
      print(paste('min dist between q is ', mindist))
    }
    ##-------------------------- gamma ---------------------------------

    # gamma_output_sim <- gamma_mcmc(Y = Y,
    #                                Z = Z_new,
    #                                q_star_1_J = q_star_1_J_new,
    #                                gamma_1_J_star = gamma_1_J_star_new,
    #                                a_gamma = a_gamma,
    #                                b_gamma = b_gamma,
    #                                X_mean = X_mean_gamma_new,
    #                                M_2 = M_2_gamma_new,
    #                                variance = variance_gamma_new,
    #                                iter_num = iter,
    #                                adaptive_prop = adaptive_prop)
    # 
    # gamma_1_J_star_new <- gamma_output_sim$gamma_1_J_star_new
    # variance_gamma_new <- gamma_output_sim$variance_gamma_new
    # M_2_gamma_new <- gamma_output_sim$M_2_gamma_new
    # X_mean_gamma_new <- gamma_output_sim$X_mean_gamma_new
    # gamma_star_count <- gamma_star_count + gamma_output_sim$gamma_count
    # acceptance_prob_list$gamma_star[iter-1] <- gamma_star_count/((iter-1)*J)
    
    gamma_output_sim <- gamma_mcmc(Y = Y,
                                   Z = Z_new,
                                   q_star_1_J = q_star_1_J_new,
                                   gamma_1_J_star = gamma_1_J_star_new,
                                   a_gamma = a_gamma,
                                   b_gamma = b_gamma,
                                   lb_gamma = lb_gamma,
                                   variance = variance_gamma_new,
                                   iter_num = gamma_star_iter_occupied)

    gamma_1_J_star_new <- gamma_output_sim$gamma_1_J_star_new
    variance_gamma_new <- gamma_output_sim$variance_gamma_new
    gamma_star_count <- gamma_star_count + gamma_output_sim$gamma_count
    acceptance_prob_list$gamma_star[iter-1] <- gamma_star_count/((iter-1)*J)
    gamma_star_iter_occupied = gamma_output_sim$iter_occupied
    print(paste0('max gamma = ', max(gamma_1_J_star_new), 'min gamma = ', min(gamma_1_J_star_new)))

    
    
    ##----------------------------------- Split and merge moves ------------------------------------
    
    # split_merge_output = split_merge_mcmc(omega_J_M = omega_J_M_new,
    #                              q_1_J_star = q_star_1_J_new,
    #                              gamma_1_J_star = gamma_1_J_star_new,
    #                              Y = Y,
    #                              Z = Z_new,
    #                              alpha_h = alpha_h,
    #                              tau = tau,
    #                              nu = nu,
    #                              a_gamma = a_gamma,
    #                              b_gamma = b_gamma,
    #                              alpha = alpha_new,
    #                              omega = omega_new,
    #                              alpha_zero = alpha_zero_new,
    #                              split_accept = split_accept,
    #                              merge_accept = merge_accept)
    # 
    # gamma_1_J_star_new <- split_merge_output$gamma_1_J_star
    # q_star_1_J_new <- split_merge_output$q_1_J_star
    # omega_J_M_new <- split_merge_output$omega_J_M
    # omega_new <- split_merge_output$omega
    # Z_new <- split_merge_output$Z
    # split_accept = split_merge_output$split_accept
    # merge_accept = split_merge_output$merge_accept
    # acceptance_prob_list$split_merge[iter-1] <- (split_accept+ merge_accept)/(iter-1)
    
    ##------------------------ Component probabilities -----------------------------------------------
    
    component_output <- component_probabilities_mcmc(omega = omega_new,
                                                     omega_J_M = omega_J_M_new,
                                                     alpha_zero = alpha_zero_new,
                                                     alpha = alpha_new,
                                                     covariance = covariance_component_new,
                                                     mean_x = mean_X_component_new,
                                                     tilde_s = tilde_s_component_new,
                                                     scaling = scaling_omega,
                                                     iter_num = iter,
                                                     adaptive_prop = adaptive_prop)
    
    omega_new <- component_output$omega_new
    tilde_s_component_new <- component_output$tilde_s_new
    mean_X_component_new <- component_output$mean_x_new
    scaling_omega <- component_output$scaling
    covariance_component_new <- component_output$covariance_new
    omega_count <- omega_count + component_output$accept
    acceptance_prob_list$omega[iter-1] <- omega_count/(iter-1)
    
    ##-------------------------- Update alpha ----------------------------------
    
    # alpha_output_sim <- alpha_mcmc(omega_J_M = omega_J_M_new,
    #                                omega = omega_new,
    #                                alpha = alpha_new,
    #                                a_alpha = a_alpha,
    #                                b_alpha = b_alpha,
    #                                X_mean = mean_X_alpha_new,
    #                                M_2 = M_2_alpha_new,
    #                                variance = variance_alpha_new,
    #                                iter_num = iter,
    #                                adaptive_prop = adaptive_prop)
    # 
    # alpha_new <- alpha_output_sim$alpha_new
    # mean_X_alpha_new <- alpha_output_sim$X_mean_new
    # M_2_alpha_new <- alpha_output_sim$M_2_new
    # variance_alpha_new <- alpha_output_sim$variance_new
    # alpha_count <- alpha_count + alpha_output_sim$accept
    # acceptance_prob_list$alpha[iter-1] <- alpha_count/(iter-1)
    # print(paste('alpha acceptance prob:', alpha_count/(iter-1), 'alpha=', alpha_new))
    
    alpha_output_sim <- alpha_mcmc(omega_J_M = omega_J_M_new,
                                   omega = omega_new,
                                   alpha = alpha_new,
                                   a_alpha = a_alpha,
                                   b_alpha = b_alpha,
                                   variance = variance_alpha_new,
                                   iter_num = iter)
    
    alpha_new <- alpha_output_sim$alpha_new
    variance_alpha_new <- alpha_output_sim$variance_new
    alpha_count <- alpha_count + alpha_output_sim$accept
    acceptance_prob_list$alpha[iter-1] <- alpha_count/(iter-1)
    print(paste('alpha acceptance prob:', alpha_count/(iter-1), 'alpha=', alpha_new))
    
    ##----------------------- Update alpha_zero --------------------------------
    
    # alpha_zero_output_sim <- alpha_zero_mcmc(omega = omega_new,
    #                                          alpha_zero = alpha_zero_new,
    #                                          a_alpha0 = a_alpha0,
    #                                          b_alpha0 = b_alpha0,
    #                                          X_mean = mean_X_alpha_zero_new,
    #                                          M_2 = M_2_alpha_zero_new,
    #                                          variance = variance_alpha_zero_new,
    #                                          iter_num = iter,
    #                                          adaptive_prop = adaptive_prop)
    # 
    # alpha_zero_new <- alpha_zero_output_sim$alpha_zero_new
    # mean_X_alpha_zero_new <- alpha_zero_output_sim$X_mean_new
    # M_2_alpha_zero_new <- alpha_zero_output_sim$M_2_new
    # variance_alpha_zero_new <- alpha_zero_output_sim$variance_new
    # alpha_zero_count <- alpha_zero_count + alpha_zero_output_sim$accept
    # acceptance_prob_list$alpha_zero[iter-1] <- alpha_zero_count/(iter-1)
    # print(paste('alpha_zero acceptance prob:', alpha_zero_count/(iter-1), 'alpha_0=', alpha_zero_new))
    
    alpha_zero_output_sim <- alpha_zero_mcmc(omega = omega_new,
                                             alpha_zero = alpha_zero_new,
                                             a_alpha0 = a_alpha0,
                                             b_alpha0 = b_alpha0,
                                             variance = variance_alpha_zero_new,
                                             iter_num = iter)
    
    alpha_zero_new <- alpha_zero_output_sim$alpha_zero_new
    variance_alpha_zero_new <- alpha_zero_output_sim$variance_new
    alpha_zero_count <- alpha_zero_count + alpha_zero_output_sim$accept
    acceptance_prob_list$alpha_zero[iter-1] <- alpha_zero_count/(iter-1)
    print(paste('alpha_zero acceptance prob:', alpha_zero_count/(iter-1), 'alpha_0=', alpha_zero_new))
    
    

    #-------------------------- Step 5: Update simulated values ------------------------
    if(update == TRUE){

      Z_output[[output_index]] <- Z_new
      omega_J_M_output[[output_index]] <- omega_J_M_new
      omega_output[[output_index]] <- omega_new
      alpha_output[output_index] <- alpha_new
      alpha_zero_output[output_index] <- alpha_zero_new
      q_star_1_J_output[[output_index]] <- q_star_1_J_new
      allocation_prob_output[[output_index]] <- allocation.prob
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
                  'q_star_1_J_output' = q_star_1_J_output,
                  'gamma_star_1_J_output' = gamma_star_1_J_output




  )

  return(my_list)

}




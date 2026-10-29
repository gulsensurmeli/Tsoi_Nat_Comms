
#' Obtain posterior samples of the mouse specific component probabilities
#'

mcmc_run_omega_JM <- function(mcmc_run_all_output,
                              mcmc_run_post_output,
                              J,
                              thinning,
                              burn_in,
                              number_iter,
                              axis.text.font = 10,
                              axis.title.font = 10,
                              adaptive_prop = 0.1){
  
  
  
  labels <- mcmc_run_post_output$cluster.label.summary
  Z <- mcmc_run_post_output$Z
  
  # Posterior mean
  mean_alpha <- mean(mcmc_run_all_output$alpha_output)
  mean_alpha_zero <- mean(mcmc_run_all_output$alpha_zero_output)
  
  
  # Number of occupied components
#  J <- max(unlist(Z))
  
  # M
  M <- length(Z)
  
  # Initialize omega
  omega_hat = function(z=z,J=J){
    omega = rep(1,J)
    tbz = table(z)
    omega[as.numeric(names(tbz))] = omega[as.numeric(names(tbz))] + tbz
    return(omega/sum(omega))
  }
  
  omega_new <- omega_hat(unlist(Z), J=J)
#  omega_new <- rep(1/J, J)
  
  mean_X_component_new <- matrix(log(omega_new[1:J-1]/omega_new[J]),
                                 nrow = 1)
  tilde_s_component_new <- t(mean_X_component_new)%*%mean_X_component_new
  covariance_component_new <- matrix(0,
                                     nrow = J-1,
                                     ncol = J-1)
  
  # Set initial index
  output_index <- 0
  
  omega_J_M_output <- NULL
  omega_output <- NULL
  
  # Iteration starts with number 2
  for(iter in 2:number_iter){
    
    if(iter >= burn_in & iter%%thinning == 0){
      
      output_index <- output_index + 1
      update <- TRUE
    }else{
      
      update <- FALSE
    }
    
    if(iter %% 100 == 0) print(paste("iter:", iter))
    
    ##--------------------- omega_J_M -----------------------------
    
    dataset_specfic_output <- dataset_specific_mcmc(Z = Z,
                                                    omega = omega_new,
                                                    alpha = mean_alpha)
    
    omega_J_M_new <- dataset_specfic_output
    
    ##------------------------ omega ------------------------------
    
    component_output <- component_probabilities_mcmc(omega = omega_new,
                                                     omega_J_M = omega_J_M_new,
                                                     alpha_zero = mean_alpha_zero,
                                                     alpha = mean_alpha,
                                                     covariance = covariance_component_new,
                                                     mean_x = mean_X_component_new,
                                                     tilde_s = tilde_s_component_new,
                                                     iter_num = iter,
                                                     adaptive_prop = adaptive_prop)
    
    omega_new <- component_output$omega_new
    tilde_s_component_new <- component_output$tilde_s_new
    mean_X_component_new <- component_output$mean_x_new
    covariance_component_new <- component_output$covariance_new
    
    ##-- Update
    if(isTRUE(update)){
      
      omega_J_M_output[[output_index]] <- omega_J_M_new
      omega_output[[output_index]] <- omega_new
    }
  }
  
  ##------------------- Posterior mean -----------------------------
  
  mean_omega_J_M <- matrix(0, nrow = J, ncol = M)
  
  for(j in 1:J){
    for(m in 1:M){
      
      mean_omega_J_M[j,m] <- mean(sapply(1:length(omega_J_M_output), function(i) omega_J_M_output[[i]][j,m]))
    }
  }
  
  omega_J_matrix <- do.call(rbind, omega_output)
  mean_omega_J <- colSums(omega_J_matrix)/nrow(omega_J_matrix)
  
  ##------------------- 90 percent credible interval ---------------
  
  lower_omega_J_M <- matrix(0, nrow = J, ncol = M)
  upper_omega_J_M <- matrix(0, nrow = J, ncol = M)
  
  for(j in 1:J){
    for(m in 1:M){
      
      lower_omega_J_M[j,m] <- quantile(sapply(1:length(omega_J_M_output), function(i) omega_J_M_output[[i]][j,m]), 0.05)
      upper_omega_J_M[j,m] <- quantile(sapply(1:length(omega_J_M_output), function(i) omega_J_M_output[[i]][j,m]), 0.95)
    }
  }
  
  
  ##--------------- Summary table -------------------
  
  loop.result <- lapply(1:M, function(m){
    
    data.frame(mouse = factor(m, levels = 1:M),
               mean_estimate = mean_omega_J_M[,m],
               lower_estimate = lower_omega_J_M[,m],
               upper_estimate = upper_omega_J_M[,m],
               cluster = paste('cluster', 1:J))
  })
  
  omega_JM_estimate_df <- do.call(rbind, loop.result)
  
  
  
  
  
  #-------------------------------------- Plot for all clusters ---------------------------------------
  
  omega_JM_plot <- ggplot(data = omega_JM_estimate_df,
                          mapping = aes(x = mouse,
                                        y = mean_estimate,
                                        fill = mouse))+
    
    geom_bar(
      stat = 'identity',
      position = position_dodge(),
      alpha = 0.5)+
    geom_errorbar(mapping = aes(x = mouse, ymin = lower_estimate, ymax = upper_estimate), width=.2,
                  position=position_dodge(.9))+
    theme_bw()+
    ylab('Mouse-specific component probability')+
    facet_wrap(~factor(cluster,
                       levels = paste('cluster', 1:J)),
               scale = 'free')+
    theme(legend.position = "none") +
    theme(axis.text=element_text(size=axis.text.font),
          axis.title=element_text(size=axis.title.font),
          plot.title = element_text(size=12))
  
  omega_JM_estimate_df$mouse <- paste('mouse', omega_JM_estimate_df$mouse)
  
  
  
  #--------------------------------- Return items -----------------------------------
  
  return(list('omega_J_M_output' = omega_J_M_output,
              'mean_omega_J_M' = mean_omega_J_M,
              'lower_omega_J_M' = lower_omega_J_M,
              'upper_omega_J_M' = upper_omega_J_M,
              
              'omega_JM_estimate_df' = omega_JM_estimate_df,
              
              'omega_JM_plot' = omega_JM_plot,
              'mean_omega_J' = mean_omega_J,
              'omega_output' = omega_output)
         )
  
}

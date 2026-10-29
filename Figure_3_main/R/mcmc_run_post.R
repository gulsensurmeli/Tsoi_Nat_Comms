
#' MCMC sampling of gamma and projection strength
#'

mcmc_run_post <- function(mcmc_run_all_output,
                          Z,
                          thinning = 5,
                          burn_in,
                          number_iter,
                          Y,
                          regions.name = NULL,
                          a_gamma,
                          b_gamma,
                          adaptive_prop = 0.1,
                          axis.title.size = 10,
                          axis.text.size = 10,
                          plot.title.size = 10,
                          cred.int = 0.9,
                          epsilon =0.01,
                          prob_thresh = 0.5){
  
  
  
  # Dimensions
  J <- length(unique(unlist(Z)))
  R <- mcmc_run_all_output$R
  M <- length(Y)
  C <- sapply(1:M,
              function(m) ncol(Y[[m]]))
  
  if(is.null(regions.name)){
    
    regions.name <- paste('region', 1:R)
  }
  
  # Hyperparameters for q
  # Compute empirical mean of projection strength
  # data.j.prior <- matrix(do.call(cbind, Y), nrow = R)
  # pdata.j.prior <- apply(data.j.prior, 2, function(x){x/sum(x)})
  # data.j.mean <- rowMeans(pdata.j.prior)
  # mx <- matrix(data.j.mean/sum(data.j.mean), nrow = 1)
  # a_h = mean(unlist(mcmc_run_all_output$a_h_output))
  # mean_alpha_h <- a_h*mx
  
  # # Posterior mean of alpha h
   alpha_h_all <- do.call(rbind,
                          mcmc_run_all_output$alpha_h_output)
   mean_alpha_h <- colSums(alpha_h_all)/nrow(alpha_h_all)
  
  # Initial value of q star
  q_star_1_J_new <- lapply(1:J,
                           function(j){
                             
                             if(length(which(unlist(Z)==j)) != 0){
                               
                               data.j.prior <- matrix(do.call(cbind, Y)[,which(unlist(Z)==j)],
                                                      nrow = R)
                               
                               pdata.j.prior <- apply(data.j.prior, 2, function(x){x/sum(x)})
                               data.j.mean <- rowMeans(pdata.j.prior)
                               
                               mx <- matrix(data.j.mean/sum(data.j.mean),
                                            nrow = 1)
                               
                             }else{
                               
                               mx <- matrix(rdirichlet(1, mean_alpha_h), nrow = 1, ncol = R)
                             }
                             
                             mx
                             
                           })
  
  q_star_1_J_new <- do.call(rbind, q_star_1_J_new)
  q_star_1_J_new <- ifelse(q_star_1_J_new == 0, .Machine$double.eps,  q_star_1_J_new)
  q_star_1_J_new <-  q_star_1_J_new/rowSums(q_star_1_J_new)
  
  # Initial covariance structure
  mean_x_q_new <- lapply(1:J, function(j) matrix(log(q_star_1_J_new[j,1:(R-1)]/q_star_1_J_new[j,R]),
                                                 nrow = 1))
  tilde_s_q_new <- lapply(1:J, function(j) t(mean_x_q_new[[j]])%*%mean_x_q_new[[j]])
  covariance_q_new <- lapply(1:J, function(j) matrix(0, nrow = R, ncol = R))
  
  
  # Set gamma empirically based on Z_new
  gamma_1_J_star_new  <- lapply(1:J,
                                function(j){
                                  
                                  if(length(which(unlist(Z)==j)) > 1){
                                    
                                    data.j.prior <- matrix(do.call(cbind, Y)[,which(unlist(Z)==j)],
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
  
  
  # gamma star covariance structure
  X_mean_gamma_new <- log(gamma_1_J_star_new)
  M_2_gamma_new <- rep(0,J)
  variance_gamma_new <- rep(0,J)
  
  
  # projection probabilities
  q_star_1_J_output <- NULL
  
  # gamma star
  gamma_star_1_J_output <- NULL
  
  # Set initial index
  output_index <- 0
  
  # Iteration starts with number 2
  for(iter in 2:number_iter){
    
    if(iter >= burn_in & iter%%thinning == 0){
      
      output_index <- output_index + 1
      update <- TRUE
    }else{
      
      update <- FALSE
    }
    
    ##-- Print for every 100th iteration
    if(iter %% 100 == 0) print(paste("iter:", iter))
    
    # Posterior sample of gamma
    gamma_output_sim <- gamma_mcmc(Y = Y,
                                   Z = Z,
                                   q_star_1_J = q_star_1_J_new,
                                   gamma_1_J_star = gamma_1_J_star_new,
                                   a_gamma = a_gamma,
                                   b_gamma = b_gamma,
                                   X_mean = X_mean_gamma_new,
                                   M_2 = M_2_gamma_new,
                                   variance = variance_gamma_new,
                                   iter_num = iter,
                                   adaptive_prop = adaptive_prop)
    
    # Update covariance structure of gamma
    gamma_1_J_star_new <- gamma_output_sim$gamma_1_J_star_new
    variance_gamma_new <- gamma_output_sim$variance_gamma_new
    M_2_gamma_new <- gamma_output_sim$M_2_gamma_new
    X_mean_gamma_new <- gamma_output_sim$X_mean_gamma_new
    
    # Posterior sample of projection probability
    q_output_sim <- q_star_mcmc(Y = Y,
                                Z = Z,
                                q_star_1_J = q_star_1_J_new,
                                gamma_1_J_star = gamma_1_J_star_new,
                                alpha_h = mean_alpha_h,
                                covariance = covariance_q_new,
                                mean_x = mean_x_q_new,
                                tilde_s = tilde_s_q_new,
                                iter_num = iter,
                                adaptive_prop = adaptive_prop)
    
    # Update covariance structure of q
    covariance_q_new <- q_output_sim$covariance_new
    mean_x_q_new <- q_output_sim$mean_x_new
    tilde_s_q_new <- q_output_sim$tilde_s_new
    q_star_1_J_new <- q_output_sim$q_star_1_J_new
    
    # If update is TRUE
    if(isTRUE(update)){
      
      q_star_1_J_output[[output_index]] <- q_star_1_J_new
      gamma_star_1_J_output[[output_index]] <- gamma_1_J_star_new
    }
    
  }
  
  ##-- Compute posterior mean
  proj_prob_mean <- Reduce("+", q_star_1_J_output)/length(q_star_1_J_output)
  gamma_mean <- Reduce("+", gamma_star_1_J_output)/length(gamma_star_1_J_output)
  
  ##-- Credible intervals
  proj_prob_lower <- matrix(0, nrow = J, ncol = R)
  proj_prob_upper <- matrix(0, nrow = J, ncol = R)
  proj_prob_med <- matrix(0, nrow = J, ncol = R)
  
  for(j in 1:J){
    for(r in 1:R){
      
      proj_prob_lower[j,r] <- quantile(sapply(1:length(q_star_1_J_output), function(i) q_star_1_J_output[[i]][j,r]), (1-cred.int)/2)
      proj_prob_upper[j,r] <- quantile(sapply(1:length(q_star_1_J_output), function(i) q_star_1_J_output[[i]][j,r]), 1-(1-cred.int)/2)
      proj_prob_med[j,r] <- quantile(sapply(1:length(q_star_1_J_output), function(i) q_star_1_J_output[[i]][j,r]), 0.5)
      
    }
  }
  
  
  
  #----------------------- Total number of projecting regions ---------------------------------
  
  neuron_projection_df <- NULL
  neuron_projection_df$q_star_1_J_output <- q_star_1_J_output
  neuron_projection_df$gamma_star_1_J_output <- gamma_star_1_J_output
  neuron_projection_df$Z <- Z
  neuron_projection_df$J <- J
  
  q_tilde_001 <- q_tilde(neuron_projection_df = neuron_projection_df,
                         epsilon = epsilon)
  q_tilde_nonzero <- q_tilde2(neuron_projection_df = neuron_projection_df,
                         n = 100)
  
  unicast.index <- (1:J)[apply(q_tilde_001,
                               1,
                               function(x) length(which(x>prob_thresh)) == 1)]
  
  
  broadcast.index <- (1:J)[apply(q_tilde_001,
                                 1,
                                 function(x) length(which(x>prob_thresh)) == R)]
  
  bicast.and.more.index <- (1:J)[-unique(c(unicast.index,
                                           broadcast.index))]
  
  
  # Cluster labels
  cluster.label.summary <- data.frame(cluster = 1:J,
                                      
                                      unicast = ifelse((1:J) %in% unicast.index, 'Yes', 'No'),
                                      bicast.and.more = ifelse((1:J) %in% bicast.and.more.index, 'Yes', 'No'),
                                      broadcast = ifelse((1:J) %in% broadcast.index, 'Yes', 'No')
                                      
  )
  
  # Relabel clusters based on the strength of projection
  strong.proj <- apply(proj_prob_mean,
                       1,
                       function(x) which(x == max(x)))
  
  
  df0 <- lapply(1:R,
                function(r){
                  
                  # Clusters with the strongest projection to region r
                  strong.proj.r <- which(strong.proj==r)
                  
                  
                  if(length(strong.proj.r) != 0){
                    
                    qj_estimate_r <- matrix(proj_prob_mean[strong.proj.r,],
                                            nrow = length(strong.proj.r))
                    
                    output <- strong.proj.r[order(-qj_estimate_r[,r])]
                  }else{
                    
                    output <- NULL
                    
                  }
                  }
                
                )
    
    
    # Reorder original label
    original.label <- unlist(df0)
    
    
    # Relabel the summary table
    cluster.label.summary <- cluster.label.summary[original.label,]
    cluster.label.summary$cluster <- 1:J
    
    
    # Change to the new labeling system
    q_star_1_J_output <- lapply(1:length(q_star_1_J_output),
                                function(t) q_star_1_J_output[[t]][original.label,])
    
    gamma_star_1_J_output <- lapply(1:length(gamma_star_1_J_output),
                                    function(t) gamma_star_1_J_output[[t]][original.label])
    
    proj_prob_mean <- proj_prob_mean[original.label,]
    proj_prob_lower <- proj_prob_lower[original.label,]
    proj_prob_upper <- proj_prob_upper[original.label,]
    proj_prob_med <- proj_prob_med[original.label,]
    
    gamma_mean <- gamma_mean[original.label]
    
    
    
    Z <- lapply(1:M,
                function(m){
                  
                  sapply(1:C[m], function(c) which(original.label == Z[[m]][c]))
                })
    
    
    # Summary of projection probabilities
    loop.result <- lapply(1:J, function(j){
      
      data.frame(cluster = paste('cluster', j),
                 region = regions.name,
                 projection.mean = proj_prob_mean[j,],
                 projection.lower = proj_prob_lower[j,],
                 projection.upper = proj_prob_upper[j,],
                 projection.med = proj_prob_med[j,])
    })
    
    estimated.projection.df <- do.call(rbind, loop.result)
    
    
    
    
  
  
  
  #---------------------------------------- Plot of projection strengths ------------------------------------
  
  
  
  cluster_label_vec <- sapply(1:J,
                              function(j){
                                
                                colnames(cluster.label.summary)[2:4][which(cluster.label.summary[j,2:4]=='Yes')[1]]
                              })
  
  # Label the number of projecting regions
  estimated.projection.df$class <- rep(cluster_label_vec, each = R)
  estimated.projection.df$class <- factor(estimated.projection.df$class,
                                          levels = colnames(cluster.label.summary)[2:4])
  
  estimated.projection.df$cluster <- factor(estimated.projection.df$cluster,
                                            levels = paste('cluster', 1:J))
  
  # All clusters
  estimated.pp.plot <- ggplot2::ggplot(estimated.projection.df,
                                       mapping = aes(x = factor(region, levels = regions.name),
                                                     y = projection.med,
                                                     group = cluster,
                                                     color = class)) +
    geom_line()+
    geom_point()+
    geom_errorbar(aes(ymin = projection.lower,
                      ymax = projection.upper),
                  width = 0.1)+
    theme_bw()+
    ylim(c(0,1))+
    ylab('projection strength')+
    xlab('region')+
    facet_wrap(~cluster)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    theme(axis.text=element_text(size=axis.text.size),
          axis.title=element_text(size=axis.title.size),
          plot.title = element_text(size=plot.title.size))
  
  
  
  
  #--------------------------------------- Plot of q_tilde ----------------------------------
  
  q_tilde_001 <- q_tilde_001[original.label,]
  q_tilde_nonzero <- q_tilde_nonzero[original.label,]
  
  colnames(q_tilde_001) <- regions.name
  rownames(q_tilde_001) <- paste('cluster', 1:J)
  colnames(q_tilde_nonzero) <- regions.name
  rownames(q_tilde_nonzero) <- paste('cluster', 1:J)
  
  
  q_tilde_001_table <- data.frame(region = rep(regions.name,
                                               each = J),
                                  cluster = rep(paste('cluster', 1:J),
                                                R),
                                  probability = as.vector(q_tilde_001))
  
  q_tilde_001_table$region <- factor(q_tilde_001_table$region,
                                     levels = regions.name)
  
  q_tilde_001_table$cluster <- factor(q_tilde_001_table$cluster,
                                      levels = paste('cluster', 1:J))
  
  q_tilde_plot <- q_tilde_001_table %>%
    ggplot(mapping = aes(x = cluster,
                         y = region,
                         fill = probability))+
    geom_tile()+
    theme_bw()+
    scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                         values = scales::rescale(c(0, 0.25, 1)),
                         limits = c(0,1))+
    xlab('cluster')+
    ylab('regions')+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  ##-- Output values
  return(list('J' = J,
              'q_star_1_J_output' = q_star_1_J_output,
              'proj_prob_mean' = proj_prob_mean,
              'proj_prob_med' = proj_prob_med,
              'proj_prob_lower' = proj_prob_lower,
              'proj_prob_upper' = proj_prob_upper,
              'Z' = Z,
              'gamma_mean' = gamma_mean,
              'gamma_star_1_J_output' = gamma_star_1_J_output,
              'estimated.projection.df' = estimated.projection.df,
              'cluster.label.summary' = cluster.label.summary,
              'estimated.pp.plot' = estimated.pp.plot,
              'q_tilde_001' = q_tilde_001,
              'q_tilde_nonzero' = q_tilde_nonzero,
              'q_tilde_plot' = q_tilde_plot)
  )
  }


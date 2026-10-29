#' Posterior predictive check with multiple replicated datasets
#'

ppc_f <- function(mcmc_run_all_output,
                  Y,
                  N,
                  regions.name = NULL){



  length_per_chain <- length(mcmc_run_all_output$Z_output)

  C <- mcmc_run_all_output$C
  R <- mcmc_run_all_output$R
  M <- mcmc_run_all_output$M

  if(is.null(regions.name)){

    regions.name <- paste('region', 1:R)
  }


  # Sum of counts for each neuron cell
  N_CM <- lapply(1:M,
                 function(m) colSums(Y[[m]]))

  # theta
  target.index <- sample(length_per_chain,
                         size = N,
                         replace = FALSE)


  theta <- lapply(1:N,
                  function(n){

                    list(q_j_star = mcmc_run_all_output$q_star_1_J_output[[target.index[n]]],
                         gamma_j_star = mcmc_run_all_output$gamma_star_1_J_output[[target.index[n]]],
                         Z = mcmc_run_all_output$Z_output[[target.index[n]]])
                  })

  # Simulate replicated data Y
  Y_rep <- lapply(1:N,
                  function(n){

                    Y_rep_n <- lapply(1:M,
                                      function(m){

                                        do.call(cbind, lapply(1:C[m],
                                                              function(c){

                                                                prob.c <- rdirichlet(n = 1,
                                                                                     theta[[n]]$q_j_star[theta[[n]]$Z[[m]][c],]*theta[[n]]$gamma_j_star[theta[[n]]$Z[[m]][c]])

                                                                rmultinom(n = 1, size = N_CM[[m]][c], prob = prob.c)
                                                              }))
                                      })

                    Y_rep_n
                  })



  #------------------------------------------ Some statistics --------------------------------------

  Y_prop <- lapply(1:M,
                   function(m){

                     matrix(as.vector(Y[[m]])/rep(colSums(Y[[m]]), each = R),
                            nrow = R)
                   })

  Y_rep_prop <- lapply(1:N,
                       function(n){

                         lapply(1:M,
                                function(m){

                                  matrix(as.vector(Y_rep[[n]][[m]])/rep(colSums(Y_rep[[n]][[m]]), each = R),
                                         nrow = R)
                                })

                       })

  Y_prop_df <- data.frame(Region = rep(regions.name, sum(C)),
                          data = 'observed',
                          value = as.vector(do.call(cbind,Y_prop)))

  Y_rep_prop_df <- do.call(rbind, lapply(1:N,
                                         function(n){

                                           data.frame(Region = rep(regions.name, sum(C)),
                                                      data = paste('replicated', n),
                                                      value = as.vector(do.call(cbind,Y_rep_prop[[n]])))

                                         }))

  # Convert into factor to keep brain region ordering
  Y_prop_df$Region <- factor(Y_prop_df$Region,
                             levels = regions.name)

  Y_rep_prop_df$Region <- factor(Y_rep_prop_df$Region,
                                 levels = regions.name)

  # Box plot
  non.zero.plot <- rbind(Y_prop_df,
                         Y_rep_prop_df) %>%
    filter(value != 0) %>%
    ggplot(mapping = aes(x = Region,
                         y = value,
                         fill = data))+

    geom_boxplot(outlier.shape = NA)+
    theme_bw()+
    ylab('projection strength')

  zero.prop <- rbind(Y_prop_df,
                     Y_rep_prop_df) %>%
    filter(value == 0) %>%
    group_by(Region, data) %>%
    summarise(N = n())

  zero.plot <- ggplot(zero.prop, aes(x = Region,
                                     y = N,
                                     fill = data))+
    geom_bar(stat = 'identity',
             position = position_dodge())+
    theme_bw()+
    ylab('Number of zeros')


  # Return replicated Y
  return(list('Y_prop' = Y_prop,
              'Y_rep_prop' = Y_rep_prop,
              'theta' = theta,
              'non.zero.plot' = non.zero.plot,
              'zero.plot' = zero.plot))


}

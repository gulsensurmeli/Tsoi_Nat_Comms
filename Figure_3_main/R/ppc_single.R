#' Posterior predictive check with single replicated data
#'

ppc_single_f <- function(mcmc_run_all_output,
                         Y,
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
                         size = 1,
                         replace = FALSE)


  theta <- list(q_j_star = mcmc_run_all_output$q_star_1_J_output[[target.index]],
                gamma_j_star = mcmc_run_all_output$gamma_star_1_J_output[[target.index]],
                Z = mcmc_run_all_output$Z_output[[target.index]])

  # Simulate replicated data Y

  Y_rep <- lapply(1:M,
                  function(m){

                    do.call(cbind, lapply(1:C[m],
                                          function(c){

                                            prob.c <- rdirichlet(n = 1,
                                                                 theta$q_j_star[theta$Z[[m]][c],]*theta$gamma_j_star[theta$Z[[m]][c]])

                                            rmultinom(n = 1, size = N_CM[[m]][c], prob = prob.c)
                                          }))
                  })




  #------------------------------------------ Some statistics --------------------------------------

  Y_prop <- lapply(1:M,
                   function(m){

                     matrix(as.vector(Y[[m]])/rep(colSums(Y[[m]]), each = R),
                            nrow = R)
                   })

  Y_rep_prop <- lapply(1:M,
                       function(m){

                         matrix(as.vector(Y_rep[[m]])/rep(colSums(Y_rep[[m]]), each = R),
                                nrow = R)
                       })

  Y_prop_df <- data.frame(Region = rep(regions.name, sum(C)),
                          data = 'observed',
                          mouse = rep(paste('brain', 1:M), c(R*C)),
                          value = as.vector(do.call(cbind,Y_prop)))

  Y_rep_prop_df <- data.frame(Region = rep(regions.name, sum(C)),
                              data = 'replicated',
                              mouse = rep(paste('brain', 1:M), c(R*C)),
                              value = as.vector(do.call(cbind,Y_rep_prop)))

  # For each of the mouse, draw histograms
  plot.out <- lapply(1:M,
                     function(m){


                       plot.out.m <- rbind(Y_prop_df,
                                           Y_rep_prop_df) %>%
                         filter(mouse == paste('brain', m)) %>%
                         ggplot(mapping = aes(x = value,
                                              fill = data))+
                         geom_histogram(position = 'dodge')+
                         theme_bw()+
                         facet_wrap(~Region, nrow = 2, scales = 'free')+
                         ggtitle(paste('Brain', m))+
                         xlab('projection strength')

                       plot.out.m

                     })

  return(plot.out)


}

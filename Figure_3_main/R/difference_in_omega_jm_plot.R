#' Function to plot histograms to compare mouse specific component probabilities
#'

difference_in_omega_jm_plot <- function(difference_in_omega_jm_output,
                                        mcmc_run_omega_output,
                                        N){


  df0 <- difference_in_omega_jm_output$significant_obs
  stopifnot(nrow(df0) >= N)


  omega_jm_trace <- mcmc_run_omega_output$omega_J_M_output
  trace.length <- length(omega_jm_trace)

  df_new <- lapply(sample(1:nrow(df0),
                          size = N,
                          replace = FALSE),

                   function(i){


                     data.j.trace <- lapply(1:trace.length,
                                            function(s) omega_jm_trace[[s]][df0$cluster[i],c(df0$data1[i],
                                                                                             df0$data2[i])])

                     data.j.trace <- do.call(rbind,
                                             data.j.trace)

                     data.frame(cluster = df0$cluster[i],
                                data = c(rep(df0$data1[i], trace.length),
                                         rep(df0$data2[i], trace.length)),
                                omega_jm = as.vector(data.j.trace),
                                group = paste('cluster', df0$cluster[i],
                                              'mouse', df0$data1[i],
                                              'and', df0$data2[i]))
                   })

  df_new <- do.call(rbind,
                    df_new)

  df_new$mouse <- as.factor(df_new$data)

  df_new %>%
    ggplot(mapping = aes(x = omega_jm,
                         color = mouse))+
    geom_histogram(fill = 'white',
                   alpha = 0.5,
                   position = 'dodge')+
    facet_wrap(~group,
               scales = "free")+
    theme_bw()




}

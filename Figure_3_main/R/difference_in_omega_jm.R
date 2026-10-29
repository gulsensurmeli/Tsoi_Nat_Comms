#' Difference in mouse specific component probabilities
#'

difference_in_omega_jm <- function(mcmc_run_omega_output){

  omega_jm_trace <- mcmc_run_omega_output$omega_J_M_output
  trace.length <- length(omega_jm_trace)

  J <- nrow(omega_jm_trace[[1]])

  # Number of mouse
  M <- ncol(omega_jm_trace[[1]])

  # Unique set of cluster
  unique.j <- 1:J

  # All possible combinations of d
  combo <- expand.grid(1:M, 1:M)
  combo <- combo[!duplicated(t(apply(combo, 1, sort))),]
  combo <- combo[!t(apply(combo, 1, function(x) x[1] == x[2])),]

  # data to save
  df0 <- lapply(unique.j,

                function(j){

                  df00 <- lapply(1:nrow(combo),

                                 function(i){

                                   data1.index <- combo[i,1]
                                   data2.index <- combo[i,2]

                                   data.j.trace <- lapply(1:trace.length,
                                                          function(s) omega_jm_trace[[s]][j,c(data1.index,
                                                                                              data2.index)])

                                   data.j.trace <- do.call(rbind,
                                                           data.j.trace)

                                   prob <- length(which(data.j.trace[,1] > data.j.trace[,2]))/trace.length


                                   data.frame(cluster = j,
                                              data1 = data1.index,
                                              data2 = data2.index,
                                              prob = prob)
                                 })

                  df00 <- do.call(rbind,
                                  df00)
                })

  df0 <- do.call(rbind,
                 df0)

  df0_all <- df0

  #--------------- Plotting
  df1_all <- data.frame(cluster = factor(paste('cluster', df0_all$cluster),
                                         levels = paste('cluster', 1:J)),
                        mouse.pair = paste('(', df0_all$data1, ',', df0_all$data2, ')', sep = ''),
                        probability = df0_all$prob)


  probability_plot <- df1_all %>%
    ggplot(mapping = aes(x = cluster,
                         y = mouse.pair,
                         fill = probability))+
    geom_tile()+
    theme_bw()+
    scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                         values = scales::rescale(c(0, 0.25, 1)),
                         limits = c(0,1))+
    xlab('cluster')+
    ylab('mouse pairs')+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  # Select the cluster and data with different p_jd
  df0_significant <- df0[which(df0$prob > 0.95 | df0$prob < 0.05),]

  # Plot for each data
  plot.for.each.data <- lapply(1:M,
                               function(m){
                                 
                                 df0_all_m <- df0_all %>%
                                   filter(data1 == m | data2 == m)
                                 
                                 # Always put m in front
                                 for(i in 1:nrow(df0_all_m)){
                                   
                                   if(df0_all_m$data1[i] != m){
                                     
                                     df0_all_m$data2[i] <- df0_all_m$data1[i]
                                     df0_all_m$data1[i] <- m
                                     df0_all_m$prob[i] <- 1-df0_all_m$prob[i]
                                   }
                                   
                                 }
                                 
                                 df1_all_m <- data.frame(cluster = factor(paste('cluster', df0_all_m$cluster),
                                                                          levels = paste('cluster', 1:J)),
                                                         mouse.pair = paste('(', df0_all_m$data1, ',', df0_all_m$data2, ')', sep = ''),
                                                         probability = df0_all_m$prob)
                                 
                                 
                                 probability_plot_m <- df1_all_m %>%
                                   ggplot(mapping = aes(x = cluster,
                                                        y = mouse.pair,
                                                        fill = probability))+
                                   geom_tile()+
                                   theme_bw()+
                                   scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                                                        values = scales::rescale(c(0, 0.25, 1)),
                                                        limits = c(0,1))+
                                   xlab('cluster')+
                                   ylab('mouse pairs')+
                                   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
                                 
                                 probability_plot_m
                               })
  
  return(list('all_obs' = df0_all,
              'plot.for.each.data' = plot.for.each.data,
              'significant_obs' = df0_significant,
              'probability_plot' = probability_plot))

}

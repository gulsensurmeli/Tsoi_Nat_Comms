#' Heat-map to show projection strength
#'

pp.standard.ordering <- function(Y,
                                 Z,
                                 regions.name = NULL,
                                 title = ''){


  J <- max(unlist(Z))
  R <- nrow(Y[[1]])

  if(is.null(regions.name)){

    regions.name <- paste('region', 1:R)
  }

  Y.cbind <- do.call(cbind, Y)
  Y.prop.cbind <- matrix(as.vector(Y.cbind)/rep(colSums(Y.cbind), each = R),
                         nrow = R)

  
  # Re-ordered projection strength
  df0 <- lapply(1:J,
                function(j){

                  mx0 <- matrix(Y.prop.cbind[,which(unlist(Z)==j)],
                                nrow = R)
                  
                  strong.proj.j <- which(rowSums(mx0) == max(rowSums(mx0)))[1]
                  
                  mx0[,order(mx0[strong.proj.j,])]
                })
  

  df0 <- do.call(cbind, df0)
  
  N.j <- sapply(1:J,
                function(j) length(which(unlist(Z)==j)))

  # Convert to data frame for plotting
  Y.prop.df <- data.frame(region = rep(regions.name,
                                       ncol(Y.prop.cbind)),

                          neuron = rep(1:ncol(Y.prop.cbind),
                                       each = length(regions.name)),

                          projection.strength = as.vector(df0))


  Y.prop.df %>%
    ggplot(mapping = aes(x = factor(region, levels = regions.name),
                         y = neuron,
                         fill = projection.strength))+
    geom_tile()+
    theme_bw()+
    scale_fill_gradientn(colours = c('white', 'yellow', 'red'),
                         values = scales::rescale(c(0, 0.25, 1)),
                         limits = c(0,1))+
    guides(fill=guide_legend(title="projection strength"))+
    xlab('region')+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          plot.title = element_text(size=12),
          plot.background = element_blank() ,
          panel.grid.major = element_blank() ,
          panel.grid.minor = element_blank() ,
          panel.border = element_blank() ,
          panel.background = element_blank()) +
    #draws x and y axis line  
    theme(axis.line = element_line(color = 'black'))+
    ylab('neurons')+

    geom_hline(yintercept = cumsum(N.j)[-length(cumsum(N.j))],
               color = 'blue',
               linetype = 'dashed',
               linewidth = 0.05)+

    ggtitle(title)


}


pp.standard.ordering.EC <- function(Y,
                                  Z,
                                  regions.name = NULL,
                                  title = '',
                                  EC_label){


  J <- max(unlist(Z))
  R <- nrow(Y[[1]])

  if(is.null(regions.name)){

    regions.name <- paste('region', 1:R)
  }

  Y.cbind <- do.call(cbind, Y)
  Y.prop.cbind <- matrix(as.vector(Y.cbind)/rep(colSums(Y.cbind), each = R),
                         nrow = R)

  EC_label_combine <- unlist(EC_label)
  
  N.j <- sapply(1:J,
                function(j) length(which(unlist(Z)==j)))
  

  # Re-ordered projection strength
  df0 <- lapply(1:J,
                function(j){

                  # Projection strengths of neurons in cluster j
                  mx0 <- matrix(Y.prop.cbind[,which(unlist(Z)==j)],
                                nrow = R)

                  # EC labels
                  EC_label_j <- EC_label_combine[which(unlist(Z)==j)]

                  # Projection strengths of LEC neurons in cluster j
                  mx0_LEC <- matrix(mx0[,which(EC_label_j == 'LEC')],
                                    nrow = R)

                  # Projection strengths of MEC neurons in cluster j
                  mx0_MEC <- matrix(mx0[,which(EC_label_j == 'MEC')],
                                    nrow = R)

                  strong.proj.j <- which(rowSums(mx0) == max(rowSums(mx0)))[1]

                  # Reorder projecting strengths for LEC and MEC
                  mx0_LEC <- matrix(mx0_LEC[,order(mx0_LEC[strong.proj.j,])],
                                    nrow = R)

                  mx0_MEC <- matrix(mx0_MEC[,order(mx0_MEC[strong.proj.j,])],
                                    nrow = R)

                  cbind(mx0_LEC, -mx0_MEC)

                })


  df0 <- do.call(cbind, df0)

  # Convert to data frame for plotting
  Y.prop.df <- data.frame(region = rep(regions.name,
                                       ncol(Y.prop.cbind)),

                          neuron = rep(1:ncol(Y.prop.cbind),
                                       each = length(regions.name)),

                          projection.strength = as.vector(df0))


  Y.prop.df %>%
    ggplot(mapping = aes(x = factor(region, levels = regions.name),
                         y = neuron,
                         fill = projection.strength))+
    geom_tile()+
    theme_bw()+
    scale_fill_gradientn(colours = c('#CB6400','white', '#16697A'),
                         values = scales::rescale(c(-1, 0, 1)),
                         limits = c(-1,1))+
    guides(fill=guide_legend(title="projection strength"))+
    xlab('region')+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          plot.title = element_text(size=12),
          plot.background = element_blank() ,
          panel.grid.major = element_blank() ,
          panel.grid.minor = element_blank() ,
          panel.border = element_blank() ,
          panel.background = element_blank()) +
    #draws x and y axis line
    theme(axis.line = element_line(color = 'black'))+
    ylab('neurons')+

    geom_hline(yintercept = cumsum(N.j)[-length(cumsum(N.j))],
               color = 'blue',
               linetype = 'dashed',
               linewidth = 0.05)+

    ggtitle(title)+
    scale_fill_gradient2(high = '#16697A', mid = "white", low = '#DB6400')


}


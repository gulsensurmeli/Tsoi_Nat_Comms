
#' Plot to show the number of neurons in each cluster, and the contribution of each mouse in each cluster
#'

#-- General function
opt.clustering.frequency <- function(clustering,
                                     main = ''){

  cluster.names <- paste( 1:length(unique(unlist(clustering))))

  loop.result <- lapply(1:length(clustering), function(m){

    data.frame(cluster = factor(cluster.names,
                                levels = cluster.names),
               mouse = paste(m),
               counts = as.vector(table(factor(clustering[[m]],
                                               levels = 1:length(cluster.names)))))
  })

  # Frequency table
  z.frequency <- do.call(rbind, loop.result)

  # Bar chart
  z.frequency %>%
    ggplot(mapping = aes(x = cluster, y = counts, fill = mouse))+
    geom_bar(stat="identity")+
    theme_bw()+
    ylab('Number of neurons')+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(main)
}


# For large clusters
opt.clustering.frequency.large <- function(clustering,
                                           main = '',
                                           large.cluster.index){
  
  cluster.names <- paste(1:length(unique(unlist(clustering))))
  
  loop.result <- lapply(1:length(clustering), function(m){
    
    data.frame(cluster = factor(cluster.names,
                                levels = cluster.names),
               mouse = paste(m),
               counts = as.vector(table(factor(clustering[[m]],
                                               levels = 1:length(cluster.names)))))
  })
  
  # Frequency table
  z.frequency <- do.call(rbind, loop.result) %>%
    filter(cluster %in% paste(large.cluster.index))
  
  # Bar chart
  z.frequency %>%
    ggplot(mapping = aes(x = cluster, y = counts, fill = mouse))+
    geom_bar(stat="identity")+
    theme_bw()+
    ylab('Number of neurons')+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    ggtitle(main)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}



##-- For the EC data
opt.clustering.frequency1 <- function(clustering,
                                      main='',
                                      EC_label){
  
  
  cluster.names <- paste(1:length(unique(unlist(clustering))))
  
  
  # Table of LEC and MEC
  LEC_df <- data.frame(cluster = factor(cluster.names,
                                        levels = cluster.names),
                       EC_label = 'LEC',
                       counts = as.vector(table(factor(unlist(clustering)[which(unlist(EC_label) == 'LEC')],
                                                       levels = 1:length(cluster.names)))))
  
  MEC_df <- data.frame(cluster = factor(cluster.names,
                                        levels = cluster.names),
                       EC_label = 'MEC',
                       counts = as.vector(table(factor(unlist(clustering)[which(unlist(EC_label) == 'MEC')],
                                                       levels = 1:length(cluster.names)))))
  
  df <- rbind(LEC_df, MEC_df)
  
  # Plot to show proportion
  
  ggplot(df, mapping = aes(x = cluster, y = counts, fill = EC_label))+
    geom_bar(stat="identity")+
    theme_bw()+
    ylab('Number of LEC and MEC neurons')+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
}


opt.clustering.frequency2 <- function(clustering,
                                      main='',
                                      EC_label){
  
  
  cluster.names <- paste(1:length(unique(unlist(clustering))))
  
  
  # Table of LEC and MEC
  LEC_df <- data.frame(cluster = factor(cluster.names,
                                        levels = cluster.names),
                       EC_label = 'LEC',
                       counts = as.vector(table(factor(unlist(clustering)[which(unlist(EC_label) == 'LEC')],
                                                       levels = 1:length(cluster.names)))))
  
  MEC_df <- data.frame(cluster = factor(cluster.names,
                                        levels = cluster.names),
                       EC_label = 'MEC',
                       counts = as.vector(table(factor(unlist(clustering)[which(unlist(EC_label) == 'MEC')],
                                                       levels = 1:length(cluster.names)))))
  
  df <- rbind(LEC_df, MEC_df)
  
  df2 <- df %>%
    group_by(cluster) %>%
    summarise(sum_counts = sum(counts))
  
  # Plot to show proportion
  ggplot(df, 
         aes(x = cluster, y = counts, fill = EC_label)) +
    geom_col(position = "fill")+
    theme_bw()+
    ylab('Proportion of LEC and MEC neurons')+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    annotate(geom = 'text', 
             label = df2$sum_counts,
             x = (1:nrow(df2) - 0.5), y = 0.95, hjust = 0, vjust = 1, size = 3,
             angle = 90)
  
}
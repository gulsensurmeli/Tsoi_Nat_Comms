set.seed(100)

# Suppose there are four mice
M <- 4

# And 100 neurons in mouse 1 and 4, 200 neurons in mouse 2 and 3
C <- c(100,200,200,100)

# and 3 brain regions
R <- 3

# There are three clusters
J <- 3

# with projection strengths
q <- matrix(c(0.5,0.4,0.1,
              0.2,0.2,0.6,
              0.1,0.5,0.4),
            
            byrow = TRUE,
            nrow = 3,
            ncol = 3)

N_data <- lapply(1:M,
                 function(m) sample(x = 20:50,
                                    size = C[m],
                                    replace = TRUE))

# Mouse 1: LEC only
# Mouse 2: neuron 1:100 LEC, neuron 101:200 MEC
# Mouse 3: neuron 1:100 LEC, neuron 101:200 MEC
# Mouse 4: MEC only

delta <- matrix(c(0.2,0.5,0.3,
                  0.5,0.4,0.1,
                  0.5,0.2,0.3,
                  0.2,0.3,0.5),
                
                nrow = J,
                ncol = M)


# Gamma
gamma <- rep(10,J)

# Effect of changing from LEC to MEC
beta <- c(-2, 2, 0)

# w_jm(xi = 0)
w_jm_LEC <- matrix(c(0.2,0.5,0.3,
                     0.5,0.4,0.1,
                     0.5,0.2,0.3,
                     0,0,0),
                   
                   nrow = J,
                   ncol = M)

w_jm_LEC_df <- data.frame(w_jm = as.vector(w_jm_LEC),
                          cluster = rep(1:J, M),
                          mouse = rep(1:M, each = J))

w_jm_LEC_df$cluster <- factor(w_jm_LEC_df$cluster,
                              levels = 1:J)

w_jm_LEC_df$mouse <- factor(w_jm_LEC_df$mouse,
                            levels = 1:M)

w_jm_LEC_df %>%
  ggplot()+
  geom_bar(mapping = aes(x = cluster,
                         y = w_jm,
                         fill = mouse),
           stat = 'identity',
           position = 'dodge')+
  theme_bw()+
  ggtitle('LEC allocation proportion')

# w_jm(xi = 1)
w_j1_MEC <- rep(0,3)

w_j2_MEC <- (delta[,2]*exp(beta))/sum(delta[,2]*exp(beta))

w_j3_MEC <- (delta[,3]*exp(beta))/sum(delta[,3]*exp(beta))

w_j4_MEC <- (delta[,4]*exp(beta))/sum(delta[,4]*exp(beta))

w_jm_MEC <- matrix(c(w_j1_MEC,
                     w_j2_MEC,
                     w_j3_MEC,
                     w_j4_MEC),
                   
                   nrow = J,
                   ncol = M)


w_jm_MEC_df <- data.frame(w_jm = as.vector(w_jm_MEC),
                          cluster = rep(1:J, M),
                          mouse = rep(1:M, each = J))

w_jm_MEC_df$cluster <- factor(w_jm_MEC_df$cluster,
                              levels = 1:J)

w_jm_MEC_df$mouse <- factor(w_jm_MEC_df$mouse,
                            levels = 1:M)

w_jm_MEC_df %>%
  ggplot()+
  geom_bar(mapping = aes(x = cluster,
                         y = w_jm,
                         fill = mouse),
           stat = 'identity',
           position = 'dodge')+
  theme_bw()+
  ggtitle('MEC allocation proportion')


# Allocations

# Data 1: Only LEC
Z_1 <- sample(1:J,
              size = 100,
              replace = TRUE,
              prob = w_jm_LEC[,1])


# Data 2: LEC and MEC
Z_2 <- c(sample(1:J,
                size = 100,
                replace = TRUE,
                prob = w_jm_LEC[,2]),
         
         sample(1:J,
                size = 100,
                replace = TRUE,
                prob = w_jm_MEC[,2]))

# Data 3: LEC and MEC
Z_3 <- c(sample(1:J,
                size = 100,
                replace = TRUE,
                prob = w_jm_LEC[,3]),
         
         sample(1:J,
                size = 100,
                replace = TRUE,
                prob = w_jm_MEC[,3]))

# Data 4: MEC only
Z_4 <- sample(1:J,
              size = 100,
              replace = TRUE,
              prob = w_jm_MEC[,4])


Z <- list(Z_1,
          Z_2,
          Z_3,
          Z_4)


# Simulated Y
Y <- lapply(1:M,
            function(m){
              
              Y_m <- lapply(1:C[m],
                            function(c){
                              
                              matrix(rmultinom(n = 1,
                                               size = N_data[[m]][c],
                                               prob = q[Z[[m]][c],]),
                                     ncol = 1)
                            })
              
              do.call(cbind, Y_m)
            })


# Projection strength of every neuron
Y.prop <- NULL
for(m in 1:M){
  
  Y.prop[[m]] <- apply(Y[[m]], 2, function(x) x/sum(x))
}


Y.prop.bind <- do.call(cbind, Y.prop)
Z.bind <- unlist(Z)
Y.prop.by.j <- lapply(1:3, function(j) Y.prop.bind[, which(Z.bind==j)])
Y.prop.by.j <- do.call(cbind, Y.prop.by.j)

Y.prop.df <- data.frame(region = rep(1:J, sum(C)),
                        neuron = rep(1:sum(C), each = J),
                        projection.strength = as.vector(Y.prop.by.j))


cluster.length <- sapply(1:max(Z.bind), function(j) length(which(unlist(Z)==j)))



# Plot of neuron-wise pp by true Z
Y.prop.df %>%
  ggplot(mapping = aes(x = region, y = neuron, fill = projection.strength))+
  geom_tile()+
  theme_bw()+
  scale_fill_gradient(low = 'white',  high = 'red', limits = c(0,1))+
  guides(fill=guide_legend(title="Projection probability"))+
  geom_hline(yintercept = cumsum(cluster.length)[-J],
             color = 'blue',
             linetype = 'dashed')+
  ggtitle('Projecion strength of each neuron')



# x
x <- NULL
x[[1]] <- rep(0, 100)
x[[2]] <- c(rep(0, 100), rep(1, 100))
x[[3]] <- c(rep(0, 100), rep(1, 100))
x[[4]] <- rep(1, 100)

# Run MCMC

mcmc_all_test <- mcmc_run_all(Y = Y,
                              x = x,
                              J = 5,
                              print_Z = TRUE,
                              iter_update = 100,
                              a_gamma = 10,
                              b_gamma = 1,
                              a0 = 1,
                              b0 = 1,
                              a_alpha = 1,
                              b_alpha = 1,
                              s = 5)


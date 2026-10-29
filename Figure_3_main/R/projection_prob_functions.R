## Compute conditional probabilities of projection
library(coda)
library(latex2exp)

# Compute the marginal probability of zero counts for each region for a DirMult mixture
ppzero <- function(n, alpha, w){
  R = dim(alpha)[2]
  J = length(w)
  pp = matrix(0,R,1)
  for (j in 1:J){
    pp = pp+w[j]*apply(matrix(1:R,R,1), 1, function(r){ddirmultinomial(c(0,n),n,c(alpha[j,r],sum(alpha[j,])-alpha[j,r]))})
  }
  return(pp)
}

# Compute the marginal probability of zero counts for each pair of region for a DirMult mixture
ppzero_pair <- function(n, alpha, w){
  R = dim(alpha)[2]
  J = length(w)
  pp = matrix(0,R,R)
  for (j in 1:J){
    for (r1 in 1:R){
      pp[r1,] = pp[r1,]+w[j]*apply(matrix(1:R,R,1), 1, function(r2){ddirmultinomial(c(0,0,n),n,c(alpha[j,r1],alpha[j,r2],sum(alpha[j,])-alpha[j,r1]-alpha[j,r2]))})
    }
  }
  return(pp)
}

# Compute the marginal probability of non-zero counts for each region and the conditional
# proability of non-zero counts for region one given another
projection_prob= function(m,N,mcmc_all_out){
  
  R= mcmc_all_out$R
  I = length(mcmc_all_out$alpha_output)
  
  #output
  pprob_pair_array = array(0,dim= c(I,R,R))
  pprob_mat = matrix(0,I,R)
  
  for(i in 1:I){
     wjm = mcmc_all_out$omega_J_M_output[[i]][,m]
     qstarjm = mcmc_all_out$q_star_1_J_output[[i]]
     gammastarjm = mcmc_all_out$gamma_star_1_J_output[[i]]
     pzero_i = ppzero(N,qstarjm*gammastarjm, wjm)
     pzero_pair_i = ppzero_pair(N,qstarjm*gammastarjm, wjm)
     pprob_mat[i,] = 1- pzero_i
     pprob_pair_temp =  (1-(matrix(pzero_i, R,R) + t(matrix(pzero_i, R,R)) - pzero_pair_i))/matrix(1-pzero_i, R,R)
     diag(pprob_pair_temp) = 1
     pprob_pair_array[i,,] = pprob_pair_temp 
  }
  pprob = apply(pprob_mat,2,mean)
  pprob_pair = apply(pprob_pair_array,c(2,3),mean)
  pprob_ci = apply(pprob_mat,2,function(x){HPDinterval(as.mcmc(x))})
  pprobpair_ci = apply(pprob_pair_array,c(2,3),function(x){HPDinterval(as.mcmc(x))})
  return(list(pprob = pprob, pprob_pair = pprob_pair, pprob_ci = pprob_ci, pprob_pair_ci = pprobpair_ci))
}


# CI for the difference in the proability of non-zero counts for region one given another
# between two mice
projection_prob_d= function(m1,m2,N,mcmc_all_out){
  
  R= mcmc_all_out$R
  I = length(mcmc_all_out$alpha_output)
  
  #output
  pprob_pair_array_1 = array(0,dim= c(I,R,R))
  pprob_mat_1 = matrix(0,I,R)
  
  pprob_pair_array_2 = array(0,dim= c(I,R,R))
  pprob_mat_2 = matrix(0,I,R)
  
  for(i in 1:I){
    wjm1 = mcmc_all_out$omega_J_M_output[[i]][,m1]
    wjm2 = mcmc_all_out$omega_J_M_output[[i]][,m2]
    qstarjm = mcmc_all_out$q_star_1_J_output[[i]]
    gammastarjm = mcmc_all_out$gamma_star_1_J_output[[i]]
    pzero_i1 = ppzero(N,qstarjm*gammastarjm, wjm1)
    pzero_i2 = ppzero(N,qstarjm*gammastarjm, wjm2)
    pzero_pair_i1 = ppzero_pair(N,qstarjm*gammastarjm, wjm1)
    pzero_pair_i2 = ppzero_pair(N,qstarjm*gammastarjm, wjm2)
    pprob_mat_1[i,] = 1- pzero_i1
    pprob_mat_2[i,] = 1- pzero_i2
    pprob_pair_temp1 =  (1-(matrix(pzero_i1, R,R) + t(matrix(pzero_i1, R,R)) - pzero_pair_i1))/matrix(1-pzero_i1, R,R)
    pprob_pair_temp2 =  (1-(matrix(pzero_i2, R,R) + t(matrix(pzero_i2, R,R)) - pzero_pair_i2))/matrix(1-pzero_i2, R,R)
    diag(pprob_pair_temp1) = 1
    diag(pprob_pair_temp2) = 1
    pprob_pair_array_1[i,,] = pprob_pair_temp1
    pprob_pair_array_2[i,,] = pprob_pair_temp2 
  }
  pprob_d = apply(pprob_mat_1-pprob_mat_2,2,mean)
  pprob_pair_d = apply(pprob_pair_array_1-pprob_pair_array_2,c(2,3),mean)
  pprob_d_ci = apply(pprob_mat_1-pprob_mat_2,2,function(x){HPDinterval(as.mcmc(x))})
  pprobpair_d_ci = apply(pprob_pair_array_1-pprob_pair_array_2,c(2,3),function(x){HPDinterval(as.mcmc(x))})
  return(list(pprob = pprob_d, pprob_pair = pprob_pair_d, pprob_ci = pprob_d_ci, pprob_pair_ci = pprobpair_d_ci))
}

# Compute expected number of counts/N for each region for a DirMult mixture
epstrength <- function(alpha, w){
  R = dim(alpha)[2]
  J = length(w)
  eps = matrix(0,R,1)
  for (j in 1:J){
    eps = eps+w[j]*alpha[j,]/sum(alpha[j,])
  }
  return(eps)
}

# Sample from DirMult mixture
rpstrength <- function(n, alpha, w){
  R = dim(alpha)[2]
  J = length(w)
  j = sample(c(1:J),1, prob = w)
  ps = rdirmnom(1,n,alpha[j,])/n
  return(ps)
}


# Compute the variance of counts/N for each region for a DirMult mixture
vpstrength <- function(n,alpha, w){
  R = dim(alpha)[2]
  J = length(w)
  eps = matrix(0,R,1)
  vps = matrix(0,R,1)
  for (j in 1:J){
    gamma0 = sum(alpha[j,])
    eps = eps+w[j]*alpha[j,]/gamma0
    vps = vps+w[j]*(1/n*(alpha[j,]/gamma0)*(1-alpha[j,]/gamma0)*(n+gamma0)/(1+gamma0)+(alpha[j,]/gamma0)^2)
  }
  vps = vps - eps^2
  return(vps)
}


# Compute the covariance of y1/N and y2?N counts for each pair of region for a DirMult mixture
cpstrength <- function(n, alpha, w){
  R = dim(alpha)[2]
  J = length(w)
  cps = matrix(0,R,R)
  eps = matrix(0,R,1)
  for (j in 1:J){
    gamma0 = sum(alpha[j,])
    eps = eps+w[j]*alpha[j,]/gamma0
    for (r1 in 1:R){
      cps[r1,] =  cps[r1,] + w[j]*(-1/n*(alpha[j,r1]/gamma0)*(alpha[j,]/gamma0)*(n+gamma0)/(1+gamma0)+(alpha[j,r1]/gamma0)*(alpha[j,]/gamma0))
      }
  }
  cps = cps - eps%*%t(eps)
  diag(cps) = vpstrength(n,alpha, w)
  return(cps)
}

# Compute the posterior expected number of counts/N for each region 
post_epstrength = function(m,N,mcmc_all_out){
  
  R= mcmc_all_out$R
  I = length(mcmc_all_out$alpha_output)
  
  #output
  eps_mat = matrix(0,I,R)
  rps_mat = matrix(0,I,R)
  
  for(i in 1:I){
    wjm = mcmc_all_out$omega_J_M_output[[i]][,m]
    qstarjm = mcmc_all_out$q_star_1_J_output[[i]]
    gammastarjm = mcmc_all_out$gamma_star_1_J_output[[i]]
    eps_mat[i,] = epstrength(qstarjm*gammastarjm, wjm)
    rps_mat[i,] = rpstrength(N,qstarjm*gammastarjm, wjm)
  }
  eps = apply(eps_mat,2,mean)
  eps_ci = apply(eps_mat,2,function(x){HPDinterval(as.mcmc(x))})
  rps_ci = apply(rps_mat,2,function(x){HPDinterval(as.mcmc(x))})
  return(list(eps = eps, eps_ci = eps_ci, rps_ci = rps_ci))
}

# Compute the posterior covariance of counts/N for each pair of regions
post_cpstrength= function(m,N,mcmc_all_out){
  
  R= mcmc_all_out$R
  I = length(mcmc_all_out$alpha_output)
  
  #output
  cps_array = array(0,dim= c(I,R,R))
  eps_mat = matrix(0,I,R)
  eps2_array = array(0,dim= c(I,R,R))
  
  for(i in 1:I){
    wjm = mcmc_all_out$omega_J_M_output[[i]][,m]
    qstarjm = mcmc_all_out$q_star_1_J_output[[i]]
    gammastarjm = mcmc_all_out$gamma_star_1_J_output[[i]]
    eps_mat[i,] = epstrength(qstarjm*gammastarjm, wjm)
    eps2_array[i,,] = eps_mat[i,]%*%t(eps_mat[i,]) 
    cps_array[i,,] = cpstrength(N,qstarjm*gammastarjm, wjm) 
  }
  eps = apply(eps_mat,2,mean)
  eps2 = apply(eps2_array,c(2,3),mean)
  cps = apply(cps_array,c(2,3),mean) + eps2 - eps%*%t(eps)
  cps_ci = apply(cps_array,c(2,3),function(x){HPDinterval(as.mcmc(x))})
  return(list(cps = cps, cps_ci = cps_ci))
}




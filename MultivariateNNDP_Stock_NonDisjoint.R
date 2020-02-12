##Stock Multivariate NNDP Code
##Pass cross validated answers into this code

#Load Packages

library(FastKNN)
library(FNN)
library(benchden)
library(dirichletprocess)
library(ks)
library(mvtnorm)
library(Rcpp)
library(plyr)
library(dplyr)

#Stock Multivariate NNDP Code

mNNDP<-function(x, MC, k, inputpt, mu0, nu0, gamma0, psi0, alpha)
{
  ## Declarations
  n = dim(x)[1]
  p = dim(x)[2]
  s = dim(inputpt)[1]
  f_stor = matrix(0, nrow = s, ncol = MC)
  
  ##KNN
  
  time_KNN_1 = Sys.time()
  
  KNN_full = get.knn(x, k = (k-1))
  knn_list = KNN_full$nn.index
  knn_metric = KNN_full$nn.dist
  
  all_nbd = cbind(1:n, knn_list)
  
  nbd_mean = matrix(0, nrow = n, ncol = p)
  nbd_var = array(0, dim = c(p, p, n))
  nbd_range = rep(0, n)
  
  i=1
  
  for(i in 1:n)
  {
    nbd_mean[i,] = colMeans(x[all_nbd[i,],])
    nbd_var[,,i] = var(x[all_nbd[i,],])
    nbd_range[i] = knn_metric[i,(k-1)]
  }
  
  time_KNN_2 = Sys.time()
  
  time_KNN = time_KNN_2 - time_KNN_1
  
  ####Now finding posterior hyperparameters
  
  ###nu_n
  
  nu_n = nu0 + k
  gamma_n = gamma0 + k 
  
  ###mu_i
  
  priorwt = nu0/nu_n
  datawt = k/nu_n
  
  mu_vec = (priorwt*t(replicate(n,mu0))) + (datawt*nbd_mean)
  
  ###psi_n
  
  part1 = replicate(n, psi0)
  part2 = (k-1)*nbd_var
  part3 = array(0, dim = c(p,p,n))
  
  i=1
  
  for(i in 1:n)
  {
    part3[,,i] = ((k*nu0)/(nu_n))*(nbd_mean[i,]-mu0)%*%t(nbd_mean[i,]-mu0)
  }
  
  psi = (part1 + part2 + part3)
  
  ###Disjointed neighbourhoods
  
  #Disjointification
  
  # disjoint_nbds = rep(0,n)
  # length_vec = NULL
  # lost_nbd_index = NULL
  # 
  # i=1
  # 
  # disjoint_nbds[1:k] = all_nbd[1,]
  # length_vec = k
  # 
  # for(i in 2:n)
  # {
  #   temp = setdiff(all_nbd[i,],disjoint_nbds)
  #   
  #   if(length(temp)>=1)
  #   {
  #     disjoint_nbds[(sum(length_vec)+1):(sum(length_vec)+length(temp))] = temp
  #     
  #     length_vec = c(length_vec, length(temp))
  #     
  #   }else
  #   {
  #     lost_nbd_index = c(lost_nbd_index, i)
  #   }
  # }
  # 
  # lenvec = rep(0,n)
  # lenvec[-lost_nbd_index] = length_vec
  
  wtvec = rep((alpha+1),n)
  
  #estimate = m.kerfn(eval = inputpt, loc = mu_vec, bw = psi*((nu0+k+1)/((gamma0+k+1-p)*(nu0+k)))
  #, df = (gamma0+k+1-p), wt = wtvec/sum(wtvec))
  
  psi_inv = array(apply(psi,3,solve),c(p,p,n))
  psi_inv_list = alply(psi_inv, 3)
  
  ####Now start MC
  
  i=1
  j=1
  
  time_MC_1 = Sys.time()
  
  #f_stor = fdraws(n, p, s, MC, wtvec, mu_vec, psi_inv, gamma_n, nu_n, inputpt)
  
  f_stor = fdraws(n, p, s, MC, wtvec, mu_vec, psi_inv_list, gamma_n, nu_n, inputpt)
  
  time_MC_2 = Sys.time()
  time_MC = time_MC_2 - time_MC_1
  time_vec = c(time_KNN, time_MC)
  
  nbd_stat = list("nbd_mean" = nbd_mean, "nbd_var" = nbd_var, "nbd_range" = nbd_range, "psi" = psi, "nbd_wt" = wtvec/(n*(alpha+1)))
  
  output = list("f_stor" = f_stor, "k" = k, "n" = n, "MC" = MC, "inputpt" = inputpt, 
                "data" = x, "hyp_par" = list(mu0, nu0, psi0, alpha), 
                "time" = time_vec, "nbd_stat" = nbd_stat)
  
  return(output)
  
}

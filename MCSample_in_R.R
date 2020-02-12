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
  
  disjoint_nbds = rep(0,n)
  length_vec = NULL
  lost_nbd_index = NULL
  
  i=1
  
  disjoint_nbds[1:k] = all_nbd[1,]
  length_vec = k
  
  for(i in 2:n)
  {
    temp = setdiff(all_nbd[i,],disjoint_nbds)
    
    if(length(temp)>=1)
    {
      disjoint_nbds[(sum(length_vec)+1):(sum(length_vec)+length(temp))] = temp
      
      length_vec = c(length_vec, length(temp))
      
    }else
    {
      lost_nbd_index = c(lost_nbd_index, i)
    }
  }
  
  lenvec = rep(0,n)
  lenvec[-lost_nbd_index] = length_vec
  
  wtvec = (alpha + lenvec)
  
  #estimate = m.kerfn(eval = inputpt, loc = mu_vec, bw = psi*((nu0+k+1)/((gamma0+k+1-p)*(nu0+k)))
  #, df = (gamma0+k+1-p), wt = wtvec/sum(wtvec))
  
  ####Now start MC
  
  i=1
  j=1
  
  time_MC_1 = Sys.time()
  
  #For first method
  
  #psi_inv = alply(array(apply(psi,3,solve),c(p,p,n)),3)
  psi_inv = array(apply(psi,3,solve),c(p,p,n))
  
  for(j in 1:MC)
  {
    
    ##Dirichlet vector generation
    
    gam = rgamma(n = n, shape = wtvec, rate=1)
    pi_vec = gam/sum(gam)
    
    ##Neighbourhood mean and variance
    
    sig = array(apply(array(apply(psi_inv, 3, rWishart,
          n = 1, df = gamma_n),c(p,p,n)), 3, solve),c(p,p,n))
    
    theta = t(apply(sig/nu_n, 3, rmvnorm, n = 1, mean = rep(0,p)))
    theta = theta + mu_vec
    
    sig = alply(sig, 3)
    theta = alply(as.array(theta), 1)
    
    ####Leader combines
    
    evalmat = mapply(dmvnorm, theta, sig, MoreArgs = list("x" = inputpt))
    
    f_stor[,j] = evalmat%*%pi_vec
    
    print(j)
    
  }
  
  time_MC_2 = Sys.time()
  time_MC = time_MC_2 - time_MC_1
  time_vec = c(time_KNN, time_MC)
  
  nbd_stat = list("nbd_mean" = nbd_mean, "nbd_var" = nbd_var, "nbd_range" = nbd_range, 
                  "psi" = psi, "nbd_wt" = wtvec/(n*(alpha+1)))
  
  output = list("f_stor" = f_stor, "estimate"=estimate, "k" = k, "n" = n, "MC" = MC, 
                "inputpt" = inputpt, "data" = x, "hyp_par" = list(mu0, nu0, psi0, alpha), 
                "time" = time_vec, "nbd_stat" = nbd_stat)
  
  return(output)
  
}

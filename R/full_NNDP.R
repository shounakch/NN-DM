#### NN-DM Pseudo-Posterior Mean ####

NNDM_postMean<-function(x, k, inputpt, mu0, nu0, gamma0, psi0)
{
  
  ## Load packages ##
  
  library(FNN)
  library(plyr)
  library(dplyr)
  
  ## Turn everything into matrix ##
  
  x = as.matrix(x)
  inputpt = as.matrix(inputpt)
  
  ## Declarations
  n = dim(x)[1]
  p = dim(x)[2]
  s = dim(inputpt)[1]
  f_hat = rep(0, s)
  
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
    
    if(p == 1)
    {
      
      nbd_mean[i,] = as.vector(mean(x[all_nbd[i,],]))
      nbd_var[,,i] = as.matrix(var(x[all_nbd[i,],]))
      nbd_range[i] = knn_metric[i,(k-1)]
      
    }else {
      
      nbd_mean[i,] = colMeans(x[all_nbd[i,],])
      nbd_var[,,i] = var(x[all_nbd[i,],])
      nbd_range[i] = knn_metric[i,(k-1)]
      
    }
    
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
  
  if(p == 1) 
  {
    
    mu_vec = (priorwt*rep(as.numeric(mu0), n)) + (datawt*nbd_mean)
    
  } else {
    
    mu_vec = (priorwt*t(replicate(n,mu0))) + (datawt*nbd_mean)
    
  }
  
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
  
  psi_inv = array(apply(psi,3,solve),c(p,p,n))
  psi_inv_list = alply(psi_inv, 3)
  
  #### Now define the bandwidth matrices ####
  
  bw_array = ((nu_n + 1) / (nu_n * (gamma_n - p + 1))) * psi
  bw_inv_array = (nu_n * (gamma_n - p + 1) / (nu_n + 1)) * psi_inv
  bw_inv_list = lapply(psi_inv_list, "*", (nu_n * (gamma_n - p + 1) / (nu_n + 1)))
  bw_inv_list = lapply(bw_inv_list, as.matrix)
  
  df_n = gamma_n - p + 1
  
  # middle_matrix = matrix(0, nrow = s, ncol = n)
  # 
  # for(j in 1:n)
  # {
  #   
  #   K_j = bw_inv_array[,,j]
  #   det_Kj = ifelse(p == 1, K_j, det(K_j))
  #   mu_j = mu_vec[j,]
  #   
  #   for(i in 1:s)
  #   {
  #     
  #     qf_ij = as.numeric(t(inputpt[i,] - mu_j) %*% K_j %*% (inputpt[i,] - mu_j))
  #     
  #     dens_ij = (det_Kj^(0.5)) * ((1 + (qf_ij / df_n))^(-0.5*(df_n + p)))
  #     
  #     middle_matrix[i,j] = dens_ij
  #     
  #   }
  #   
  # }
  
  middle_matrix = NNDM_pmean(n, p, s, mu_vec, bw_inv_list, 
                             df_n, nu_n, inputpt)
  
  norm_const = gamma(0.5*(df_n + p)) / (((df_n * pi)^(0.5*p)) * gamma(0.5*df_n))
  
  f_hat = norm_const * rowMeans(middle_matrix)
  
  output = list("NNDM_PPM" = f_hat, "k" = k, "n" = n, "inputpt" = inputpt, 
                "data" = x, 
                "hyp_par" = list("mu0"=mu0, 
                                 "nu0"=nu0, 
                                 "gamma0"=gamma0, 
                                 "psi0"=psi0))
  
  return(output)
  
}

#### Univariate NN-DM ####

uNNDM<-function(x, 
                MC = 1000, 
                k = 10, 
                inputpt, 
                mu0 = 0,
                nu0 = 0.001, 
                gamma0 = 1, 
                delta0sq = "CV", 
                alpha = "RR", 
                samples = TRUE)
{
  
  #### Load packages ####
  
  library(FNN)
  library(plyr)
  library(dplyr)
  
  #### Define LOOCV and sampler function ####
  
  univ_LOOCV<-function(x, k, mu0, nu0, phi0seq, gamma0)
  {
    
    n = length(x)
    
    i=1
    L = length(phi0seq)
    crit = rep(0, L)
    
    #Grand Neighbourhood
    
    grand_nbd = cbind(1:n, get.knn(x, k = k)$nn.index)
    grand_means = as.vector(sapply(1:n,function(i)(mean(x[grand_nbd[i,]]))))
    grand_vars = as.vector(sapply(1:n,function(i)(var(x[grand_nbd[i,]])*k)))
    
    nun = k + nu0
    gamman = k + gamma0
    pwt = nu0/nun
    dwt = k/nun
    
    A = matrix(0, nrow = n, ncol = n)
    
    A = adj_mat(grand_nbd)
    
    m1 = Tmats_univ(n, L, k, grand_nbd, x, A, phi0seq, grand_means, grand_vars, mu0, nu0, gamma0)
    
    crit=rep(0,L)
    for(l in 1:L)
    {
      crit[l] = mean(log(rowMeans(m1[,,l])))
    }
    
    return(crit)
    
  }
  
  uNNDP<-function(x, MC, k, inputpt, mu0, nu0, gamma0, phi0sq, alpha)
  {
    ## Declarations
    n = length(x)
    s = length(inputpt)
    f_stor = matrix(0, nrow = s, ncol = MC)
    
    ##KNN
    
    time_KNN_1 = Sys.time()
    
    KNN_full = get.knn(x, k = (k-1))
    knn_list = KNN_full$nn.index
    knn_metric = KNN_full$nn.dist
    
    all_nbd = cbind(1:n, knn_list)
    
    nbd_mean = rep(0,n)
    nbd_var = rep(0, n)
    nbd_range = rep(0, n)
    
    for(i in 1:n)
    {
      nbd_mean[i] = mean(x[all_nbd[i,]])
      nbd_var[i] = var(x[all_nbd[i,]])
    }
    
    time_KNN_2 = Sys.time()
    
    time_KNN = time_KNN_2 - time_KNN_1
    
    ####Now finding posterior hyperparameters
    
    ###nu_n
    
    nun = nu0 + k
    
    ###gamma_n
    
    gamman = gamma0 + k
    
    ###mu_i
    
    priorwt = nu0/nun
    datawt = k/nun
    
    mu_vec = (priorwt*mu0) + (datawt*nbd_mean)
    
    ###phi_i_sq
    
    part1 = (gamma0*phi0sq)
    part2 = (k-1)*nbd_var
    part3 = ((k*nu0)/(nun))*(mu0-nbd_mean)^2
    
    phivec = (part1+part2+part3)/gamman
    
    rate_par = (phivec*gamman)/2
    
    ### Dirichlet Weights
    
    wtvec = rep((alpha + 1),n)
    
    ###Start MC
    
    time_MC_1 = Sys.time()
    
    f_stor = fdraws_univ(n, s, MC, wtvec, mu_vec, rate_par, gamman, nun, inputpt)
    
    time_MC_2 = Sys.time()
    time_MC = time_MC_2 - time_MC_1
    time_vec = c(time_KNN, time_MC)
    
    nbd_stat = list("nbd_mean" = nbd_mean, "nbd_var" = nbd_var, "nbd_range" = nbd_range, 
                    "phivec" = phivec)
    
    output = list("f_stor" = f_stor, "k" = k, "n" = n, "MC" = MC, "inputpt" = inputpt, 
                  "data" = x, 
                  "hyp_par" = list("mu0"=mu0, 
                                   "nu0"=nu0, 
                                   "gamma0"=gamma0, 
                                   "delta0sq"=phi0sq, 
                                   "alpha"=alpha), 
                  "time" = time_vec, "nbd_stat" = nbd_stat)
    
    return(output)
    
  }
  
  #### Turn everything into matrix ####
  
  x = as.numeric(x)
  inputpt = as.numeric(inputpt)
  n = dim(x)[1]
  
  ### First choose delta0sq ###
  
  if(mode(delta0sq) == "character")
  {
    
    if(delta0sq == "CV")
    {
      
      print("Performing Leave One Out Cross Validation")
      
      LOOCVfun<-function(delta0sq)
      {
        
        return(-univ_LOOCV(x = x,
                           k = k,
                           mu0 = mu0,
                           nu0 = nu0,
                           phi0seq = delta0sq,
                           gamma0 = gamma0))
        
      }
      
      optim_mod = optimize(f = LOOCVfun,
                           lower = 0.001,
                           upper = 100,
                           maximum = FALSE)
      
      delta0sq_opt = optim_mod$minimum
      
      delta0sq = delta0sq_opt
      
    } 
    
  }
  
  ### Next choose alpha ###
  
  if(mode(alpha) == "character")
  {
    
    psi0 = gamma0 * delta0sq
    
    alpha_num = psi0
    H = (alpha_num * (nu0 + k + 1)/((nu0 + k) * (gamma0 + k)))
    sig_f0 = var(x)
    alpha = H * (sig_f0)^(-1) * ((nu0 + k)^(-1))
    
  }
  
  ########         Next run NN-DM         #######
  
  ### Run only pseudo-posterior mean if samples == FALSE
  
  if(samples == FALSE)
  {
    
    return(NNDM_postMean(x = x, 
                          k = k, 
                          inputpt = inputpt, 
                          mu0 = mu0, 
                          nu0 = nu0, 
                          gamma0 = gamma0, 
                          psi0 = delta0sq * gamma0))
    
  }else
  {
    
    return(uNNDP(x = x,
                 MC = MC,
                 k = k,
                 inputpt = inputpt,
                 mu0 = mu0,
                 nu0 = nu0,
                 gamma0 = gamma0,
                 phi0sq = delta0sq,
                 alpha = alpha))
    
  }
  
}
                                  
#### Multivariate NN-DM ####                                  

mNNDM<-function(x, 
                MC = 1000, 
                k = 10, 
                inputpt, 
                mu0,
                nu0, 
                gamma0, 
                psi0 = "CV", 
                alpha = "RR", 
                samples = TRUE)
{
  
  #### Load packages ####
  
  library(FNN)
  library(plyr)
  library(dplyr)
  
  #### Define cross validation function ####
  
  mult.LOOCV<-function(x, k, mu0, nu0, phi0seq, gamma0, cutoff)
  {
    
    #resample cutoff points out of n for LOOCV
    
    N = dim(as.matrix(x))[1]
    
    if(N<=cutoff){
      
      x = as.matrix(x)
      
      n = dim(x)[1]
      p = dim(x)[2]
      
      #L = length(gamma0seq)
      
      L = length(phi0seq)
      crit = rep(0, L)
      
      #Grand Neighbourhood
      
      grand_nbd = cbind(1:n, get.knn(x, k = k)$nn.index) #Take an extra neighbour
      grand_means = t(sapply(1:n,function(i)(colMeans(as.matrix(x[grand_nbd[i,],])))))
      grand_vars = t(sapply(1:n,function(i)(var(as.matrix(x[grand_nbd[i,],]))*(k/(k+1)))))
      
      gvars = array(0, c(p,p,n))
      i=1
      
      for(i in 1:n)
      {
        gvars[,,i] = (k+1)*matrix(grand_vars[i,],p,p)
      }
      
      ##Each row of grand_vars can be made into a p*p variance matrix
      ##Get adjacency matrix
      
      A = matrix(0, nrow = n, ncol = n)
      
      A=adj_mat(grand_nbd)
      
      t1=Sys.time()
      m1 = Tmats(n, p, L, k, grand_nbd, x, A, phi0seq, grand_means, gvars, mu0, nu0, gamma0)
      t2=Sys.time()
      t2-t1
      
      # t1=Sys.time()
      # m1 = Tmats_gamma0(n, p, L, k, grand_nbd, x, A, gamma0seq, grand_means, gvars, mu0, nu0, phi0sq)
      # t2=Sys.time()
      # t2-t1
      
      crit=rep(0,L)
      for(l in 1:L)
      {
        crit[l] = mean(log(rowMeans(m1[,,l])))
      }
      
      return(crit)
    }else
    {
      ind1 = sample(N, cutoff, replace = FALSE)
      dat = x[ind1,]
      knew = ceiling((k/N)*cutoff)
      return(mult.LOOCV(dat, knew, mu0, nu0, phi0seq, gamma0))
    }
    
  }
  
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
    
    ##Dirichlet Weights
    
    wtvec = rep((alpha+1),n)
    
    psi_inv = array(apply(psi,3,solve),c(p,p,n))
    psi_inv_list = alply(psi_inv, 3)
    
    ####Now start MC sampling
    
    time_MC_1 = Sys.time()
    
    f_stor = fdraws(n, p, s, MC, wtvec, mu_vec, psi_inv_list, gamma_n, nu_n, inputpt)
    
    time_MC_2 = Sys.time()
    time_MC = time_MC_2 - time_MC_1
    
    time_vec = c(time_KNN, time_MC)
    
    nbd_stat = list("nbd_mean" = nbd_mean, "nbd_var" = nbd_var, "nbd_range" = nbd_range, "psi" = psi, "nbd_wt" = wtvec/(n*(alpha+1)))
    
    output = list("f_stor" = f_stor, "k" = k, "n" = n, "MC" = MC, "inputpt" = inputpt, 
                  "data" = x, 
                  "hyp_par" = list("mu0"=mu0, 
                                   "nu0"=nu0, 
                                   "psi0"=psi0, 
                                   "alpha"=alpha), 
                  "time" = time_vec, "nbd_stat" = nbd_stat)
    
    return(output)
    
  }
  
  
  #### Turn everything into matrix ####
  
  x = as.matrix(x)
  inputpt = as.matrix(inputpt)
  n = dim(x)[1]
  p = dim(x)[2]
  
  ### First choose psi0 ###
  
  if(mode(psi0) == "character")
  {
    
    if(psi0 == "CV")
    {
      
      print("Performing Leave One Out Cross Validation")
      
      LOOCVfun<-function(delta0sq)
      {
        
        return(-mult.LOOCV(x = x,
                           k = k,
                           mu0 = mu0,
                           nu0 = nu0,
                           phi0seq = delta0sq,
                           gamma0 = gamma0,
                           cutoff = 1000))
        
      }
      
      optim_mod = optimize(f = LOOCVfun,
                           lower = 0.001,
                           upper = 100,
                           maximum = FALSE)
      
      delta0sq_opt = optim_mod$minimum
      
      psi0 = ((gamma0 - p + 1)*delta0sq_opt) * diag(p)
      
    } 
    
  }
  
  ### Next choose alpha ###
  
  if(mode(alpha) == "character")
  {
    
    alpha_num = mean(diag(psi0))
    H = (alpha_num * (nu0 + k + 1)/((nu0 + k) * (gamma0 + k - p + 1))) * diag(p)
    sig_f0 = var(x)
    alpha = det(H) * ((det(sig_f0))^(-1)) * ((nu0 + k)^(-1))
    
  }
  
  ########         Next run NN-DM         #######
  
  ### Run only pseudo-posterior mean if samples == FALSE
  
  if(samples == FALSE)
  {
    
    return(NNDM_postMean(x = x, 
                          k = k, 
                          inputpt = inputpt, 
                          mu0 = mu0, 
                          nu0 = nu0, 
                          gamma0 = gamma0, 
                          psi0 = psi0))
    
  }else
  {
    
    return(mNNDP(x = x,
                 MC = MC,
                 k = k,
                 inputpt = inputpt,
                 mu0 = mu0,
                 nu0 = nu0,
                 gamma0 = gamma0,
                 psi0 = psi0,
                 alpha = alpha))
    
  }
  
}

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
    nbd_range[i] = knn_metric[i,(k-1)]
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
  # 
  # wtvec = (alpha + lenvec)
  
  wtvec = rep((alpha + 1),n)
  
  time_MC_1 = Sys.time()
  
  ####Now start MC
  
  # i=1
  # j=1
  # 
  # 
  # for(j in 1:MC)
  # {
  # 
  #   ##Neighbourhood mean and variance
  # 
  #   sig_sq = 1/rgamma(n = n, shape = gamman/2, rate = (gamman*phivec)/2)
  #   theta = rnorm(n = n, mean = mu_vec, sd=sqrt(sig_sq/nun))
  # 
  #   ##Dirichlet vector generation
  # 
  #   gam = rgamma(n = n, shape = wtvec, rate=1)
  #   pi_vec = gam/sum(gam)
  # 
  #   ##Leader combines
  # 
  #   f_stor[,j] = mapply(FUN=dnorm,theta,sqrt(sig_sq),MoreArgs=list("x"=inputpt))%*%pi_vec
  # 
  #   if((j%%500)==0)
  #   {
  #     print(j)
  #   }
  # }
  
  f_stor = fdraws(n, s, MC, wtvec, mu_vec, rate_par, gamman, nun, inputpt)
  
  time_MC_2 = Sys.time()
  time_MC = time_MC_2 - time_MC_1
  time_vec = c(time_KNN, time_MC)
  
  nbd_stat = list("nbd_mean" = nbd_mean, "nbd_var" = nbd_var, "nbd_range" = nbd_range, 
                  "phivec" = phivec)
  
  output = list("f_stor" = f_stor, "k" = k, "n" = n, "MC" = MC, "inputpt" = inputpt, 
                "data" = x, "hyp_par" = c(mu0, nu0, gamma0, phi0sq, alpha), 
                "time" = time_vec, "nbd_stat" = nbd_stat)
  
  return(output)
  
}
### NNDP simulations ###

library(benchden)
library(dirichletprocess)
library(ks)
library(FastKNN)
library(FNN)
library(MCMCpack)
library(abind)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(EMMIXcskew)
library(EMMIXskew)
library(doParallel)

#Move all the source files mentioned below to same folder

##Base 1,5,9,17

dnums = c(4,8,12,20)
setwd("~/NNDP/NNDP Multivariate Codes/Updated Codes/New_Simulations_NonDisjoint_NNDP") #Set correct directory where source files are situated
Rcpp::sourceCpp("Faster_LOOCV_Rcpp.cpp")
Rcpp::sourceCpp("MCSampleRcpp.cpp")
source("Mult_LOOCV_Rcpp.R")
source("MultivariateNNDP_Stock_NonDisjoint.R")
source("NNDP_multivariate_simulations.R")

sample_sizes = 200
n_replicates = 20
n_cores = 6

NNDP_final_run = function(n, n_replicates, true_den_num, n_cores)
{
  x_samp = rnd_sim_densities(1, 1, true_den_num)$X
  p = length(x_samp)
  dir.name = paste("n=", n, "_p=", p, "_dnum=", true_den_num, sep = "") ## Add dimension
  dir.create(dir.name)
  
  #setwd(dir.name)
  
  cl = makeForkCluster(n_cores)
  registerDoParallel(cl)
  NNDP_single_run<-function(i, n, true_den_num)
  {
    set.seed((2001+i))
    print(i)
    
    n=n
    k=10 ## Fix k for multivariate k = 10, n = 1000
    mu0=rep(0,p)
    nu0=0.001
    gamma0=p
    alpha = 0.1
    MC = 1000
    
    whole = rnd_sim_densities(n, 500, true_den_num) ## Change to rnd_sim_densities()
    
    x = whole$X
    
    #inputpt = seq(-3,3,length.out = 200) # change appro..
    inputpt = whole$X_new
    trueval = whole$f0
    f0 = trueval
    
    LOOCVfun<-function(phi0sq)
    {
      return(-mult.LOOCV(x, k, mu0, nu0, phi0sq, gamma0))
    }
    
    t.optim1 = proc.time()
    L1 = optim(1, LOOCVfun, lower=0.001, method = "L-BFGS-B")
    t.optim2 = proc.time()
    t.optim = (t.optim2 - t.optim1)[3]
    
    phi0sq = L1$par
    
    psi0=diag(rep((gamma0-p+1)*phi0sq,p))
    
    ##Now do NNDP
    
    t.NNDP1 = proc.time()
    N1 = mNNDP(x, MC, k, inputpt, mu0, nu0, gamma0, psi0, alpha)
    t.NNDP2 = proc.time()
    
    t.NNDP = (t.NNDP2 - t.NNDP1)[3]
    
    f_stor = N1$f_stor
    
    #kde_res = kde(x, H = Hlscv(x), eval.points = inputpt)$estimate
    
    # t.DP1 = proc.time()
    # dpobj = DirichletProcessMvnormal(x)
    # dpDraw = Fit(dpobj, 2000, updatePrior = FALSE, progressBar = TRUE)
    # dpEst = PosteriorFrame(dpDraw, inputpt, ci_size = 0.05, ndraws = 1000)
    # t.DP2 = proc.time()
    # t.DP = (t.DP2 - t.DP1)[3]
    
    #credible_bands1 = apply(res_nndp, 1, quantile, c(0.025, 0.975))
    #post_mean = apply(res_nndp, 1, mean)
    
    output1 = cbind(inputpt, f0, f_stor)
    #output2 = c(t.NNDP + t.optim, t.DP)
    
    colnames(output1) = c(rep("inputpt", p), "f_0(x)", rep("f_hat", MC))
    
    fname1 = paste(dir.name, "/", "run", i, ".csv", sep = "")  ## Add dimension
    #fname2 = paste(dir.name, "/", "run", i, "_time", ".txt", sep = "")
    
    write.csv(output1, fname1)
    #write.csv(output2, fname2)
  }
  parSapply(cl, 1:n_replicates, NNDP_single_run, n = n, true_den_num = true_den_num)
  stopCluster(cl)
}

for(i in 1:length(dnums))
{
  
  t.full1=proc.time()
  
  print(i)
  print(dnums)
  
  NNDP_final_run(sample_sizes[1], n_replicates, dnums[i], n_cores)
  setwd("~/NNDP/NNDP Multivariate Codes/Updated Codes/New_Simulations_NonDisjoint_NNDP")
  
  t.full2=proc.time()
  t.full=(t.full2-t.full1)[3]
  
  print(t.full)
}

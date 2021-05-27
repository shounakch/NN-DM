### NNDP simulations ###

library(benchden)
library(RcppGSL)
library(dirichletprocess)
library(ks)
library(FastKNN)
library(FNN)
library(tensorA)
library(MCMCpack)
library(abind)
library(mvtnorm)
library(CholWishart)
library(expm)
library(matrixStats)
library(Rcpp)
library(doParallel)
library(dirichletprocess)

dnums = c(23,27)
#, 4, 5, 6, 10, 11, 12, 20, 22, 23, 27)
setwd("/home/grad/sc602/NNDP/NNDP Univariate Codes/Updated Codes/Final_Simulations") #Set correct directory where source files are situated
Rcpp::sourceCpp("LOOCV_Univ_Faster.cpp")
Rcpp::sourceCpp("NNDP_MCSample.cpp")
source("LOOCV_Univ_phi0.R")
source("NNDP_Univariate_NonDisjoint.R")

sample_sizes = c(1000)
n_replicates = 20
n_cores = 6

##Choice of n=500, k = 15, alpha = 0.4

NNDP_final_run = function(n, n_replicates, true_den_num, n_cores)
{
  dir.name = paste("New_NonDisjoint_n=", n, "_dnum=", true_den_num, sep = "") ## Add dimension
  dir.create(dir.name)
  cl = makeForkCluster(n_cores)
  registerDoParallel(cl)
  #setwd(dir.name)
  NNDP_single_run = function(i, n, true_den_num)
  {
    
    print(i)
    set.seed((2001+i))
    
    n = n
    k = ceiling(n^(1/3)) ## Fix k for multivariate k = 10, n = 1000
    mu0 = 0
    nu0 = 0.001
    gamma0 = 1
    #alpha = 1/(100+log(n+1))
    alpha = 0.005
    MC = 1000
   
    x = rberdev(n, dnum=true_den_num) ## Change to rnd_sim_densities()
    inputpt = rberdev(500, dnum = true_den_num)
    #inputpt = c(-3,-2,-1,0,1,2,3)
    trueval = dberdev(inputpt, dnum = true_den_num)
    f0 = trueval
    
    LOOCVfun<-function(phi0sq)
    {
      return(-univ_LOOCV(x, k, mu0, nu0, phi0sq, gamma0))
    }
   
    t.optim1 = proc.time()
    L1 = optim(1, LOOCVfun, lower=0.001, method = "L-BFGS-B")
    t.optim2 = proc.time()
    t.optim = (t.optim2 - t.optim1)[3]
   
    phi0sq = L1$par
   
    ##Now do NNDP
   
    t.NNDP1 = proc.time()
    N1 = uNNDP(x, MC, k, inputpt, mu0, nu0, gamma0, phi0sq, alpha)
    t.NNDP2 = proc.time()
   
    t.NNDP = (t.NNDP2 - t.NNDP1)[3]
   
    f_stor = N1$f_stor
    
    ##KDE
    
    kde_res = kde(x, h = hpi(x), eval.points = inputpt)$estimate
    
    ##DP 
    
    t.DP1 = proc.time()
    dpobj = DirichletProcessGaussian(x)
    #dpobj = UpdateAlpha(dpobj1)
    dpDraw = Fit(dpobj, 2500, updatePrior = FALSE, progressBar = TRUE)
    dpEst = PosteriorFrame(dpDraw, inputpt, ci_size = 0.05, ndraws = 1000)
    t.DP2 = proc.time()
    t.DP = (t.DP2 - t.DP1)[3]
   
    output1 = cbind(inputpt, f0, f_stor, dpEst[,1], dpEst[,2], dpEst[,3], kde_res)
    output2 = c((t.NNDP + t.optim), t.DP)
   
    colnames(output1) = c("x", "f_0(x)", rep("f_hat", MC), "DP_mean", "DP_2.5", "DP_97.5", "KDE_Plugin")
    
    fname1 = paste(dir.name, "/", "run", i, ".csv", sep = "")  ## Add dimension
    fname2 = paste(dir.name, "/", "run", i, "_time", ".txt", sep = "")
   
    write.csv(output1, fname1)
    write.csv(output2, fname2)
  }
  parSapply(cl, 1:n_replicates, NNDP_single_run, n = n, true_den_num = true_den_num)
  stopCluster(cl)
}

for(i in 1:length(dnums))
{
  print(i)
  print(dnums)
  print(sample_sizes[1])
  #print(alpha)
  tfull1 = proc.time()
  NNDP_final_run(sample_sizes[1], n_replicates, dnums[i], n_cores)
  tfull2 = proc.time()
  print((tfull2-tfull1)[3])
  setwd("/home/grad/sc602/NNDP/NNDP Univariate Codes/Updated Codes/Final_Simulations")
}
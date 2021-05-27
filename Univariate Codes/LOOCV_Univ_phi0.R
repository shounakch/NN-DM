##R file for Univariate Leave One Out Cross Validation - need to write Rcpp file.
##Load Faster LOOCV File from Multivariate LOOCV Folder before running this. 

##Load Packages

library(FastKNN)
library(FNN)
library(benchden)
library(dirichletprocess)
library(ks)
library(mvtnorm)
library(tensorA)
library(abind)
library(MCMCpack)
library(CholWishart)
library(expm)
library(Rcpp)
library(matrixStats)

##Criterion function as a function of the hyperparameters and the reduced neighbourhoods

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
  
  # f<-function(l)
  # {
  #   phi0sq = phi0seq[l]
  #   
  #   Tvec = rep(0,n)
  #   i=1
  #   j=1
  #   
  #   for(i in 1:n) #Remove X_{i}#Then exactly j!=i st A[i,j]=1 contains i, others don't
  #   {
  #     condn = which(A[i,]==1)
  #     
  #     dummyvec = rep(0,n)
  #     dummyvec[condn] = x[i]
  #     dummyvec[-condn] = x[grand_nbd[,(k+1)]][-condn]
  #     
  #     nbd.mean.new = ((k+1)*grand_means/k)-(dummyvec/k)
  #     nbd.var.new = ((k+1)*grand_vars/k)-((k+1)*(grand_means-dummyvec)^2/k^2)
  #     
  #     newloc = ((pwt*mu0)+(dwt*nbd.mean.new))[-i]
  #     newphisq = (((gamma0*phi0sq)+(k*nbd.var.new)
  #                  +(k*nu0*(nbd.mean.new-mu0)^2/nun))/gamman)[-i]
  #     
  #     newbw = sqrt(newphisq*(nun+1)/nun)
  #     
  #     Tvec[i] = mean(dt(x = (x[i]-newloc)/newbw,df=gamman)/newbw)
  #     
  #   }
  #   
  #   return(mean(log(Tvec)))
  # }
  
  t1=Sys.time()
  m1 = Tmats_univ(n, L, k, grand_nbd, x, A, phi0seq, grand_means, grand_vars, mu0, nu0, gamma0)
  t2=Sys.time()
  t2-t1
  
  crit=rep(0,L)
  for(l in 1:L)
  {
    crit[l] = mean(log(rowMeans(m1[,,l])))
  }
  
  #crit = sapply(1:L,f)
  
  return(crit)
}

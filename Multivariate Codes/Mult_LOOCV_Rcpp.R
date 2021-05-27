##Multivariate Leave One Out Cross Validation R file
##Call the Faster_LOOCV_Rcpp.cpp file from Rcpp

#Load Packages

library(FastKNN)
library(FNN)
library(dirichletprocess)
library(ks)
library(mvtnorm)
library(Rcpp)
library(plyr)
library(dplyr)

mult.LOOCV<-function(x, k, mu0, nu0, phi0seq, gamma0)
{
  N = dim(as.matrix(x))[1]
  cutoff = 1000
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
### Declare libraries ###

library(FNN)
library(dirichletprocess)
library(MCMCpack)
library(ks)
library(EMMIXcskew)
library(EMMIXskew)
library(mvtnorm)

rnd_sim_densities<-function(n, ntest, dnum)
{
  if(dnum == 1)
  {
    mu1 = c(-2, -2)
    mu2 = c(2, 2)
    rho = 0.8
    Sig = rho*matrix(1, 2, 2) + (1-rho)*diag(2)
    mixture_probs = c(0.4, 0.6)
    X =  rmvnorm.mixt(n, rbind(mu1, mu2) , rbind(Sig, Sig) , props = mixture_probs)
    X_new = rmvnorm.mixt(ntest, rbind(mu1, mu2) , rbind(Sig, Sig) , props = mixture_probs)
    f0 = dmvnorm.mixt(X_new, rbind(mu1, mu2) , rbind(Sig, Sig) , props = mixture_probs)
  }
  
  if(dnum == 2)
  {
    mu = c(0,0)
    rho = 0.8
    Sig = rho*matrix(1, 2, 2) + (1-rho)*diag(2)
    skew_param = c(0.5, 0.5)
    X = rdmsn(n, 2, mu, Sig, skew_param)
    X_new = rdmsn(ntest, 2, mu, Sig, skew_param)
    f0 = ddmsn(X_new, ntest,2, mu, Sig, skew_param)
  }
  
  if(dnum == 3)
  {
    mu = c(1,1)
    rho = 0.8
    Sig = rho*matrix(1, 2, 2) + (1-rho)*diag(2)
    df_param = 10
    X = rmvt(n, Sig, df_param, delta = mu, type = "shifted")
    X_new = rmvt(ntest, Sig, df_param, delta = mu, type = "shifted")
    f0 = dmvt(X_new, mu, Sig, df_param, log = F)
  }
  
  if(dnum == 4)
  {
    mu1 = c(-2, -2)
    mu2 = c(2,2)
    rho = 0.8
    Sig = rho*matrix(1, 2, 2) + (1-rho)*diag(2)
    skew_param = c(0.5, 0.5)
    mixture_probs = c(0.25, 0.75)
    df_param = 10
    X = matrix(0, n, 2)
    for(i in 1:n)
    {
      u = runif(1)
      if(u<mixture_probs[1])
      {
        X[i,] = rdmst(1, 2, mu1, Sig, df_param, skew_param)
      }
      else
      {
        X[i,] = rdmst(1, 2, mu2, Sig, df_param, skew_param)
      }
    }
    X_new = matrix(0, ntest, 2)
    for(i in 1:ntest)
    {
      u = runif(1)
      if(u<mixture_probs[1])
      {
        X_new[i,] = rdmst(1, 2, mu1, Sig, df_param, skew_param)
      }
      else
      {
        X_new[i,] = rdmst(1, 2, mu2, Sig, df_param, skew_param)
      }
    }
    f0 = 0.25*ddmst(X_new, ntest, 2, mu1, Sig, df_param, skew_param) + 0.75*ddmst(X_new,ntest, 2, mu2, Sig, df_param, skew_param)
  }
  
  if(dnum == 5)
  {
    mu1 = c(-2, -2, -2)
    mu2 = c(2, 2, 2)
    rho = 0.8
    Sig = rho*matrix(1, 3, 3) + (1-rho)*diag(3)
    mixture_probs = c(0.4, 0.6)
    X =  rmvnorm.mixt(n, rbind(mu1, mu2) , rbind(Sig, Sig) , props = mixture_probs)
    X_new = rmvnorm.mixt(ntest, rbind(mu1, mu2) , rbind(Sig, Sig) , props = mixture_probs)
    f0 = dmvnorm.mixt(X_new, rbind(mu1, mu2) , rbind(Sig, Sig) , props = mixture_probs)
  }
  
  if(dnum == 6)
  {
    mu = c(0,0, 0)
    rho = 0.8
    Sig = rho*matrix(1, 3, 3) + (1-rho)*diag(3)
    skew_param = c(0.5, 0.5, 0.5)
    X = rdmsn(n, 3, mu, Sig, skew_param)
    X_new = rdmsn(ntest, 3, mu, Sig, skew_param)
    f0 = ddmsn(X_new, ntest, 3, mu, Sig, skew_param)
  }
  
  if(dnum == 7)
  {
    mu = c(1,1, 1)
    rho = 0.8
    Sig = rho*matrix(1, 3, 3) + (1-rho)*diag(3)
    df_param = 10
    X = rmvt(n, Sig, df_param, delta = mu, type = "shifted")
    X_new = rmvt(ntest, Sig, df_param, delta = mu, type = "shifted")
    f0 = dmvt(X_new, mu, Sig, df_param, log = F)
  }
  
  if(dnum == 8)
  {
    mu1 = c(-2, -2, -2)
    mu2 = c(2,2, 2)
    rho = 0.8
    Sig = rho*matrix(1, 3, 3) + (1-rho)*diag(3)
    skew_param = c(0.5, 0.5, 0.5)
    mixture_probs = c(0.25, 0.75)
    df_param = 10
    X = matrix(0, n, 3)
    for(i in 1:n)
    {
      u = runif(1)
      if(u<mixture_probs[1])
      {
        X[i,] = rdmst(1, 3, mu1, Sig, df_param, skew_param)
      }
      else
      {
        X[i,] = rdmst(1, 3, mu2, Sig, df_param, skew_param)
      }
    }
    X_new = matrix(0, ntest, 3)
    for(i in 1:ntest)
    {
      u = runif(1)
      if(u<mixture_probs[1])
      {
        X_new[i,] = rdmst(1, 3, mu1, Sig, df_param, skew_param)
      }
      else
      {
        X_new[i,] = rdmst(1, 3, mu2, Sig, df_param, skew_param)
      }
    }
    f0 = 0.25*ddmst(X_new, ntest, 3, mu1, Sig, df_param, skew_param) + 0.75*ddmst(X_new, ntest, 3, mu2, Sig, df_param, skew_param)
  }
  
  if(dnum == 9)
  {
    mu1 = c(-2, -2, -2, -2)
    mu2 = c(2, 2, 2, 2)
    rho = 0.8
    Sig = rho*matrix(1, 4, 4) + (1-rho)*diag(4)
    mixture_probs = c(0.4, 0.6)
    X =  rmvnorm.mixt(n, rbind(mu1, mu2) , rbind(Sig, Sig) , props = mixture_probs)
    X_new = rmvnorm.mixt(ntest, rbind(mu1, mu2) , rbind(Sig, Sig) , props = mixture_probs)
    f0 =  dmvnorm.mixt(X_new, rbind(mu1, mu2) , rbind(Sig, Sig) , props = mixture_probs)
  }
  
  if(dnum == 10)
  {
    mu = c(0,0, 0, 0)
    rho = 0.8
    Sig = rho*matrix(1, 4, 4) + (1-rho)*diag(4)
    skew_param = c(0.5, 0.5, 0.5, 0.5)
    X = rdmsn(n, 4, mu, Sig, skew_param)
    X_new = rdmsn(ntest, 4, mu, Sig, skew_param)
    f0 = ddmsn(X_new, ntest, 4, mu, Sig, skew_param)
  }
  
  if(dnum == 11)
  {
    mu = c(1,1, 1, 1)
    rho = 0.8
    Sig = rho*matrix(1, 4, 4) + (1-rho)*diag(4)
    df_param = 10
    X = rmvt(n, Sig, df_param, delta = mu, type = "shifted")
    X_new = rmvt(ntest, Sig, df_param, delta = mu, type = "shifted")
    f0 = dmvt(X_new, mu, Sig, df_param, log = F)
  }
  
  if(dnum == 12)
  {
    mu1 = c(-2, -2, -2, -2)
    mu2 = c(2,2, 2, 2)
    rho = 0.8
    Sig = rho*matrix(1, 4, 4) + (1-rho)*diag(4)
    skew_param = c(0.5, 0.5, 0.5, 0.5)
    mixture_probs = c(0.25, 0.75)
    df_param = 10
    X = matrix(0, n, 4)
    for(i in 1:n)
    {
      u = runif(1)
      if(u<mixture_probs[1])
      {
        X[i,] = rdmst(1, 4, mu1, Sig, df_param, skew_param)
      }
      else
      {
        X[i,] = rdmst(1, 4, mu2, Sig, df_param, skew_param)
      }
    }
    X_new = matrix(0, ntest, 4)
    for(i in 1:ntest)
    {
      u = runif(1)
      if(u<mixture_probs[1])
      {
        X_new[i,] = rdmst(1, 4, mu1, Sig, df_param, skew_param)
      }
      else
      {
        X_new[i,] = rdmst(1, 4, mu2, Sig, df_param, skew_param)
      }
    }
    f0 = 0.25*ddmst(X_new, ntest, 4, mu1, Sig, df_param, skew_param) + 0.75*ddmst(X_new, ntest, 4, mu2, Sig, df_param, skew_param)
  }
  
  if(dnum == 13)
  {
    mu1 = c(-2, -2, -2, -2, -2)
    mu2 = c(2, 2, 2, 2, 2)
    rho = 0.8
    Sig = rho*matrix(1, 5, 5) + (1-rho)*diag(5)
    mixture_probs = c(0.4, 0.6)
    X =  rmvnorm.mixt(n, rbind(mu1, mu2) , rbind(Sig, Sig) , props = mixture_probs)
    X_new =  rmvnorm.mixt(ntest, rbind(mu1, mu2) , rbind(Sig, Sig) , props = mixture_probs)
    f0 = dmvnorm.mixt(X_new, rbind(mu1, mu2) , rbind(Sig, Sig) , props = mixture_probs)
  }
  
  if(dnum == 14)
  {
    mu = c(0,0, 0, 0, 0)
    rho = 0.8
    Sig = rho*matrix(1, 5, 5) + (1-rho)*diag(5)
    skew_param = c(0.5, 0.5, 0.5, 0.5, 0.5)
    X = rdmsn(n, 5, mu, Sig, skew_param)
    X_new = rdmsn(ntest, 5, mu, Sig, skew_param)
    f0 = ddmsn(X_new, ntest, 5, mu, Sig, skew_param)
  }
  
  if(dnum == 15)
  {
    mu = c(1,1, 1, 1, 1)
    rho = 0.8
    Sig = rho*matrix(1, 5, 5) + (1-rho)*diag(5)
    df_param = 10
    X = rmvt(n, Sig, df_param, delta = mu, type = "shifted")
    X_new = rmvt(ntest, Sig, df_param, delta = mu, type = "shifted")
    f0 = dmvt(X_new, mu, Sig, df_param, log = F)
  }
  
  if(dnum == 16)
  {
    mu1 = c(-2, -2, -2, -2, -2)
    mu2 = c(2,2, 2, 2, 2)
    rho = 0.8
    Sig = rho*matrix(1, 5, 5) + (1-rho)*diag(5)
    skew_param = c(0.5, 0.5, 0.5, 0.5, 0.5)
    mixture_probs = c(0.25, 0.75)
    df_param = 10
    X = matrix(0, n, 5)
    for(i in 1:n)
    {
      u = runif(1)
      if(u<mixture_probs[1])
      {
        X[i,] = rdmst(1, 5, mu1, Sig, df_param, skew_param)
      }
      else
      {
        X[i,] = rdmst(1, 5, mu2, Sig, df_param, skew_param)
      }
    }
    X_new = matrix(0, ntest, 5)
    for(i in 1:ntest)
    {
      u = runif(1)
      if(u<mixture_probs[1])
      {
        X_new[i,] = rdmst(1, 5, mu1, Sig, df_param, skew_param)
      }
      else
      {
        X_new[i,] = rdmst(1, 5, mu2, Sig, df_param, skew_param)
      }
    }
    f0 = 0.25*ddmst(X_new, ntest, 5, mu1, Sig, df_param, skew_param) + 0.75*ddmst(X_new, ntest, 5, mu2, Sig, df_param, skew_param)
  }
  
  if(dnum == 17)
  {
    mu1 = c(-2, -2, -2, -2, -2, -2)
    mu2 = c(2, 2, 2, 2, 2, 2)
    rho = 0.8
    Sig = rho*matrix(1, 6, 6) + (1-rho)*diag(6)
    mixture_probs = c(0.4, 0.6)
    X =  rmvnorm.mixt(n, rbind(mu1, mu2) , rbind(Sig, Sig) , props = mixture_probs)
    X_new =  rmvnorm.mixt(ntest, rbind(mu1, mu2) , rbind(Sig, Sig) , props = mixture_probs)
    f0 =  dmvnorm.mixt(X_new, rbind(mu1, mu2) , rbind(Sig, Sig) , props = mixture_probs)
  }
  if(dnum == 18)
  {
    mu = c(0,0, 0, 0, 0,0)
    rho = 0.8
    Sig = rho*matrix(1, 6, 6) + (1-rho)*diag(6)
    skew_param = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
    X = rdmsn(n, 6, mu, Sig, skew_param)
    X_new = rdmsn(ntest, 6, mu, Sig, skew_param)
    f0 = ddmsn(X_new, ntest, 6, mu, Sig, skew_param)
  }
  
  if(dnum == 19)
  {
    mu = c(1,1, 1, 1, 1, 1)
    rho = 0.8
    Sig = rho*matrix(1, 6, 6) + (1-rho)*diag(6)
    df_param = 10
    X = rmvt(n, Sig, df_param, delta = mu, type = "shifted")
    X_new = rmvt(ntest, Sig, df_param, delta = mu, type = "shifted")
    f0 = dmvt(X_new, mu, Sig, df_param, log = F)
  }
  
  if(dnum == 20)
  {
    mu1 = c(-2, -2, -2, -2, -2, -2)
    mu2 = c(2,2, 2, 2, 2, 2)
    rho = 0.8
    Sig = rho*matrix(1, 6, 6) + (1-rho)*diag(6)
    skew_param = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
    mixture_probs = c(0.25, 0.75)
    df_param = 10
    X = matrix(0, n, 6)
    for(i in 1:n)
    {
      u = runif(1)
      if(u<mixture_probs[1])
      {
        X[i,] = rdmst(1, 6, mu1, Sig, df_param, skew_param)
      }
      else
      {
        X[i,] = rdmst(1, 6, mu2, Sig, df_param, skew_param)
      }
    }
    X_new = matrix(0, ntest, 6)
    for(i in 1:ntest)
    {
      u = runif(1)
      if(u<mixture_probs[1])
      {
        X_new[i,] = rdmst(1, 6, mu1, Sig, df_param, skew_param)
      }
      else
      {
        X_new[i,] = rdmst(1, 6, mu2, Sig, df_param, skew_param)
      }
    }
    f0 = 0.25*ddmst(X_new, ntest, 6, mu1, Sig, df_param, skew_param) + 0.75*ddmst(X_new, ntest, 6, mu2, Sig, df_param, skew_param)
  }
  
  if(dnum == 21)
  {
    mu1 = rep(-2, 7)
    mu2 = rep(2, 7)
    rho = 0.8
    Sig = rho*matrix(1, 7, 7) + (1-rho)*diag(7)
    mixture_probs = c(0.4, 0.6)
    X =  rmvnorm.mixt(n, rbind(mu1, mu2) , rbind(Sig, Sig) , props = mixture_probs)
    X_new =  rmvnorm.mixt(ntest, rbind(mu1, mu2) , rbind(Sig, Sig) , props = mixture_probs)
    f0 =  dmvnorm.mixt(X_new, rbind(mu1, mu2) , rbind(Sig, Sig) , props = mixture_probs)
  }
  
  if(dnum == 22)
  {
    mu = rep(0, 7)
    rho = 0.8
    Sig = rho*matrix(1, 7, 7) + (1-rho)*diag(7)
    skew_param = rep(0.5, 7)
    X = rdmsn(n, 7, mu, Sig, skew_param)
    X_new = rdmsn(ntest, 7, mu, Sig, skew_param)
    f0 = ddmsn(X_new, ntest, 7, mu, Sig, skew_param)
  }
  
  if(dnum == 23)
  {
    mu = rep(1, 7)
    rho = 0.8
    Sig = rho*matrix(1, 7, 7) + (1-rho)*diag(7)
    df_param = 10
    X = rmvt(n, Sig, df_param, delta = mu, type = "shifted")
    X_new = rmvt(ntest, Sig, df_param, delta = mu, type = "shifted")
    f0 = dmvt(X_new, mu, Sig, df_param, log = F)
  }
  
  if(dnum == 24)
  {
    mu1 = c(-2, -2, -2, -2, -2, -2, - 2)
    mu2 = c(2,2, 2, 2, 2, 2, 2)
    rho = 0.8
    Sig = rho*matrix(1, 7, 7) + (1-rho)*diag(7)
    skew_param = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
    mixture_probs = c(0.25, 0.75)
    df_param = 10
    X = matrix(0, n, 7)
    for(i in 1:n)
    {
      u = runif(1)
      if(u<mixture_probs[1])
      {
        X[i,] = rdmst(1, 7, mu1, Sig, df_param, skew_param)
      }
      else
      {
        X[i,] = rdmst(1, 7, mu2, Sig, df_param, skew_param)
      }
    }
    X_new = matrix(0, ntest, 7)
    for(i in 1:ntest)
    {
      u = runif(1)
      if(u<mixture_probs[1])
      {
        X_new[i,] = rdmst(1, 7, mu1, Sig, df_param, skew_param)
      }
      else
      {
        X_new[i,] = rdmst(1, 7, mu2, Sig, df_param, skew_param)
      }
    }
    f0 = 0.25*ddmst(X_new, ntest, 7, mu1, Sig, df_param, skew_param) + 0.75*ddmst(X_new,ntest, 7, mu2, Sig, df_param, skew_param)
  }
  
  return(list("X" = X, "X_new" = X_new, "f0" = f0))
}

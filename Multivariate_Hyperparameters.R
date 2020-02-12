setwd("~/NNDP/NNDP Multivariate Codes/Updated Codes") #Set correct directory where source files are situated
Rcpp::sourceCpp("Faster_LOOCV_Rcpp.cpp")
Rcpp::sourceCpp("MCSampleRcpp.cpp")
source("Mult_LOOCV_Rcpp.R")
source("MultivariateNNDP_Stock.R")
source("NNDP_multivariate_simulations.R")

n = 1000
#ntest = 500
# dnum = 6
# whole = rnd_sim_densities(n, ntest, dnum)
# x = whole$X
# p = dim(x)[2]
# inputpt = whole$X_new
# trueval = whole$f0

p = 2
x = rmvt(n, delta = rep(0,p), sigma = diag(p), df = 10)
#inputpt = rmvt(ntest, delta = rep(0,p), sigma = diag(p), df = 5)
#trueval = dmvt(inputpt, delta = rep(0,p), sigma = diag(p), log = FALSE, df = 5)

MC = 2000
k = 10

G = 50
x1 = seq(-3,3,length.out=G)
x2 = seq(-3,3,length.out=G)

inputpt = cbind(c(mapply(rep,x1,G)),c(t(mapply(rep,x2,G))))

mu0 = rep(0,p)
nu0 = 0.0001
alpha = 0.1
gamma0 = p-1+1

LOOCVfun<-function(phi0sq)
{
  return(-mult.LOOCV(x, k, mu0, nu0, phi0sq, gamma0))
}

t.optim1 = Sys.time()
L1 = optim(1, LOOCVfun, lower=0.001, method = "L-BFGS-B")
t.optim2 = Sys.time()

t.optim = t.optim2 - t.optim1

phi0sq = L1$par

psi0=diag(rep((gamma0-p+1)*phi0sq,p))

##USE OPTIM INSTEAD!!!!!
##WHY NOT DO CROSS VALIDATION JUST OVER k?

##Do NNDP

t.NNDP1 = Sys.time()
mod1 = mNNDP(x, MC, k, inputpt, mu0, nu0, gamma0, psi0, alpha)
t.NNDP2 = Sys.time()
t.NNDP = t.NNDP2-t.NNDP1

f_stor=mod1$f_stor

t.KDE1 = Sys.time()
f_kde = kde(x, H=Hlscv(x), eval.points = inputpt)$estimate
t.KDE2 = Sys.time()
t.KDE = t.KDE2 - t.KDE1

fNNDPmean = rowMeans(f_stor)
fNNDPmedian = apply(f_stor,1,median)
f_kde = f_kde

mean(abs(trueval-f_kde)/trueval)
mean(abs(trueval-fNNDPmean)/trueval)
mean(abs(trueval-fNNDPmedian)/trueval)

z = NULL
z1 = NULL
z2 = NULL
z3 = NULL
for(i in 1:G)
{
  indices = seq((i-1)*G+1, i*G)
  z = rbind(z, trueval[indices])
  z1 = rbind(z1, f_kde[indices])
  z2 = rbind(z2, fNNDPmean[indices])
  z3 = rbind(z3, fNNDPmedian[indices])
}

persp(x1, x2, z3, theta = 30, phi = 30, expand = 0.5, col = "lightblue", ltheta = 120, shade = 0.75, 
      ticktype = "detailed", xlab = "X1", ylab = "X2", zlab = "f0" )

delta = x1[2]-x1[1]
f_new = rowSums(z2)*delta
f_new_kde = rowSums(z1)*delta

plot(x1, f_new, type="l")
lines(x1, f_new_kde, col="red")
lines(x1, rowSums(z)*delta, col="blue")

est = rowMeans(f_stor)
low = apply(f_stor, 1, quantile, 0.025)
high = apply(f_stor, 1, quantile, 0.975)

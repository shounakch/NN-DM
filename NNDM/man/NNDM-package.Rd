\name{NNDM-package}
\alias{NNDM-package}
\alias{NNDM}
\docType{package}
\title{
\packageTitle{NNDM}
}
\description{
\packageDescription{NNDM}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{NNDM}
\packageIndices{NNDM}
This package has three primary functions.
uNNDM: Provides univariate implementation of the method.
mNNDM: Provides multivariate implementation of the method.
NNDM_postMean: Provides pseudo-posterior mean of the estimator without performing Monte Carlo sampling.
}
\author{
\packageAuthor{NNDM}

Maintainer: \packageMaintainer{NNDM}
}
\references{
Chattopadhyay, Shounak, Antik Chakraborty, and David B. Dunson. "Nearest neighbor Dirichlet process." arXiv preprint arXiv:2003.07953 (2020).
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }
\seealso{
\code{\link[<pkg>:<pkg>-package]{uNNDM}}
\code{\link[<pkg>:<pkg>-package]{mNNDM}}
\code{\link[<pkg>:<pkg>-package]{NNDM_postMean}}
}
\examples{
### Univariate ###

  n = 500
x = rnorm(n)
inputpt = seq(-4, 4, by = 0.01)
f0 = dnorm(inputpt)
res = uNNDM(x = x, inputpt = inputpt) ## phi0sq is chosen by cross validation and alpha is tuned from the data.
fmean = rowMeans(res$f_stor)
plot(inputpt, fmean, type = "l", ylim = c(0, 0.5), main= "Black:Pseudo-posterior mean, Red: True")
lines(inputpt, f0, col = "red")

### Multivariate ###

library(mvtnorm)
n = 500
p = 2
n_t = 200
x = rmvnorm(n, sigma = diag(p))
inputpt = rmvnorm(n_t, sigma = diag(p))
f0 = dmvnorm(inputpt, sigma = diag(p))
k = 10
MC = 1000
mu0 = rep(0, p)
nu0 = 0.001
gamma0 = p
res = mNNDP(x, MC, k, inputpt, mu0, nu0, gamma0) ## psi0 is chosen by cross validation and alpha is tuned from the data.
fmean = rowMeans(res$f_stor)
plot(f0, fmean, main ="Accuracy of the estimate", xlab = "True", ylab = "Estimate")
abline(a=0, b=1)
}

# Nearest Neighbor Dirichlet Mixtures (NN-DM)

There is a rich literature on Bayesian methods for density estimation, which characterize the unknown density as a mixture of kernels. Such methods have advantages in terms of providing uncertainty quantification in estimation, while being adaptive to a rich variety of densities. However, relative to frequentist locally adaptive kernel methods, Bayesian approaches can be slow and unstable to implement in relying on Markov chain Monte Carlo algorithms. To maintain most of the strengths of Bayesian approaches without the computational disadvantages, we propose a class of nearest neighbor-Dirichlet mixtures. The approach starts by grouping the data into neighborhoods based on standard algorithms. Within each neighborhood, the density is characterized via a Bayesian parametric model, such as a Gaussian with unknown parameters. Assigning a Dirichlet prior to the weights on these local kernels, we obtain a pseudo-posterior for the weights and kernel parameters. A simple and embarrassingly parallel Monte Carlo algorithm is proposed to sample from the resulting pseudo-posterior for the unknown density. Desirable asymptotic properties are shown, and the methods are evaluated in simulation studies and applied to a motivating data set in the context of classification.

Arxiv: https://arxiv.org/abs/2003.07953

# Install the package (works for Linux, macOS, and Windows)

First install the 'devtools' package. If using Windows, install Rtools from here: https://cran.r-project.org/bin/windows/Rtools/rtools40.html before installing 'devtools'. If using macOS, install gfortran from here (if not already installed): https://mac.r-project.org/tools/.

```
install.packages("devtools")
```
Before installing the 'NNDM' package, install the following three packages in R:

```
install.packages(c("FNN", "dplyr", "plyr"))
```

Now install the package from GitHub:

```
library(devtools)
devtools::install_github("shounakchattopadhyay/NN-DM")
```
# Using the package

The package has two core functions: 'uNNDM' for univariate density estimation and 'mNNDM' for multivariate density estimation. By default, both the functions carry out cross-validation to choose the neighborhood bandwidth matrix hyperparameter, and then proceeds to perform Monte Carlo sampling. One can obtain the pseudo-posterior mean directly without performing Monte Carlo sampling by setting 'samples=FALSE' in the function. Examples are given below.

## Univariate Example

Here we consider data generated from N(0,1). 

```
### Univariate Example ###

library(NNDM)

n = 500
x = rnorm(n)
inputpt = seq(-4, 4, by = 0.01)
f0 = dnorm(inputpt)
res = uNNDM(x = x, inputpt = inputpt) ## phi0sq is chosen by cross validation and alpha is tuned from the data.
fmean = rowMeans(res$f_stor)
plot(inputpt, fmean, type = "l", ylim = c(0, 0.5), main= "Black: Pseudo-posterior mean, Red: True")
lines(inputpt, f0, col = "red")

```

## Multivariate Example

Here we consider data from the bivariate standard Gaussian density.
```
### Multivariate Example ###

library(mvtnorm)
library(NNDM)

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
res = mNNDM(x = x, inputpt = inputpt, mu0 = mu0, nu0 = nu0, gamma0 = gamma0) ## psi0 is chosen by cross validation and alpha is tuned from the data.
fmean = rowMeans(res$f_stor)
plot(f0, fmean, main = "Accuracy of the estimate", xlab = "True", ylab = "Estimate")
abline(a=0, b=1)
```
## Pseudo-Posterior Mean

The function 'NNDM_postMean' obtains the pseudo-posterior mean without carrying out Monte Carlo sampling. All the hyperparameters must be supplied to this function.

```
x = cbind(rnorm(200), rnorm(200))
k = 10
inputpt = cbind(rnorm(500), rnorm(500))
f0 = dnorm(inputpt[,1])*dnorm(inputpt[,2])
mu0 = rep(0,dim(x)[2])
nu0 = 0.001
gamma0 = dim(x)[2]
delta0sq = 1
psi0 = ((gamma0-dim(x)[2]+1) * delta0sq) * diag(dim(x)[2])

f_hat = NNDM_postMean(x, k, inputpt, mu0, nu0, gamma0, psi0)
```



\name{mNNDM}
\alias{mNNDM}
\title{Multivariate implementation of the Nearest Neighbor Dirichlet Mixture}

\description{
  The function implements the fast pseudo-Bayesian density estimator proposed in Chattopadhyay et al. (2021+) for multivariate data. The method first creates nearest neighborhoods for each data point and fits a Gaussian kernel within each neighborhood. Priors are specified for parameters for each neighborhood which are then updated using data within that neighborhood. A final density estimate for a new point is formed by combining Gaussian density values of this query point for each of neighborhoods with updated paramaters using Dirichlet weights. 
}

\usage{
mNNDM(x, MC = 1000, k = 10, inputpt, mu0, nu0, gamma0, psi0 = "CV", alpha = "RR", samples = TRUE)
}

\arguments{
  \item{x}{Multivariate data on which model is fitted. Is a n by p matrix where n = number of samples and p = dimension of each sample.}
  \item{MC}{Number of Monte Carlo samples to be drawn from the pseudo-posterior}
  \item{k}{Number of nearest neighbors for forming the neighborhoods. Default value is k = 10.}
  \item{inputpt}{Query points at which density is to be estimated. Is a n_t by p matrix by n_t = number of query points.}
  \item{mu0}{Prior mean of neighborhood mean parameter in each neighborhood. Typically taken to be mu0 = 0_p, the vector of 0s in p dimensions.}
  \item{nu0}{Prior variance coefficient of neighborhood mean parameter. Typically set to nu0 = 0.001.}
  \item{gamma0}{Shape parameter of inverse Gamma prior on neighborhood variance parameter. Typically set to gamma0 = p. Must be greater than p-1.}
  \item{psi0}{Scale parameter of inverse Wishart prior on neighborhood variance parameter. Default choice is to choose psi0 by leave-one-out cross validation, denoted by "CV". User may provide their own choice of psi0. We recommend the default value psi0 = I_p, the p by p identity matrix if not using cross validation. Must be positive definite.}
  \item{alpha}{Prior concentration parameter of the Dirichlet weights. Default choice is to use alpha = "RR" which uses Rivoirard and Rousseau (2012). User may provide their own choice of alpha. Must be positive.}
  \item{samples}{Logical value set to TRUE by default, which draws pseudo-posterior samples using MC. If FALSE, returns only the pseudo-posterior mean.}
}

\details{
  A list is returned with the following components.
}

\value{
  \item{f_stor}{n_t by MC matrix containing MC number of pseudo-posterior samples for density at each query point in inputpt. Not returned if samples=FALSE.}
  \item{f_hat}{n_t by 1 vector containing the pseudo-posterior mean of fitted density estimate samples at each query point. Not returned if samples=TRUE.}
  \item{k}{Returns the choice of the number of nearest neighbors provided.}
  \item{n}{Returns the number of training samples used.}
  \item{MC}{Number of Monte Carlo iterations performed.}
  \item{inputpt}{Returns the matrix of query points at which density estimation is performed.}
  \item{data}{Returns the training data used to fit the model.}
  \item{hyp_par}{Returns a list of hyperparameters containing the values of mu0, nu0, gamma0, psi0 and alpha used in fitting the model.}
  \item{time}{Returns the runtime of several components in fitting the model. First component provides the runtime for k-NN allocation and second component provides the runtime for MC sampling.}
  \item{nbd_stat}{Returns a list of the neighborhood sample statistics used to fit the model. Consists of nbd_mean (sample neighborhood mean), nbd_var (sample neighborhood variance), nbd_range (sample neighborhood radius), psi (Sample neighborhood summary used to generate samples of neighborhood variance parameter) and nbd_wt (neighborhood weights in pseudo-posterior mean).}
}

\references{
  Chattopadhyay, S., Chakraborty, A., & Dunson, D. B. (2020). Nearest Neighbor Dirichlet Mixtures. arXiv preprint arXiv:2003.07953.
  
  Rivoirard, V., & Rousseau, J. (2012). Bernstein–von Mises theorem for linear functionals of the density. The Annals of Statistics, 40(3), 1489-1523. 
}

\author{Shounak Chattopadhyay}

\examples{
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
res = mNNDM(x = x, inputpt = inputpt, mu0 = mu0, nu0 = nu0, gamma0 = gamma0) ## psi0 is chosen by cross validation and alpha is tuned from the data.
fmean = rowMeans(res$f_stor)
plot(f0, fmean, main ="Accuracy of the estimate", xlab = "True", ylab = "Estimate")
abline(a=0, b=1)
}



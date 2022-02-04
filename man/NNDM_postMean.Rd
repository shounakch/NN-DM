\name{NNDM_postMean}
\alias{NNDM_postMean}
\title{Pseudo-posterior mean of the Nearest Neighbor Dirichlet Mixture}

\description{
  The function obtains the pseudo-posterior mean of the fast pseudo-Bayesian density estimator proposed in Chattopadhyay et al. (2021+) without performing any Monte Carlo sampling. Works for both univariate and multivariate data.
}

\usage{
NNDM_postMean(x, k, inputpt, mu0, nu0, gamma0, psi0) 
}

\arguments{
  \item{x}{Data on which model is fitted. Is a n by p matrix where n = number of samples and p = dimension of each sample.}
  \item{k}{Number of nearest neighbors for forming the neighborhoods.}
  \item{inputpt}{Query points at which density is to be estimated. Is a n_t by p matrix by n_t = number of query points.}
  \item{mu0}{Prior mean of neighborhood mean parameter in each neighborhood.}
  \item{nu0}{Prior variance coefficient of neighborhood mean parameter.}
  \item{gamma0}{Shape parameter of inverse Wishart prior on neighborhood variance parameter. Must be greater than p-1.}
  \item{psi0}{Scale parameter of inverse Wishart prior on neighborhood variance parameter. Must be positive definite.}
}

\details{
  A list is returned with the following components.
}

\value{
  \item{f_hat}{n_t by 1 vector containing the pseudo-posterior mean of fitted density estimate samples at each query point.}
  \item{k}{Returns the choice of the number of nearest neighbors provided.}
  \item{n}{Returns the number of training samples used.}
  \item{inputpt}{Returns the matrix of query points at which density estimation is performed.}
  \item{data}{Returns the training data used to fit the model.}
  \item{hyp_par}{Returns a list of hyperparameters containing the values of mu0, nu0, gamma0, and psi0 used in fitting the model.}
}

\references{
  Chattopadhyay, Shounak, Antik Chakraborty, and David B. Dunson. "Nearest neighbor Dirichlet process." arXiv preprint arXiv:2003.07953 (2020).
}

\author{Shounak Chattopadhyay}

\examples{
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
}



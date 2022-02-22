# Nearest Neighbor Dirichlet Mixtures (NN-DM)

There is a rich literature on Bayesian methods for density estimation, which characterize the unknown density as a mixture of kernels. Such methods have advantages in terms of providing uncertainty quantification in estimation, while being adaptive to a rich variety of densities. However, relative to frequentist locally adaptive kernel methods, Bayesian approaches can be slow and unstable to implement in relying on Markov chain Monte Carlo algorithms. To maintain most of the strengths of Bayesian approaches without the computational disadvantages, we propose a class of nearest neighbor-Dirichlet mixtures. The approach starts by grouping the data into neighborhoods based on standard algorithms. Within each neighborhood, the density is characterized via a Bayesian parametric model, such as a Gaussian with unknown parameters. Assigning a Dirichlet prior to the weights on these local kernels, we obtain a pseudo-posterior for the weights and kernel parameters. A simple and embarrassingly parallel Monte Carlo algorithm is proposed to sample from the resulting pseudo-posterior for the unknown density. Desirable asymptotic properties are shown, and the methods are evaluated in simulation studies and applied to a motivating data set in the context of classification.

Arxiv: https://arxiv.org/abs/2003.07953

# Install the package (works for Linux, macOS, and Windows)

First install the 'devtools' package. If using Windows, install Rtools from here: https://cran.r-project.org/bin/windows/Rtools/rtools40.html before installing 'devtools'.

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

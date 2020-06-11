# GCPBayes package


> Run a Gibbs sampler for a multivariate Bayesian sparse group selection model with Dirac, continoues and hierarchical spike prior for detecting pleiotropic effects on two traits. This package is designed for summary statistics containing estimated regression coefficients and its estimated covariance matrix.

> Here we provide a comprehensible detailed vignette to run Bayesian meta-analysis models using GCPBayes package.

> We illustrate the inference by some simulated summary statistics.

> This vignette can reproduces a part of the simulation study of our submitted paper entitled as follows:

  > Bayesian meta-analysis models to Cross Cancer Genomic Investigation of pleiotropic effects using group structure. Baghfalaki T., P.E. Sugier, T. Truong, A.T. Pettitt, . Mengersen, and B. Liquet (2020).



**Getting Started**

  > In this vignette we used the following R packages
- library(MASS)
- library(mvtnorm)
- library(invgamma)
- library(wiqid)
- library(gdata)
- library(truncnorm)
- library(bayesplot)
- library(abind)
- library(magrittr)


**Example Usage**
  > Examples from the GCPaBayes paper can be replicated by following the process described in

- CS: This analysis is presented [here](/CSExamp.md)
- DS: This analysis is presented [here](/DSExamp.md)
- HS: This analysis is presented [here](/HSExamp.md)

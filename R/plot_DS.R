#'  Trace plot for the output of DS
#'
#' @param DS All the generated results by DS function.
#' @param kappa logical; if TRUE, the trace plot of kappa is drawan.
#' @param sigma2 logical; if TRUE, the trace plot of sigma2 is drawan.
#' @param Beta1 A vector of indexes for drawing trace plot of Beta1. If NULL, the trace plot of Beta1 is not drawn.
#' @param Beta2 A vector of indexes for drawing trace plot of Beta2. If NULL, the trace plot of Beta2 is not drawn.
#' @param nrow1 Number of rows of the trace plot of Beta1.
#' @param nrow2 Number of rows of the trace plot of Beta2.
#' @example inst/examplePlotDStrace.R
#' @export
plot_trace_DS = function(DS, kappa=FALSE, sigma2=FALSE,
                         Beta1=NULL, Beta2=NULL,
                         nrow1 = 1, nrow2 = 1) {
  #require(bayesplot)
  #require(abind)
  #require(magrittr)
  #We need to reshape data to look like what bayesplot wants
  reshaped = reshape_DS(DS)
  mcmc_kappa = reshaped$mcmc_kappa
  mcmc_sigma2 = reshaped$mcmc_sigma2
  mcmc_Beta1 = reshaped$mcmc_Beta1
  mcmc_Beta2 = reshaped$mcmc_Beta2

  bayesplot::color_scheme_set("viridis")
  if (kappa) {
    p1 = bayesplot::mcmc_trace(mcmc_kappa, "kappa")
    print(p1)
  }
  if (sigma2) {
    p2 = bayesplot::mcmc_trace(mcmc_sigma2, "sigma2")
    print(p2)
  }
  if (!is.null(Beta1)) {
    pars1 = paste(Beta1)
    p3 = bayesplot::mcmc_trace(mcmc_Beta1, facet_args = list(nrow=nrow1))
    print(p3)
  }
  if (!is.null(Beta2)) {
    pars2 = paste(Beta2)
    p4 = bayesplot::mcmc_trace(mcmc_Beta2, facet_args = list(nrow=nrow2))
    print(p4)
  }
}



#' ACF plot for the output of DS
#'
#' @param DS All the generated results by DS function.
#' @param kappa logical; if TRUE, the ACF plot of kappa is drawan.
#' @param sigma2 logical; if TRUE, the ACF plot of sigma2 is drawan.
#' @param Beta1 A vector of indexes for drawing ACF plot of Beta1. If NULL, the ACF plot of Beta1 is not drawn.
#' @param Beta2 A vector of indexes for drawing ACF plot of Beta2. If NULL, the ACF plot of Beta2 is not drawn.
#' @param nrow1 Number of rows of the ACF plot of Beta1.
#' @param nrow2 Number of rows of the ACF plot of Beta2.
#' @example inst/examplePlotDSacf.R
#' @export
plot_acf_DS = function(DS, kappa=FALSE, sigma2=FALSE,
                       Beta1=NULL, Beta2=NULL,
                       nrow1 = 1, nrow2 = 1) {

  # We need to reshape data to look like what bayesplot wants
  reshaped = reshape_DS(DS)
  mcmc_kappa = reshaped$mcmc_kappa
  mcmc_sigma2 = reshaped$mcmc_sigma2
  mcmc_Beta1 = reshaped$mcmc_Beta1
  mcmc_Beta2 = reshaped$mcmc_Beta2

  if (kappa) {
    p1 = bayesplot::mcmc_acf(mcmc_kappa, "kappa")
    print(p1)
  }
  if (sigma2) {
    p2 = bayesplot::mcmc_acf(mcmc_sigma2, "sigma2")
    print(p2)
  }
  if (!is.null(Beta1)) {
    pars1 = paste(Beta1)
    p3 = bayesplot::mcmc_acf(mcmc_Beta1)
    print(p3)
  }
  if (!is.null(Beta2)) {
    pars2 = paste(Beta2)
    p4 = bayesplot::mcmc_acf(mcmc_Beta2)
    print(p4)
  }
}


#' Density plot for the output of DS
#'
#' @param DS All the generated results by DS function.
#' @param kappa logical; if TRUE, the density plot of kappa is drawan.
#' @param sigma2 logical; if TRUE, the density plot of sigma2 is drawan.
#' @param Beta1 A vector of indexes for drawing density plot of Beta1. If NULL, the density plot of Beta1 is not drawn.
#' @param Beta2 A vector of indexes for drawing density plot of Beta2. If NULL, the density plot of Beta2 is not drawn.
#' @param nrow1 Number of rows of the density plot of Beta1.
#' @param nrow2 Number of rows of the density plot of Beta2.
#' @example inst/examplePlotDSdensity.R
#' @export
plot_density_DS = function(DS, kappa=FALSE, sigma2=FALSE,
                           Beta1=NULL, Beta2=NULL,
                           nrow1 = 1, nrow2 = 1) {

  # We need to reshape data to look like what bayesplot wants
  reshaped = reshape_DS(DS)
  mcmc_kappa = reshaped$mcmc_kappa
  mcmc_sigma2 = reshaped$mcmc_sigma2
  mcmc_Beta1 = reshaped$mcmc_Beta1
  mcmc_Beta2 = reshaped$mcmc_Beta2

  nchains = dim(mcmc_kappa)[2]

  bayesplot::color_scheme_set("mix-teal-pink")
  if (kappa) {
    if (nchains==1) {
      p1 = bayesplot::mcmc_dens(mcmc_kappa, "kappa")
      print(p1)
    } else {
      p1 = bayesplot::mcmc_dens_overlay(mcmc_kappa, "kappa")
      print(p1)
    }
  }
  if (sigma2) {
    if (nchains==1) {
      p2 = bayesplot::mcmc_dens(mcmc_sigma2, "sigma2")
      print(p2)
    } else {
      p2 = bayesplot::mcmc_dens_overlay(mcmc_sigma2, "sigma2")
      print(p2)
    }
  }
  if (!is.null(Beta1)) {
    pars1 = paste(Beta1)
    if (nchains==1) {
      p3 = bayesplot::mcmc_dens(mcmc_Beta1, facet_args = list(nrow=nrow1))
      print(p3)
    } else {
      p3 = bayesplot::mcmc_dens_chains(mcmc_Beta1)
      print(p3)
    }
  }
  if (!is.null(Beta2)) {
    pars2 = paste(Beta2)
    if (nchains==1) {
      p4 = bayesplot::mcmc_dens(mcmc_Beta2, facet_args = list(nrow=nrow2))
      print(p4)
    } else {
      p4 = bayesplot::mcmc_dens_chains(mcmc_Beta2)
      print(p4)
    }
  }
}


#' Internal : Reshape DS mcmc chain
#'
#' @param DS All the generated results by DS function.
#' @export
reshape_DS = function(DS) {
  # We need to reshape data to look like what bayesplot wants

  # If only 1 chain
  if (length(DS)==1) {
    # Just reshape data
    chains = DS[[1]]$mcmcchain
    mcmc_kappa = array(chains$kappa, dim=c(length(chains$kappa),1, 1), dimnames = list(NULL, NULL, "kappa"))
    mcmc_sigma2 = array(chains$sigma2, dim=c(length(chains$sigma2),1,1), dimnames = list(NULL, NULL, "sigma2"))
    mcmc_Beta1 = array(chains$Beta1, dim=c(dim(chains$Beta1)[1],1,dim(chains$Beta1)[2]))
    mcmc_Beta2 = array(chains$Beta2, dim=c(dim(chains$Beta2)[1],1,dim(chains$Beta2)[2]))

    dimnames(mcmc_Beta1)[[3]] = paste0("Beta1[",seq(dim(mcmc_Beta1)[3]),"]")
    dimnames(mcmc_Beta2)[[3]] = paste0("Beta2[",seq(dim(mcmc_Beta2)[3]),"]")

    dimnames(mcmc_kappa)[[2]] = "chain1"
    dimnames(mcmc_sigma2)[[2]] = "chain1"
    dimnames(mcmc_Beta1)[[2]] = "chain1"
    dimnames(mcmc_Beta2)[[2]] = "chain1"

    dimnames(mcmc_kappa)[[1]] = seq(dim(mcmc_kappa)[1])
    dimnames(mcmc_sigma2)[[1]] = seq(dim(mcmc_sigma2)[1])
    dimnames(mcmc_Beta1)[[1]] = seq(dim(mcmc_Beta1)[1])
    dimnames(mcmc_Beta2)[[1]] = seq(dim(mcmc_Beta2)[1])

    names(dimnames(mcmc_kappa)) = c("iter", "chain", "param")
    names(dimnames(mcmc_sigma2)) = c("iter", "chain", "param")
    names(dimnames(mcmc_Beta1)) = c("iter", "chain", "param")
    names(dimnames(mcmc_Beta2)) = c("iter", "chain", "param")

  # If more than 1 chain
  } else if (length(DS)>1) {
    # Need to combine all chains together
    outputs = DS$Outputs
    mcmc_kappa = lapply(outputs, function(x) {
      chain = x$mcmcchain$kappa
      array(chain, dim = c(length(chain),1, 1))
    }) %>% abind::abind(., along=2)

    mcmc_sigma2 = lapply(outputs, function(x) {
      chain = x$mcmcchain$sigma2
      array(chain, dim = c(length(chain),1,1))
    }) %>% abind::abind(., along=2)

    mcmc_Beta1 = lapply(outputs, function(x) {
      chain = x$mcmcchain$Beta1
      array(chain, dim = c(dim(chain)[1],1,dim(chain)[2]))
    }) %>% abind::abind(., along=2)

    mcmc_Beta2 = lapply(outputs, function(x) {
      chain = x$mcmcchain$Beta2
      array(chain, dim = c(dim(chain)[1],1, dim(chain)[2]))
    }) %>% abind::abind(., along=2)


    dimnames(mcmc_kappa)[[3]] = "kappa"
    dimnames(mcmc_sigma2)[[3]] = "sigma2"
    dimnames(mcmc_Beta1)[[3]] = paste0("Beta1[",seq(dim(mcmc_Beta1)[3]),"]")
    dimnames(mcmc_Beta2)[[3]] = paste0("Beta2[",seq(dim(mcmc_Beta2)[3]),"]")

    dimnames(mcmc_kappa)[[2]] = paste0("chain", seq(dim(mcmc_kappa)[2]))
    dimnames(mcmc_sigma2)[[2]] = paste0("chain", seq(dim(mcmc_sigma2)[2]))
    dimnames(mcmc_Beta1)[[2]] = paste0("chain", seq(dim(mcmc_Beta1)[2]))
    dimnames(mcmc_Beta2)[[2]] = paste0("chain", seq(dim(mcmc_Beta2)[2]))

    dimnames(mcmc_kappa)[[1]] = seq(dim(mcmc_kappa)[1])
    dimnames(mcmc_sigma2)[[1]] = seq(dim(mcmc_sigma2)[1])
    dimnames(mcmc_Beta1)[[1]] = seq(dim(mcmc_Beta1)[1])
    dimnames(mcmc_Beta2)[[1]] = seq(dim(mcmc_Beta2)[1])

    names(dimnames(mcmc_kappa)) = c("iter", "chain", "param")
    names(dimnames(mcmc_sigma2)) = c("iter", "chain", "param")
    names(dimnames(mcmc_Beta1)) = c("iter", "chain", "param")
    names(dimnames(mcmc_Beta2)) = c("iter", "chain", "param")
  }

  out = list(
    mcmc_kappa=mcmc_kappa ,
    mcmc_sigma2=mcmc_sigma2,
    mcmc_Beta1=mcmc_Beta1,
    mcmc_Beta2=mcmc_Beta2
  )
  return(out)
}

#'  Trace plot for the output of CS
#'
#' @param CS All the generated results by CS function.
#' @param kappa logical; if TRUE, the trace plot of kappa is drawan.
#' @param tau2 logical; if TRUE, the trace plot of tau2 is drawan.
#' @param Beta1 A vector of indexes for drawing trace plot of Beta1. If NULL, the trace plot of Beta1 is not drawn.
#' @param Beta2 A vector of indexes for drawing trace plot of Beta2. If NULL, the trace plot of Beta2 is not drawn.
#' @param nrow1 Number of rows of the trace plot of Beta1.
#' @param nrow2 Number of rows of the trace plot of Beta2.
#' @example inst/examplePlotCStrace.R
#' @export
plot_trace_CS = function(CS, kappa=FALSE, tau2=FALSE,
                         Beta1=NULL, Beta2=NULL,
                         nrow1 = 1, nrow2 = 1) {
  # We need to reshape data to look like what bayesplot wants
  reshaped = reshape_CS(CS)
  mcmc_kappa = reshaped$mcmc_kappa
  mcmc_tau2 = reshaped$mcmc_tau2
  mcmc_Beta1 = reshaped$mcmc_Beta1
  mcmc_Beta2 = reshaped$mcmc_Beta2

  bayesplot::color_scheme_set("viridis")
  if (kappa) {
    p1 = bayesplot::mcmc_trace(mcmc_kappa, "kappa")
    print(p1)
  }
  if (tau2) {
    p2 = bayesplot::mcmc_trace(mcmc_tau2, "tau2")
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



#' ACF plot for the output of CS
#'
#' @param CS All the generated results by CS function.
#' @param kappa logical; if TRUE, the ACF plot of kappa is drawan.
#' @param tau2 logical; if TRUE, the ACF plot of tau2 is drawan.
#' @param Beta1 A vector of indexes for drawing ACF plot of Beta1. If NULL, the ACF plot of Beta1 is not drawn.
#' @param Beta2 A vector of indexes for drawing ACF plot of Beta2. If NULL, the ACF plot of Beta2 is not drawn.
#' @param nrow1 Number of rows of the ACF plot of Beta1.
#' @param nrow2 Number of rows of the ACF plot of Beta2.
#' @example inst/examplePlotCSacf.R
#' @export
plot_acf_CS = function(CS, kappa=FALSE, tau2=FALSE,
                       Beta1=NULL, Beta2=NULL,
                       nrow1 = 1, nrow2 = 1) {

  # We need to reshape data to look like what bayesplot wants
  reshaped = reshape_CS(CS)
  mcmc_kappa = reshaped$mcmc_kappa
  mcmc_tau2 = reshaped$mcmc_tau2
  mcmc_Beta1 = reshaped$mcmc_Beta1
  mcmc_Beta2 = reshaped$mcmc_Beta2

  if (kappa) {
    p1 = bayesplot::mcmc_acf(mcmc_kappa, "kappa")
    print(p1)
  }
  if (tau2) {
    p2 = bayesplot::mcmc_acf(mcmc_tau2, "tau2")
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


#' Density plot for the output of CS
#'
#' @param CS All the generated results by CS function.
#' @param kappa logical; if TRUE, the density plot of kappa is drawan.
#' @param tau2 logical; if TRUE, the density plot of tau2 is drawan.
#' @param Beta1 A vector of indexes for drawing density plot of Beta1. If NULL, the density plot of Beta1 is not drawn.
#' @param Beta2 A vector of indexes for drawing density plot of Beta2. If NULL, the density plot of Beta2 is not drawn.
#' @param nrow1 Number of rows of the density plot of Beta1.
#' @param nrow2 Number of rows of the density plot of Beta2.
#' @example inst/examplePlotCSdensity.R
#' @export
plot_density_CS = function(CS, kappa=FALSE, tau2=FALSE,
                           Beta1=NULL, Beta2=NULL,
                           nrow1 = 1, nrow2 = 1) {
  # We need to reshape data to look like what bayesplot wants
  reshaped = reshape_CS(CS)
  mcmc_kappa = reshaped$mcmc_kappa
  mcmc_tau2 = reshaped$mcmc_tau2
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
  if (tau2) {
    if (nchains==1) {
      p2 = bayesplot::mcmc_dens(mcmc_tau2, "tau2")
      print(p2)
    } else {
      p2 = bayesplot::mcmc_dens_overlay(mcmc_tau2, "tau2")
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


#' Internal : Reshape CS mcmc chain
#'
#' @param CS All the generated results by CS function.
#' @export
reshape_CS = function(CS) {
  # We need to reshape data to look like what bayesplot wants

  # If only 1 chain
  if (length(CS)==1) {
    # Just reshape data
    chains = CS[[1]]$mcmcchain
    mcmc_kappa = array(chains$kappa, dim=c(length(chains$kappa),1, 1), dimnames = list(NULL, NULL, "kappa"))
    mcmc_tau2 = array(chains$tau2, dim=c(length(chains$tau2),1,1), dimnames = list(NULL, NULL, "tau2"))
    mcmc_Beta1 = array(chains$Beta1, dim=c(dim(chains$Beta1)[1],1,dim(chains$Beta1)[2]))
    mcmc_Beta2 = array(chains$Beta2, dim=c(dim(chains$Beta2)[1],1,dim(chains$Beta2)[2]))

    dimnames(mcmc_Beta1)[[3]] = paste0("Beta1[",seq(dim(mcmc_Beta1)[3]),"]")
    dimnames(mcmc_Beta2)[[3]] = paste0("Beta2[",seq(dim(mcmc_Beta2)[3]),"]")

    dimnames(mcmc_kappa)[[2]] = "chain1"
    dimnames(mcmc_tau2)[[2]] = "chain1"
    dimnames(mcmc_Beta1)[[2]] = "chain1"
    dimnames(mcmc_Beta2)[[2]] = "chain1"

    dimnames(mcmc_kappa)[[1]] = seq(dim(mcmc_kappa)[1])
    dimnames(mcmc_tau2)[[1]] = seq(dim(mcmc_tau2)[1])
    dimnames(mcmc_Beta1)[[1]] = seq(dim(mcmc_Beta1)[1])
    dimnames(mcmc_Beta2)[[1]] = seq(dim(mcmc_Beta2)[1])

    names(dimnames(mcmc_kappa)) = c("iter", "chain", "param")
    names(dimnames(mcmc_tau2)) = c("iter", "chain", "param")
    names(dimnames(mcmc_Beta1)) = c("iter", "chain", "param")
    names(dimnames(mcmc_Beta2)) = c("iter", "chain", "param")

    # If more than 1 chain
  } else if (length(CS)>1) {
    # Need to combine all chains together
    outputs = CS$Outputs
    mcmc_kappa = lapply(outputs, function(x) {
      chain = x$mcmcchain$kappa
      array(chain, dim = c(length(chain),1, 1))
    }) %>% abind::abind(., along=2)

    mcmc_tau2 = lapply(outputs, function(x) {
      chain = x$mcmcchain$tau2
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
    dimnames(mcmc_tau2)[[3]] = "tau2"
    dimnames(mcmc_Beta1)[[3]] = paste0("Beta1[",seq(dim(mcmc_Beta1)[3]),"]")
    dimnames(mcmc_Beta2)[[3]] = paste0("Beta2[",seq(dim(mcmc_Beta2)[3]),"]")

    dimnames(mcmc_kappa)[[2]] = paste0("chain", seq(dim(mcmc_kappa)[2]))
    dimnames(mcmc_tau2)[[2]] = paste0("chain", seq(dim(mcmc_tau2)[2]))
    dimnames(mcmc_Beta1)[[2]] = paste0("chain", seq(dim(mcmc_Beta1)[2]))
    dimnames(mcmc_Beta2)[[2]] = paste0("chain", seq(dim(mcmc_Beta2)[2]))

    dimnames(mcmc_kappa)[[1]] = seq(dim(mcmc_kappa)[1])
    dimnames(mcmc_tau2)[[1]] = seq(dim(mcmc_tau2)[1])
    dimnames(mcmc_Beta1)[[1]] = seq(dim(mcmc_Beta1)[1])
    dimnames(mcmc_Beta2)[[1]] = seq(dim(mcmc_Beta2)[1])

    names(dimnames(mcmc_kappa)) = c("iter", "chain", "param")
    names(dimnames(mcmc_tau2)) = c("iter", "chain", "param")
    names(dimnames(mcmc_Beta1)) = c("iter", "chain", "param")
    names(dimnames(mcmc_Beta2)) = c("iter", "chain", "param")
  }

  out = list(
    mcmc_kappa=mcmc_kappa ,
    mcmc_tau2=mcmc_tau2,
    mcmc_Beta1=mcmc_Beta1,
    mcmc_Beta2=mcmc_Beta2
  )
  return(out)
}

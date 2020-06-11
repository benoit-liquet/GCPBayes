#' Continuous Spike
#'
#'
#' @description
#' Run a Gibbs sampler for a multivariate Bayesian sparse group selection model with continuous spike prior for detecting pleiotropic effects on two traits. This function is designed for summary statistics containing estimated regression coefficients and their estimated covariance matrices.
#'
#'
#' @details
#' Run a Gibbs sampler using a continuous spike.
#'
#'
#' @param Betah1 A numerical vector of length mg representing the regression coefficients for the first trait.
#' @param Betah2 A numerical vector of length mg representing the regression coefficients for the second trait.
#' @param Sigmah1 A mg*mg positive definite covariance matrix where it is estimated covariance matrix of Betah1.
#' @param Sigmah2 A mg*mg positive definite covariance matrix where it is estimated covariance matrix of Betah2.
#' @param kappa0 Initial value for kappa (its dimension is equal to nchains).
#' @param tau20 Initial value for tau2 (its dimension is equal to nchains).
#' @param zeta10 Initial value for zeta1 (its dimension is equal to nchains).
#' @param zeta20 Initial value for zeta2 (its dimension is equal to nchains).
#' @param mg Number of variables in the group.
#' @param niter Number of iterations for the Gibbs sampler.
#' @param burnin Number of burn-in iterations.
#' @param nthin The lag of the iterations used for the posterior analysis is defined (or thinning rate).
#' @param nchains Number of Markov chains, when nchains>1, the function calculates the Gelman-Rubin convergence statistic, as modified by Brooks and Gelman (1998).
#' @param a1,a2 Hyperparameters of kappa. Default is a1=0.1 and a2=1.
#' @param c1,c2 Hyperparameters of tau2. Default is c1=0.1 and c2=1.
#' @param sigma2 Variance of spike (multivariate normal distribution with a diagonal covariance matrix with small variance) representing the null effect distribution. Default is 10^-3.
#' @param snpnames Names of variables for the group.
#' @param genename Name of group.
#'
#'
#' @return
#' - mcmcchain: The list of simulation output for all parameters.
#' - Criteria: genename, snpnames, log10BF, lBFDR, theta, PPA1, PPA2 and the number of studies with nonzero signal by credible interval (CI).
#' - Statistics of Trait 1 for Beta_1: Summary statistics including Mean, SD, val2.5pc, Median,	val97.5pc and BGR (if nchains>1).
#' - Statistics of Trait 2  for Beta_2: Summary statistics including Mean, SD, val2.5pc, Median,	val97.5pc and BGR (if nchains>1).
#' - Other Parameters: Summary statistics for kappa and tau2 including Mean, SD, val2.5pc, Median,	val97.5pc and BGR (if nchains>1).
#'
#' @author Taban Baghfalaki.
#'
#' @references
#' T. Baghfalaki, P.E. Sugier, T. Truong, A.T. Pettitt, K. Mengersen and B. Liquet. (2020). Bayesian meta-analysis models to Cross Cancer Genomic Investigation of pleiotropic effects using group structure. *Submitted in Statistics in Medicine*.
#'
#' @example inst/exampleCS.R
#'
#' @md
#' @export


CS=function(Betah1, Betah2, Sigmah1, Sigmah2, kappa0, tau20, zeta10, zeta20, mg, niter=1000, burnin=500, nthin=2, nchains=1, a1=0.1,
            a2=1, c1=0.1, c2=1, sigma2=10^-3, snpnames, genename){


  RES1=list()
  RESBeta1= RESBeta2= matrix(0,1,mg)
  RESOthers=matrix(0,1,2)
  for(j in 1:nchains){
    RES1[[j]]=CS0(Betah1=Betah1, Betah2=Betah2, Sigmah1=Sigmah1, Sigmah2=Sigmah2, kappa0 = kappa0[j], tau20 = tau20[j],
                  zeta10=zeta10[j], zeta20 = zeta20[j],
                  mg=mg, niter=niter, burnin=burnin, nthin=nthin, a1=a1, a2=a2, c1=c1, c2=c2, sigma2=sigma2, snpnames=snpnames, genename=genename)

    RESBeta1=rbind(RESBeta1,RES1[[j]]$mcmcchain$Beta1)
    RESBeta2=rbind(RESBeta2,RES1[[j]]$mcmcchain$Beta2)
    RESOthers=rbind(RESOthers,cbind(RES1[[j]]$mcmcchain$kappa,RES1[[j]]$mcmcchain$tau2))

  }
  RESBeta1=RESBeta1[-1,]
  RESBeta2=RESBeta2[-1,]
  RESOthers=RESOthers[-1,]

  TabB1=wiqid::simpleRhat(RESBeta1,n.chains=nchains)
  TabB2=wiqid::simpleRhat(RESBeta2,n.chains=nchains)
  TabOthers=wiqid::simpleRhat(RESOthers,n.chains=nchains)

  TabOthers=t(TabOthers)
  colnames(TabOthers)=c("kappa","tau2")

  Tab=cbind(snpnames,TabB1,TabB2)
  colnames(Tab)<-c("Name of SNP", "BGR for Beta_1", "BGR for Beta_2")
  Tabb=list(Tab,TabOthers)


  RES1new=list(RES1,Tabb)
  names(RES1new) <- c("Outputs", "BGR" )

  ifelse(nchains==1,return(RES1),return(RES1new))

}



#' Internal: Continuous Spike
#'
#'
#'
#'
#' @param Betah1 A numerical vector of length mg representing the regression coefficients for the first trait.
#' @param Betah2 A numerical vector of length mg representing the regression coefficients for the second trait.
#' @param Sigmah1 A mg*mg positive definite covariance matrix where it is estimated covariance matrix of Betah1.
#' @param Sigmah2 A mg*mg positive definite covariance matrix where it is estimated covariance matrix of Betah2.
#' @param kappa0 Initial value for kappa.
#' @param tau20 Initial value for tau2.
#' @param zeta10 Initial value for zeta1.
#' @param zeta20 Initial value for zeta2.
#' @param mg Number of variables in the group.
#' @param niter Number of iterations for the Gibbs sampler.
#' @param burnin Number of burn-in iterations.
#' @param nthin The lag of the iterations used for the posterior analysis is defined (or thinning rate).
#' @param a1,a2 Hyperparameters of kappa. Default is a1=0.1 and a2=1.
#' @param c1,c2 Hyperparameters of tau2. Default is c1=0.1 and c2=1.
#' @param sigma2 Variance of spike (multivariate normal distribution with a diagonal covariance matrix with small variance) representing the null effect distribution. Default is 10^-3.
#' @param snpnames Names of variables for the group.
#' @param genename Name of group.
#'
#' @export
CS0=function(Betah1, Betah2, Sigmah1, Sigmah2, kappa0, tau20, zeta10, zeta20, mg, niter=1000, burnin=500, nthin=2, a1=0.1, a2=1, c1=0.1, c2=1, sigma2=10^-3, snpnames, genename){

  Beta1=matrix(0,niter,mg)
  Beta2=matrix(0,niter,mg)
  Beta1[1,]=Betah1
  Beta2[1,]=Betah2
  zeta1=zeta2=rep(0,niter)
  zeta1[1]=zeta10
  zeta2[1]=zeta20

  kappa=rep(0,niter);kappa[1]=kappa0
  tau2=rep(0,niter);tau2[1]=tau20

  K=2
  for(r in 2:niter){
    ##################### Betak ####
    if(zeta1[r-1]==1){
      Sigma1=MASS::ginv(MASS::ginv(Sigmah1)+(tau2[r-1]^(-1))*diag(mg))
      Mean1=Sigma1%*%MASS::ginv(Sigmah1)%*%Betah1
      Beta1[r,]=mvtnorm::rmvnorm(1, mean = Mean1, sigma =Sigma1)
    }

    if(zeta2[r-1]==1){
      Sigma2=MASS::ginv(MASS::ginv(Sigmah2)+(tau2[r-1]^(-1))*diag(mg))
      Mean2=Sigma2%*%MASS::ginv(Sigmah2)%*%Betah2
      Beta2[r,]=mvtnorm::rmvnorm(1, mean = Mean2, sigma =Sigma2)
    }

    if(zeta1[r-1]==0){
      Sigma1=MASS::ginv(MASS::ginv(Sigmah1)+sigma2^-1*diag(mg))
      Mean1=Sigma1%*%MASS::ginv(Sigmah1)%*%Betah1
      Beta1[r,]=mvtnorm::rmvnorm(1, mean = Mean1, sigma =Sigma1)
    }

    if(zeta2[r-1]==0){
      Sigma2=MASS::ginv(MASS::ginv(Sigmah2)+sigma2^-1*diag(mg))
      Mean2=Sigma2%*%MASS::ginv(Sigmah2)%*%Betah2
      Beta2[r,]=mvtnorm::rmvnorm(1, mean = Mean2, sigma =Sigma2)
    }
    ##################### kappa ####
    kappa[r]=rbeta(1,zeta1[r-1]+zeta2[r-1]+a1,K-(zeta1[r-1]+zeta2[r-1])+a2)
    ##################### zetak ####

    w1=1+(kappa[r]/(1-kappa[r]))*exp(log(sigma2/tau2[r-1])*(mg/2)-(.5/sigma2)*t(Beta1[r,])%*%Beta1[r,]*((sigma2/tau2[r-1])-1))
    pzeta1=1-(w1^-1)
    zeta1[r]=rbinom(1,1,pzeta1)

    w2=1+(kappa[r]/(1-kappa[r]))*exp(log(sigma2/tau2[r-1])*(mg/2)-(.5/sigma2)*t(Beta2[r,])%*%Beta2[r,]*((sigma2/tau2[r-1])-1))
    pzeta2=1-(w2^-1)
    zeta2[r]=rbinom(1,1,pzeta2)
    ##################### tau2  ####
    K1=zeta1[r]+zeta2[r]

    if(K1==0)(tau2[r]=invgamma::rinvgamma( 1, shape=c1, rate=c2 ))


    if(K1!=0)(tau2[r]=invgamma::rinvgamma( 1,shape=(K1*mg)/2+c1,
                                           rate=c2+(Beta1[r,]%*%Beta1[r,]+Beta2[r,]%*%Beta2[r,])/2))
  }
  #$$$$$$$$$$$$$$$$$$$$$$$$$ New Postorior samples After Burn-in $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  niter1=burnin+1
  index=niter1:niter
  indexn=index[(index %% nthin) == 0]


  tau2=tau2[indexn]
  kappa=kappa[indexn]
  zeta1=zeta1[indexn]
  zeta2=zeta2[indexn]

  Beta1=Beta1[indexn,]
  Beta2=Beta2[indexn,]
  if(mg==1)((Beta1=as.matrix(Beta1)) & (Beta2=as.matrix(Beta2)))
  PGENE1=PGENE2=rep(0,mg)
  B1CI=B2CI=c()
  for(int in 1:mg){
    B1CI=quantile(Beta1[,int],c(0.025,0.975))
    B2CI=quantile(Beta2[,int],c(0.025,0.975))
    if((0< B1CI[1]) | (0>B1CI[2]))(PGENE1[int]=1)
    if((0<B2CI[1]) | (0>B2CI[2]))(PGENE2[int]=1)
  }

  Geneplotci=PGENE1+PGENE2



  MeanBeta1=apply(Beta1,2,mean)
  MeanBeta2=apply(Beta2,2,mean)

  SDBeta1=apply(Beta1,2,sd)
  SDBeta2=apply(Beta2,2,sd)

  QBeta1=t(apply(Beta1, 2, function(x) quantile(x, c(.025, 0.5, .975))))
  QBeta2=t(apply(Beta2, 2, function(x) quantile(x, c(.025, 0.5, .975))))

  Other=cbind(kappa,tau2)
  MeanOther=apply(Other,2,mean)

  SDOther=apply(Other,2,sd)

  QOther=t(apply(Other, 2, function(x) quantile(x, c(.025, 0.5, .975))))

  #$$$$$$$$$$$$$$$$$$$$$$$$ locFDR $$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  pzeta1=pzeta2=pz=theta=rep(0,length(indexn))
  for(r in 1:length(indexn)){

    w1=1+(kappa[r]/(1-kappa[r]))*exp(log(sigma2/tau2[r])*(mg/2)-(.5/sigma2)*t(Beta1[r,])%*%Beta1[r,]*((sigma2/tau2[r])-1))
    pzeta1[r]=(w1^-1)


    w2=1+(kappa[r]/(1-kappa[r]))*exp(log(sigma2/tau2[r])*(mg/2)-(.5/sigma2)*t(Beta2[r,])%*%Beta2[r,]*((sigma2/tau2[r])-1))
    pzeta2[r]=(w2^-1)

    pz[r]=pzeta1[r]*pzeta2[r]
    theta[r]=(1-pzeta1[r])*(1-pzeta2[r])

  }
  PSI1=mean(zeta1)
  PSI2=mean(zeta2)

  locFDR=mean(pz[is.nan(pz)==FALSE])
  locFDR
  ResFDR=locFDR
  ################################################################
  BF=((1-locFDR)/(locFDR))*(a2^2/((a1+a2)^2-a2^2))
  BF

  ResBF=log10(BF)
  Pleop=mean(theta[is.nan(theta)==FALSE])
  ResPleop=Pleop



  mcmcchain=list(kappa=kappa,tau2=tau2,Beta1=Beta1,Beta2=Beta2,zeta1=zeta1,zeta2=zeta2)


  Reslast= list("Name of Gene"=genename,"Name of SNP"=snpnames,log10BF=ResBF,lBFDR=ResFDR,theta=ResPleop,
                "# studies nonzero signal by CI"=Geneplotci,
                PPA1=PSI1,PPA2=PSI2)


  Trait_1= cbind(snpnames, MeanBeta1, SDBeta1, QBeta1)
  colnames(Trait_1)= cbind("Name of SNP", "Mean", "SD", "val2.5pc",	"Median",	"val97.5pc")

  Trait_2= cbind(snpnames, MeanBeta2, SDBeta2, QBeta2)
  colnames(Trait_2)= cbind("Name of SNP", "Mean", "SD", "val2.5pc",	"Median",	"val97.5pc")



  Others= cbind(MeanOther, SDOther, QOther)
  colnames(Others)= cbind( "Mean", "SD", "val2.5pc",	"Median",	"val97.5pc")

  OUT=list(sim.matrix=mcmcchain, Criteria=Reslast, Trait_1, Trait_2, Others)

  names(OUT) <- c("mcmcchain", "Criteria", "Statistics of Trait 1 for Beta_1", "Statistics of Trait 2  for Beta_2", "Other Parameters")

  ;OUT
}



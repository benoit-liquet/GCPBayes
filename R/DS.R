#' Dirac Spike
#'
#'
#' @description
#' Run a Gibbs sampler for a multivariate Bayesian sparse group selection model with Dirac spike prior for detecting pleiotropic effects on two traits. This function is designed for summary statistics containing estimated regression coefficients and their estimated covariance matrices.
#'
#'
#' @details
#' Run a Gibbs sampler using a Dirac spike
#'
#'
#' @param Betah1 A numerical vector of length mg representing the regression coefficients for the first trait.
#' @param Betah2 A numerical vector of length mg representing the regression coefficients for the second trait.
#' @param Sigmah1 A mg*mg positive definite covariance matrix where it is estimated covariance matrix Betah1.
#' @param Sigmah2 A mg*mg positive definite covariance matrix where it is estimated covariance matrix Betah2.
#' @param kappa0 Initial value for kappa (its dimension is equal to nchains).
#' @param sigma20 Initial value for sigma2 (its dimension is equal to nchains).
#' @param mg Number of variables in the group.
#' @param niter Number of iterations for the Gibbs sampler.
#' @param burnin Number of burn-in iterations.
#' @param nthin The lag of the iterations used for the posterior analysis is defined (or thinning rate).
#' @param nchains Number of Markov chains, when nchains>1, the function calculates the Gelman-Rubin convergence statistic, as modified by Brooks and Gelman (1998).
#' @param a1,a2 Hyperparameters of kappa. Default is a1=0.1 and a2=0.1.
#' @param d1,d2 Hyperparameters of sigma2. Default is d1=0.1 and d2=0.1.
#' @param snpnames Names of variables for the group.
#' @param genename Name of group.
#'
#'
#' @return
#' - mcmcchain: The list of simulation output for all parameters.
#' - Criteria: genename, snpnames, log10BF, lBFDR, theta, PPA1, PPA2 and the number of studies with nonzero signal by credible interval (CI).
#' - Statistics of Trait 1 for Beta_1: Summary statistics including Mean, SD, val2.5pc, Median,	val97.5pc and BGR (if nchains>1).
#' - Statistics of Trait 2  for Beta_2: Summary statistics including Mean, SD, val2.5pc, Median,	val97.5pc and BGR (if nchains>1).
#' - Other Parameters: Summary statistics for kappa and sigma2 including Mean, SD, val2.5pc, Median,	val97.5pc and BGR (if nchains>1).
#'
#'
#' @author Taban Baghfalaki.
#'
#' @references
#' T. Baghfalaki, P.E. Sugier, T. Truong, A.T. Pettitt, K. Mengersen and B. Liquet. (2020). Bayesian meta-analysis models to Cross Cancer Genomic Investigation of pleiotropic effects using group structure. *Submitted in Statistics in Medicine*.
#'
#' @example inst/exampleDS.R
#'
#' @md
#' @export





DS=function(Betah1, Betah2, Sigmah1, Sigmah2, kappa0, sigma20, mg, niter=1000, burnin=500, nthin=2, nchains=1, a1=0.1, a2=0.1, d1=0.1, d2=0.1, snpnames, genename){

  RES1=list()
  RESBeta1= RESBeta2= matrix(0,1,mg)
  RESOthers=matrix(0,1,2)
  for(j in 1:nchains){
    RES1[[j]]=DS0(Betah1=Betah1, Betah2=Betah2, Sigmah1=Sigmah1, Sigmah2=Sigmah2, kappa0 = kappa0[j], sigma20 = sigma20[j],
                  mg=mg, niter=niter, burnin=burnin, nthin=nthin, a1=a1, a2=a2, d1=d1, d2=d2, snpnames=snpnames, genename=genename)

    RESBeta1=rbind(RESBeta1,RES1[[j]]$mcmcchain$Beta1)
    RESBeta2=rbind(RESBeta2,RES1[[j]]$mcmcchain$Beta2)
    RESOthers=rbind(RESOthers,cbind(RES1[[j]]$mcmcchain$kappa,RES1[[j]]$mcmcchain$sigma2))

  }
  RESBeta1=RESBeta1[-1,]
  RESBeta2=RESBeta2[-1,]
  RESOthers=RESOthers[-1,]

  TabB1=wiqid::simpleRhat(RESBeta1,n.chains=nchains)
  TabB2=wiqid::simpleRhat(RESBeta2,n.chains=nchains)
  TabOthers=wiqid::simpleRhat(RESOthers,n.chains=nchains)

  TabOthers=t(TabOthers)
  colnames(TabOthers)=c("kappa","sigma2")

  Tab=cbind(snpnames,TabB1,TabB2)
  colnames(Tab)<-c("Name of SNP", "BGR for Beta_1", "BGR for Beta_2")
  Tabb=list(Tab,TabOthers)


  RES1new=list(RES1,Tabb)
  names(RES1new) <- c("Outputs", "BGR" )

  ifelse(nchains==1,return(RES1),return(RES1new))

}





DS0=function(Betah1, Betah2, Sigmah1, Sigmah2, kappa0, sigma20, mg, niter=1000, burnin, nthin=2, a1=0.1, a2=0.1, d1=0.1, d2=0.1, snpnames, genename){

  Beta1=matrix(0,niter,mg)
  Beta2=matrix(0,niter,mg)

  Beta1[1,]=Betah1
  Beta2[1,]=Betah2
  kappa=rep(0,niter);kappa[1]=kappa0
  sigma2=rep(0,niter);sigma2[1]=sigma20

  psi1=psi2=rep(1,niter)

  MeanBeta1=MeanBeta2=MedianBeta1=MedianBeta2=c()
  Geneplotci=Geneplotmed=c()
  PROB1j=PROB2j=c()
  PROB1=PROB2=c()
  K=2
  for(r in 2:niter){
    ##################### Betak ####
    Omega1=MASS::ginv(Sigmah1)+(1/sigma2[r-1])*diag(mg)
    Mean1=MASS::ginv(Omega1)%*%MASS::ginv(Sigmah1)%*%Betah1
    kappat1=(1/(1+(kappa[r-1]/(1-kappa[r-1])* exp(.5*t(Mean1)%*%Omega1%*%Mean1-mg/2*log(sigma2[r-1])-.5*determinant(Omega1, logarithm = TRUE)$modulus))))
    Tab=rbinom(1,1,as.numeric(kappat1))
    HH1=MASS::ginv(Omega1)
    gdata::lowerTriangle(HH1) <- gdata::upperTriangle(HH1, byrow=TRUE)
    if(Tab==0)(Beta1[r,]=mvtnorm::rmvnorm(1, mean = Mean1, sigma =HH1))
    if(Tab==1)(Beta1[r,]=rep(0,mg))

    Omega2=MASS::ginv(Sigmah2)+(1/sigma2[r-1])*diag(mg)
    Mean2=MASS::ginv(Omega2)%*%MASS::ginv(Sigmah2)%*%Betah2
    kappat2=(1/(1+(kappa[r-1]/(1-kappa[r-1])* exp(.5*t(Mean2)%*%Omega2%*%Mean2-mg/2*log(sigma2[r-1])-.5*determinant(Omega2, logarithm = TRUE)$modulus))))
    Tab=rbinom(1,1,as.numeric(kappat2))
    HH2=MASS::ginv(Omega2)
    gdata::lowerTriangle(HH2) <- gdata::upperTriangle(HH2, byrow=TRUE)
    if(Tab==0)(Beta2[r,]=mvtnorm::rmvnorm(1, mean = Mean2, sigma =HH2))
    if(Tab==1)(Beta2[r,]=rep(0,mg))
    ##################### kappa ####
    K0=0
    if(max(Beta1[r,])==0)((K0=K0+1) & (psi1[r]=0))
    if(max(Beta2[r,])==0)((K0=K0+1) & (psi2[r]=0))
    kappa[r]=rbeta(1,K-K0+a1,K0+a2)

    ##################### sigma2 ####

    sigma2[r]=invgamma::rinvgamma(1,mg*(K-K0)/2+d1,(Beta1[r,]%*%Beta1[r,]+Beta2[r,]%*%Beta2[r,])/2+d2)
  }

  #$$$$$$$$$$$$$$$$$$$$$$$$$ New Postorior samples After Burn-in $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  niter1=burnin+1
  index=niter1:niter
  indexn=index[(index %% nthin) == 0]

  kappa=kappa[indexn]
  sigma2=sigma2[indexn]
  Beta1=Beta1[indexn,]
  Beta2=Beta2[indexn,]
  if(mg==1)((Beta1=as.matrix(Beta1)) & (Beta2=as.matrix(Beta2)))
  psi1=psi1[indexn]
  psi2=psi2[indexn]
  PPA1=mean(psi1)
  PPA2=mean(psi2)

  PGENE1=PGENE2=rep(0,mg)
  B1CI=B2CI=c()
  for(int in 1:mg){
    B1CI=quantile(Beta1[,int],c(0.025,0.975))
    B2CI=quantile(Beta2[,int],c(0.025,0.975))
    if((0< B1CI[1]) | (0>B1CI[2]))(PGENE1[int]=1)
    if((0<B2CI[1]) | (0>B2CI[2]))(PGENE2[int]=1)
  }

  Geneplotci=PGENE1+PGENE2

  PGENE1m=PGENE2m=rep(0,mg)
  B1CI=B2CI=c()
  for(int in 1:mg){
    B1CI=median(Beta1[,int])
    B2CI=median(Beta2[,int])
    if(B1CI!=0)(PGENE1m[int]=1)
    if(B2CI!=0)(PGENE2m[int]=1)
  }

  Geneplotmed=PGENE1m+PGENE2m


  MeanBeta1=apply(Beta1,2,mean)
  MeanBeta2=apply(Beta2,2,mean)


  QBeta1=t(apply(Beta1, 2, function(x) quantile(x, c(.025, 0.5, .975))))
  QBeta2=t(apply(Beta2, 2, function(x) quantile(x, c(.025, 0.5, .975))))


  SDBeta1=apply(Beta1,2,sd)
  SDBeta2=apply(Beta2,2,sd)


  Other=cbind(kappa,sigma2)
  MeanOther=apply(Other,2,mean)

  SDOther=apply(Other,2,sd)

  QOther=t(apply(Other, 2, function(x) quantile(x, c(.025, 0.5, .975))))

  #$$$$$$$$$$$$$$$$$$$$$$$$ locFDR, theta, BF $$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  pz1=pz2=pz=theta=rep(0,length(indexn))
  for(r in 1:length(indexn)){
    Omega1=MASS::ginv(Sigmah1)+(1/sigma2[r])*diag(mg)
    Mean1=MASS::ginv(Omega1)%*%MASS::ginv(Sigmah1)%*%Betah1
    pz1[r]=(1/(1+(kappa[r]/(1-kappa[r])* exp(.5*t(Mean1)%*%Omega1%*%Mean1-mg/2*log(sigma2[r])-.5*determinant(Omega1, logarithm = TRUE)$modulus))))


    Omega2=MASS::ginv(Sigmah2)+(1/sigma2[r])*diag(mg)
    Mean2=MASS::ginv(Omega2)%*%MASS::ginv(Sigmah2)%*%Betah2

    pz2[r]=(1/(1+(kappa[r]/(1-kappa[r])* exp(.5*t(Mean2)%*%Omega2%*%Mean2-mg/2*log(sigma2[r])-.5*determinant(Omega2, logarithm = TRUE)$modulus))))

    pz[r]=pz1[r]*pz2[r]
    theta[r]=(1-pz1[r])*(1-pz2[r])

  }


  locFDR=mean(pz[is.nan(pz)==FALSE])
  locFDR
  ResFDR=locFDR

  BF=((1-locFDR)/(locFDR))*(a2^2/((a1+a2)^2-a2^2))
  BF

  Pleop=mean(theta[is.nan(theta)==FALSE])
  ResPleop=Pleop

  ResBF=log10(BF)


  mcmcchain=list(kappa=kappa,sigma2=sigma2,Beta1=Beta1,Beta2=Beta2)


  Reslast= list("Name of Gene"=genename,"Name of SNP"=snpnames,log10BF=ResBF,lBFDR=ResFDR,theta=ResPleop,
                "# studies nonzero signal by CI"=Geneplotci,#"# studies nonzero signal by Med"=Geneplotmed,
                PPA1=PPA1,PPA2=PPA2)


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


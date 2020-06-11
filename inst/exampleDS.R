mg=10
SD=diag(0.05,mg)
corr=0.5
R=matrix(corr,mg,mg)+(1-corr)*diag(mg)
Sigma=crossprod(crossprod(SD,R),SD)
sign=rbinom(mg,1,.5)
sign[sign==0]=-1
betat=rep(1,mg)
betat=betat*sign
Betah=mvtnorm::rmvnorm(2,betat,Sigma)

snpnames=1:mg
genename="simulated_data"


Betah1=Betah[1,]; Betah2=Betah[2,];
Sigmah1=Sigma; Sigmah2=Sigma;



res1 = DS(Betah1, Betah2,
          Sigmah1, Sigmah2,
          kappa0=0.5, sigma20=1,
          mg=mg, niter=1000, burnin=500,
          nchains=1, nthin=2, a1=0.1, a2=0.1, d1=0.1, d2=0.1,
          snpnames, genename)

res2 = DS(Betah1, Betah2,
          Sigmah1, Sigmah2,
          kappa0=c(0.5,0.25,0.6), sigma20=c(1,1.2,1.5),
          mg=mg, niter=1000, burnin=500,
          nchains=3, nthin=2, a1=0.1, a2=0.1, d1=0.1, d2=0.1,
          snpnames, genename)

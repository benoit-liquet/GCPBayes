mg=10
SD=diag(0.05,mg)
corr=0.5
R=matrix(corr,mg,mg)+(1-corr)*diag(mg)
Sigma=crossprod(crossprod(SD,R),SD)
sign=rbinom(mg,1,.5)
sign[sign==0]=-1
Betat=rep(1,mg)
Betat=Betat*sign
Betah=mvtnorm::rmvnorm(2,Betat,Sigma)

snpnames=1:mg
genename="simulated_data"


#######################
Betah1=Betah[1,]; Betah2=Betah[2,];
Sigmah1=Sigma; Sigmah2=Sigma;
a1=0.1; a2=0.1; d1=0.1; d2=0.1; c1=1; c2=1; e2=1; N00=100



res1 = HS(Betah1, Betah2,
          Sigmah1, Sigmah2, kappa0 = 0.5, kappastar0=0.5, sigma20 =1, s20 =1,
          mg=mg, niter=1000, burnin=500, nthin=3, nchains=1, a1=a1, a2=a2, d1=d1,
          d2=d2, c1=c1, c2=c2, e2=e2, snpnames, genename)



res2 = HS(Betah1, Betah2,
          Sigmah1, Sigmah2, kappa0 = c(0.5,0.25,0.7), kappastar0=c(0.5,0.25,0.2),
          sigma20 =c(1,1.5,2), s20 =c(1,1.5,2),
          mg=mg, niter=1000, burnin=500, nthin=3, nchains=3, a1=a1, a2=a2, d1=d1,
          d2=d2, c1=c1, c2=c2, e2=e2, snpnames, genename)



plot_trace_HS(res1, kappa=TRUE, kappastar=TRUE,  sigma2=TRUE,
              Beta1=1:mg, Beta2=1:mg,
              nrow1 = 2, nrow2 = 2)

plot_trace_HS(res2, kappa=TRUE, kappastar=TRUE,  sigma2=TRUE,
              Beta1=1:mg, Beta2=1:mg,
              nrow1 = 2, nrow2 = 2)




























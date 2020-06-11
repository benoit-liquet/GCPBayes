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
kappa0=0.5; kappastar0=0.5; sigma20=1; mg=mg;
a1=0.1; a2=0.1; d1=0.1; d2=0.1; c1=1; c2=1; s20=1; N00=100
pvalue=2*pnorm(-abs(Betah/sqrt(diag(Sigma))))

zinit=rep(0,2)
for(j in 1:2){
  index=1:mg
  PVALUE=p.adjust(pvalue[j,])
  SIGNALS=index[PVALUE<0.05]
  modelf1=rep(0,mg)
  modelf1[SIGNALS]=1
  if(max(modelf1)==1)(zinit[j]=1)
}


res1 = CS(Betah1, Betah2,
          Sigmah1, Sigmah2,
          kappa0=0.5, tau20=1,
          zeta10=zinit[1], zeta20 =zinit[2],
          mg=mg, niter=1000, burnin=500,
          nchains=1, nthin=2, a1=.1, a2=1, c1=0.1, c2=1, sigma2=10^-3,
          snpnames, genename)

res2 = CS(Betah1, Betah2,
          Sigmah1, Sigmah2,
          kappa0=c(0.5,0.25,0.5), tau20=c(1,1.25,1.5),
          zeta10=rep(zinit[1],3), zeta20 =rep(zinit[2],3),
          mg=mg, niter=1000, burnin=500,
          nchains=3, nthin=5, a1=.1, a2=1, c1=0.1, c2=1,
          sigma2=10^-3, snpnames, genename)

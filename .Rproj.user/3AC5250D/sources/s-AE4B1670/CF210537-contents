SNVtest<-
function(tumor1, tumor2, pfreq, nrep=1000) {
  w<-tumor1+tumor2>0
  tumor1<-tumor1[w]
  tumor2<-tumor2[w]
  pfreq<-pfreq[w]
  
  n <- length(pfreq)
  if (n!=length(tumor1) |n!=length(tumor2)) stop("Number of mutations in tumor1 and tumor2 should be equal to the length of pfreq")
  if (any(!(tumor1 %in% c(0,1))) | any(!(tumor2 %in% c(0,1)))) stop("tumor1 and tumor2 should be a vectors with only values 0 and 1")
  if (any(pfreq>1 | pfreq<0)) stop("pfreq should take numerical values between 0 and 1")
  
  # test is based on loci whether either tumor is having a mutation 
  ii <- which(tumor1 + tumor2 > 0)
  if (length(ii) == 0) return(rep(NA,13)) 
  cn <- length(ii)
  
  cp <- pfreq[ii]/(2 - pfreq[ii])
  gen<-matrix(runif(nrep*cn)<=rep(cp,nrep),ncol=nrep)
  

  Lik<-
  function(ksi,x,y,p)
  {a<-which(x+y==2)
  b<-which(x+y==1)
  e<-c(a,b)
  #sum(log( (ksi*p+(1-ksi)*p^2)[a]))+sum((log(2*(1-ksi)*p*(1-p))[b]))-sum(log(ksi+(1-ksi)*(1-p[e])^2))
  sum(log( (ksi*p+(1-ksi)*p^2)[a]))+sum(log(2*(1-ksi)*(p*(1-p))[b]))-
    sum(log( (ksi*p+(1-ksi)*p^2+2*(1-ksi)*p*(1-p))[e]))
  }
  Lik2<-
  function(ksi,x,p)
  {a<-which(x==1)#matching mutations indicator
  b<-which(x==0)
  e<-c(a,b)
  #sum(log( (ksi*p+(1-ksi)*p^2)[a]))+sum((log(2*(1-ksi)*p*(1-p))[b]))-sum(log(ksi+(1-ksi)*(1-p[e])^2))
  sum(log( (ksi*p+(1-ksi)*p^2)[a]))+sum(log(2*(1-ksi)*(p*(1-p))[b]))-
    sum(log( (ksi*p+(1-ksi)*p^2+2*(1-ksi)*p*(1-p))[e]))
  }
  
  fLR<-  function(x,y,p)
  {o<-optim(0.5,function(ksi){Lik(ksi,x,y,p)},control=list(fnscale =-1),method="L-BFGS-B",lower=0.00001,upper=0.999999)
  c(o$value-Lik(0,x,y,p),o$par)
  }  
  fLR2<-
  function(x,p)
  {o<-optim(0.5,function(ksi){Lik2(ksi,x,p)},control=list(fnscale =-1),method="L-BFGS-B",lower=0.00001,upper=0.999999)
  o$value-Lik2(0,x,p)
  }
  
  #cond LR test + ref
  lr<-fLR(tumor1, tumor2, pfreq)
  f2<-function(x){fLR2(as.numeric(x),pfreq[tumor1 + tumor2 > 0])}
  ts2<-apply (gen,2,f2) 
    out <- c(sum(tumor1),sum(tumor2),sum(tumor1==tumor2 & tumor2==1) ,
           lr,pmax(sum(ts2>=lr[1]) /nrep, 1/nrep ))

  names(out) <- c("n1","n2","n_match", "LRstat","maxKsi","LRpvalue")
  
  
  out
}
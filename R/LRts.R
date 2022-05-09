LRts <-
function(tum,p, MX)
{a<-tum[1,]
e<-tum[2,]
f<-tum[3,]
g<-tum[4,]
h<-tum[5,]
fpg<-f+g
m<-0.5
psq<-p^2
p1p<-p*(1-p)
p12<-(1-p)^2
ema<-e-a

lcpmm<-function(y){
x<-y[1]
m<-y[2]
sum(log(x*p+(m^2+(1-m)^2)*(1-x)^2*psq/(1-x*p))*a + 
log( (2*m*(1-m))*(1-x)^2*psq/(1-x*p))*(ema)  +
log (2*(1-x)*p1p/(1-x*p))*(fpg) + log(p12/(1-x*p))*h)
}

 lcpm<-function(x){
 prod((x*p+(m^2+(1-m)^2)*(1-x)^2*psq/(1-x*p))^a * 
 ( (2*m*(1-m))*(1-x)^2*p^2/(1-x*p))^(ema)  * 
 (2*(1-x)*p1p/(1-x*p))^(fpg) * (p12/(1-x*p))^h)
}

 lcpmm2<-function(y){
x<-0
m<-y
sum( log(x*p+(m^2+(1-m)^2)*(1-x)^2*psq/(1-x*p))*a + 
log( (2*m*(1-m))*(1-x)^2*p^2/(1-x*p))*(ema)  +
log (2*(1-x)*p1p/(1-x*p))*(fpg) +log (p12/(1-x*p))*h)
}

lcpm2<-function(y){
x<-y[1]
m<-y[2]
prod( (x*p+(m^2+(1-m)^2)*(1-x)^2*psq/(1-x*p))^a * 
( (2*m*(1-m))*(1-x)^2*p^2/(1-x*p))^(ema)  * 
(2*(1-x)*p1p/(1-x*p))^(fpg) * (p12/(1-x*p))^h)
}



if (sum(a)==0 & sum(a-e)==0)
  {
  m<-0.5;
  opt<-optimize(lcpm,interval=c(0,1),maximum=TRUE)
  ts2<-2*log(opt$objective/lcpm(0))
   return(c(max(0,ts2),0.5,opt$maximum,0.5))
 }  
else {
mles<-NA
try( mles<- optim(c(0.5,0.75),lcpmm,lower=c(0,0.5),upper=c(0.99999999,MX),
method="L-BFGS-B",control=list(fnscale=-1)) ,silent=TRUE)
if (is.null(mles$par) ) {print(tum); print("numerical optimization didn't converge")} 
else
{

m0hat<-optimize(lcpmm2,interval=c(0.5,MX),maximum=TRUE)$maximum
ts2<-2*(mles$value-log(lcpm2(c(0,m0hat))))
  }


 return(c(max(0,ts2),m0hat,mles$par))
}
}


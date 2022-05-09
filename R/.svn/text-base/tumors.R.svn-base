tumors <-
function(p,cc,m)
{
J<-length(p)
if (cc==0) t0<-rep(0,J) else {t0<-runif(J)
t0[t0<=p*cc]<-1
t0[t0<1]<-0
t0[t0==1 & runif(J)<=m]<-2
}
t1<-t0
t2<-t0
x<-runif(J)
t2[t2==0 & x<=(1-cc)*p*m/(1-cc*p)]<-2
t2[t2==0 & x>(1-cc)*p*m/(1-cc*p)& x<=(1-cc)*p/(1-cc*p)]<-1
x<-runif(J)
t1[t1==0 & x<=(1-cc)*p*m/(1-cc*p)]<-2
t1[t1==0 & x>(1-cc)*p*m/(1-cc*p)& x<=(1-cc)*p/(1-cc*p)]<-1

a<-rep(0,J)
e<-rep(0,J)
f<-rep(0,J)
g<-rep(0,J)
h<-rep(0,J)

e[t1>0 & t2>0]<-1
a[t1==t2 & t1>0 & t2>0]<-1
f[t1==0 & t2>0]<-1
g[t1>0 & t2==0]<-1
h[t1==0 & t2==0]<-1

return(rbind(a,e,f,g,h))
}

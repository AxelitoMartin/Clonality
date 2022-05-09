LRtesting3or4tumors <-
function(LOHtable,ptlist,refLOHtable=NULL, pfreq=NULL,noloh,loh1,loh2,Nsim=100,m=0.5)
{ 



 ############### generating tumors from top 2
tumors4top2<-function(p,cc=0,c2=0,c3=0,m=0.5,rettype="matrix")
{
J<-length(p)
if (cc==0) t0<-rep(0,J) 
else 
{
t0<-runif(J)
t0[t0<=p*cc]<-1
t0[t0<1]<-0
t0[t0==1 & runif(J)<=m]<-2
}
t1<-t0
t2<-t0
t3<-t0 
t4<-t0

x<-runif(J)
t1[t1==0 & x<=(1-cc)*p*m/(1-cc*p)]<-2
t1[t1==0 & x>(1-cc)*p*m/(1-cc*p)& x<=(1-cc)*p/(1-cc*p)]<-1

if (c2==0 & c3==0)
{

x<-runif(J)
t3[t3==0 & x<=(1-cc)*p*m/(1-cc*p)]<-2
t3[t3==0 & x>(1-cc)*p*m/(1-cc*p)& x<=(1-cc)*p/(1-cc*p)]<-1
x<-runif(J)
t2[t2==0 & x<=(1-cc)*p*m/(1-cc*p)]<-2
t2[t2==0 & x>(1-cc)*p*m/(1-cc*p)& x<=(1-cc)*p/(1-cc*p)]<-1
x<-runif(J)
t4[t4==0 & x<=(1-cc)*p*m/(1-cc*p)]<-2
t4[t4==0 & x>(1-cc)*p*m/(1-cc*p)& x<=(1-cc)*p/(1-cc*p)]<-1
}


else 
if (c2!=0 & c3==0)

{

x<-runif(J)
t2[t2==0 & x<=c2*p*m/(1-cc*p)]<-2
t2[t2==0 & x>c2*p*m/(1-cc*p)& x<=c2*p/(1-cc*p)]<-1
t4<-t3<-t2

x<-runif(J)
t3[t3==0 & x<=(1-cc-c2)*p*m/((1-(c2+cc)*p))]<-2
t3[t3==0 & x>(1-cc-c2)*p*m/((1-(c2+cc)*p))& x<=(1-cc-c2)*p/((1-(c2+cc)*p))]<-1
x<-runif(J)
t2[t2==0 & x<=(1-cc-c2)*p*m/((1-(c2+cc)*p))]<-2
t2[t2==0 & x>(1-cc-c2)*p*m/((1-(c2+cc)*p))& x<=(1-cc-c2)*p/((1-(c2+cc)*p))]<-1
x<-runif(J)
t4[t4==0 & x<=(1-cc-c2)*p*m/((1-(c2+cc)*p))]<-2
t4[t4==0 & x>(1-cc-c2)*p*m/((1-(c2+cc)*p))& x<=(1-cc-c2)*p/((1-(c2+cc)*p))]<-1


}

else 
if (c2==0 & c3!=0)

{
x<-runif(J)
t2[t2==0 & x<=(1-cc)*p*m/((1-(cc)*p))]<-2
t2[t2==0 & x>(1-cc)*p*m/((1-(cc)*p))& x<=(1-cc)*p/((1-(cc)*p))]<-1

x<-runif(J)
t3[t3==0 & x<=c3*p*m/(1-cc*p)]<-2
t3[t3==0 & x>c3*p*m/(1-cc*p)& x<=c3*p/(1-cc*p)]<-1
t4<-t3

x<-runif(J)
t3[t3==0 & x<=(1-cc-c2-c3)*p*m/((1-(c3+c2+cc)*p))]<-2
t3[t3==0 & x>(1-cc-c2-c3)*p*m/((1-(c3+c2+cc)*p))& x<=(1-cc-c2-c3)*p/((1-(c3+c2+cc)*p))]<-1

x<-runif(J)
t4[t4==0 & x<=(1-cc-c2-c3)*p*m/((1-(c3+c2+cc)*p))]<-2
t4[t4==0 & x>(1-cc-c2-c3)*p*m/((1-(c3+c2+cc)*p))& x<=(1-cc-c2-c3)*p/((1-(c3+c2+cc)*p))]<-1


}



else 
if (c2!=0 & c3!=0)

{

x<-runif(J)
t2[t2==0 & x<=c2*p*m/(1-cc*p)]<-2
t2[t2==0 & x>c2*p*m/(1-cc*p)& x<=c2*p/(1-cc*p)]<-1
t4<-t3<-t2

x<-runif(J)
t2[t2==0 & x<=(1-cc-c2)*p*m/((1-(c2+cc)*p))]<-2
t2[t2==0 & x>(1-cc-c2)*p*m/((1-(c2+cc)*p))& x<=(1-cc-c2)*p/((1-(c2+cc)*p))]<-1

x<-runif(J)
t3[t3==0 & x<=c3*p*m/(1-(cc+c2)*p)]<-2
t3[t3==0 & x>c3*p*m/(1-(cc+c2)*p)& x<=c3*p/(1-(cc+c2)*p)]<-1
t4<-t3


x<-runif(J)
t3[t3==0 & x<=(1-cc-c2-c3)*p*m/((1-(c3+c2+cc)*p))]<-2
t3[t3==0 & x>(1-cc-c2-c3)*p*m/((1-(c3+c2+cc)*p))& x<=(1-cc-c2-c3)*p/((1-(c3+c2+cc)*p))]<-1

x<-runif(J)
t4[t4==0 & x<=(1-cc-c2-c3)*p*m/((1-(c3+c2+cc)*p))]<-2
t4[t4==0 & x>(1-cc-c2-c3)*p*m/((1-(c3+c2+cc)*p))& x<=(1-cc-c2-c3)*p/((1-(c3+c2+cc)*p))]<-1


}


if (rettype=="matrix") return(cbind(t1,t2,t3,t4))
else return(paste(t1,t2,t3,t4,sep=""))
}



########## probabilities  of outcomes in topology 2

prob4.full.top2<-function(pi,cc,c2,c3,m)
{cbind(

#1111   mean(t %in% c("1111"))*2 
cc*pi+(1-cc)*pi*(c2)*pi * (m^2+(1-m)^2)/ ((1-cc*pi) )+
(1-cc)*pi*(1-cc-c2)*pi*c3*pi * (m^3+(1-m)^3)/ ((1-cc*pi)*(1-(cc+c2)*pi) )+
(1-cc)*pi*(1-cc-c2)*pi*(1-cc-c2-c3)^2*pi^2* (m^4+(1-m)^4)/ ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c2+c3)*pi))
,
#######


#1112,1121  mean(t %in% c("1112","1121"))*2 
(1-cc)*pi*(1-cc-c2)*pi *(1-cc-c2-c3)^2*pi^2 * 2*(m^3*(1-m)+m*(1-m)^3)/ ((1-cc*pi)*(1-(cc+c2+c3)*pi) *(1-(cc+c2)*pi))
,

#1211 mean(t %in% c("1211")) *2
(1-cc)*pi*(1-cc-c2)*pi *c3*pi * (m^2*(1-m)+m*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2)*pi))+
(1-cc)*pi*(1-cc-c2)*pi *(1-cc-c2-c3)^2*pi^2 * (m^3*(1-m)+m*(1-m)^3)/ ((1-cc*pi)*(1-(cc+c2+c3)*pi) *(1-(cc+c2)*pi))
,

#2111  mean(t %in% c("2111")) *2
(1-cc)*pi*c2*pi * 2*(m*(1-m))/ ((1-cc*pi))+
(1-cc)*pi*(1-cc-c2)*pi *c3*pi * (m^2*(1-m)+m*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2)*pi))+
(1-cc)*pi*(1-cc-c2)*pi *(1-cc-c2-c3)^2*pi^2 * (m^3*(1-m)+m*(1-m)^3)/ ((1-cc*pi)*(1-(cc+c2+c3)*pi) *(1-(cc+c2)*pi))
,

#######
#1122  mean(t %in% c("1122")) *2
(1-cc)*pi*(1-cc-c2)*pi *c3*pi * (m^2*(1-m)+m*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2)*pi))+
(1-cc)*pi*(1-cc-c2)*pi *(1-cc-c2-c3)^2*pi^2 * 2*(m^2*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2+c3)*pi) *(1-(cc+c2)*pi))
,


#1212,1221   mean(t %in% c("1212","1221")) *2
(1-cc)*pi*(1-cc-c2)*pi *(1-cc-c2-c3)^2*pi^2 * 4*(m^2*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2+c3)*pi) *(1-(cc+c2)*pi))
,



#######
#1110,1101   mean(t %in% c("1110","1101"))*2
(1-pi)*(1-cc)*pi*(1-cc-c2)*pi *(1-cc-c2-c3)*pi * 2*(m^3+(1-m)^3)/ ((1-cc*pi)*(1-(cc+c2+c3)*pi) *(1-(cc+c2)*pi))
,

#0111,  mean(t %in% c("0111")) *2
(1-pi)*c2*pi / ((1-cc*pi) )+
(1-pi)*(1-cc-c2)*pi*c3*pi * (m^2+(1-m)^2)/ ((1-cc*pi) *(1-(cc+c2)*pi))+
(1-pi)*(1-cc-c2)*pi *(1-cc-c2-c3)^2*pi^2 * (m^3+(1-m)^3)/ ((1-cc*pi)*(1-(cc+c2+c3)*pi) *(1-(cc+c2)*pi))
,

#1011   mean(t %in% c("1011")) *2
(1-pi)*(1-cc)*pi*c3*pi * (m^2+(1-m)^2)/ ((1-cc*pi) *(1-(cc+c2)*pi))+
(1-pi)*(1-cc)*pi *(1-cc-c2-c3)^2*pi^2 * (m^3+(1-m)^3)/ ((1-cc*pi)*(1-(cc+c2+c3)*pi) *(1-(cc+c2)*pi))
,


#######

#1012,1021 mean(t %in% c("1012","1021")) *2
(1-pi)*(1-cc)*pi *(1-cc-c2-c3)^2*pi^2 *2* (m^2*(1-m)+m*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2+c3)*pi) *(1-(cc+c2)*pi))
,

#0112,0121 mean(t %in% c("0112","0121")) *2
(1-pi)*(1-cc-c2)*pi *(1-cc-c2-c3)^2*pi^2 *2* (m^2*(1-m)+m*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2+c3)*pi) *(1-(cc+c2)*pi))
,


#1201,1210,2101,2110 mean(t %in% c("1201","1210","2101","2110"))  *2
(1-pi)*(1-cc)*pi*(1-cc-c2)*pi *(1-cc-c2-c3)*pi *4* (m^2*(1-m)+m*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2+c3)*pi) *(1-(cc+c2)*pi))
,

#0122   mean(t %in% c("0122")) *2
(1-pi)*(1-cc-c2)*pi*c3*pi * 2*(m*(1-m))/ ((1-cc*pi) *(1-(cc+c2)*pi))+
(1-pi)*(1-cc-c2)*pi *(1-cc-c2-c3)^2*pi^2 * (m^2*(1-m)+m*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2+c3)*pi) *(1-(cc+c2)*pi))
,

#1022 mean(t %in% c("1022")) *2
(1-pi)*(1-cc)*pi*c3*pi * 2*(m*(1-m))/ ((1-cc*pi) *(1-(cc+c2)*pi))+
(1-pi)*(1-cc)*pi *(1-cc-c2-c3)^2*pi^2 * (m^2*(1-m)+m*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2+c3)*pi) *(1-(cc+c2)*pi))
,

#2210, 2201   mean(t %in% c("2210", "2201")) *2
(1-pi)*(1-cc)*pi *(1-cc-c2)*pi *(1-cc-c2-c3)*pi * 2*(m^2*(1-m)+m*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2+c3)*pi) *(1-(cc+c2)*pi))
,

#######

#0011   mean(t %in% c("0011")) *2
(1-pi)^2*c3*pi / ((1-cc*pi) *(1-(cc+c2)*pi))+
(1-pi)^2*(1-cc-c2-c3)^2*pi^2 * (m^2+(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2+c3)*pi) *(1-(cc+c2)*pi))
,

#1100 mean(t %in% c("1100")) *2 
(1-pi)^2*(1-cc)*pi*(1-cc-c2)*pi *(m^2+(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c2+c3)*pi) )
,



#1010,1001   mean(t %in% c("1010","1001")) *2 
(1-pi)^2*(1-cc)*pi*(1-cc-c2-c3)*pi *2*(m^2+(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c2+c3)*pi) )
,


#0110,0101   mean(t %in% c("0110","0101")) *2 
(1-pi)^2*(1-cc-c2)*pi*(1-cc-c2-c3)*pi *2*(m^2+(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c2+c3)*pi) )
,

#######

#0012   mean(t %in% c("0012")) *2 
(1-pi)^2*(1-cc-c2-c3)^2*pi^2 *2*(m*(1-m))/ ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c2+c3)*pi) )
,

#1200  mean(t %in% c("1200")) *2 
(1-pi)^2*(1-cc)*pi*(1-cc-c2)*pi *2*(m*(1-m))/ ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c2+c3)*pi) )
,

#0102,0120 mean(t %in% c("0102","0120")) *2
(1-pi)^2*(1-cc-c2)*pi*(1-cc-c2-c3)*pi *4*(m*(1-m))/ ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c2+c3)*pi) )
,


#1002,1020  mean(t %in% c("1002","1020")) *2
(1-pi)^2*(1-cc)*pi*(1-cc-c2-c3)*pi *4*(m*(1-m))/ ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c2+c3)*pi) )
,


#######

#0001,0010  mean(t %in% c("0001","0010"))*2  

(1-pi)^3*(1-cc-c2-c3)*pi *2/ ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c2+c3)*pi) )
,


#1000, mean(t %in% c("1000" ))*2 
(1-pi)^3*(1-cc)*pi / ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c2+c3)*pi) )
,

# 0100 mean(t %in% c("0100" ))*2 
(1-pi)^3*(1-cc-c2)*pi / ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c2+c3)*pi) )
,

#######

#0000  mean(t %in% c("0000")) 
(1-pi)^4/ ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c2+c3)*pi) )

)
}


###############outcomes in topology 2

categ.full.top2<-list(

#1111   
c("1111","2222")
,
#######
#1112,1121  
c("1112","1121","2221","2212")
,
#1211 
c("1211","2122")
,
#2111 
c("2111","1222")
,
#######
#1122  
c("1122","2211")
,
#1212,1221 
 c("1212","1221","2121","2112")
,
#######
#1110,1101  
 c("1110","1101","2220","2202")
,
#0111,  
 c("0111","0222")
,
#1011   
 c("1011","2022")
,
#######

#1012,1021 
 c("1012","1021","2021","2012")
,
#0112,0121 
c("0112","0121","0221","0212")
,
#1201,1210,2101,2110 
 c("1201","1210","2101","2110","2102","2120","1202","1220")
,
#0122   
c("0122","0211")
,
#1022 
c("1022","2011")
,
#2210, 2201  
 c("2210", "2201","1120","1102")
,
#######
#0011  
 c("0011","0022")
,
#1100 
 c("1100","2200")
,

#1010,1001   
c("1010","1001","2020","2002")
,
#0110,0101   
c("0110","0101","0220","0202")
,
#######

#0012   
 c("0012","0021")
,
#1200  
 c("1200","2100")
,
#0102,0120 
 c("0102","0120","0201","0210")
,

#1002,1020  
c("1002","1020","2001","2010")
,

#######
#0001,0010  
 c("0001","0010","0002","0020")
,

#1000, 
 c("1000" ,"2000")
,
# 0100 
 c("0100","0200" )
,
#######

#0000
 c("0000"))
 
 
 outcome4.full.top2<-function(t)
{out<-NULL
for (i in 1:length(t))
out<-rbind(out,grepl(t[i],categ.full.top2))
out
}



lik4.full.top2<-function(tum,cc,c2,c3,pi,m=0.5)
{if (cc<0 | c2<0 |c3<0 | cc+c2>=1 | cc+c2+c3>=1) return(-10000)
else {
t<-apply(tum,1,paste,collapse="")
probs<-prob4.full.top2(pi,cc,c2,c3,m)
ev<-outcome4.full.top2(t)
max(sum(ev*log(probs)),-10000)
}}

   
   
##############generate tumors from topology 1
tumors4top1<-function(p,cc=0,c2=0,c3=0,m=0.5,rettype="matrix")
{
J<-length(p)
if (cc==0) t0<-rep(0,J) else 
{
t0<-runif(J)
t0[t0<=p*cc]<-1
t0[t0<1]<-0
t0[t0==1 & runif(J)<=m]<-2
}
t1<-t0
t2<-t0
t3<-t0 
t4<-t0

if (c2==0 )
{
x<-runif(J)
t2[t2==0 & x<=(1-cc)*p*m/(1-cc*p)]<-2
t2[t2==0 & x>(1-cc)*p*m/(1-cc*p)& x<=(1-cc)*p/(1-cc*p)]<-1
x<-runif(J)
t1[t1==0 & x<=(1-cc)*p*m/(1-cc*p)]<-2
t1[t1==0 & x>(1-cc)*p*m/(1-cc*p)& x<=(1-cc)*p/(1-cc*p)]<-1
}


else 
{

x<-runif(J)
t2[t2==0 & x<=c2*p*m/(1-cc*p)]<-2
t2[t2==0 & x>c2*p*m/(1-cc*p)& x<=c2*p/(1-cc*p)]<-1
t1<-t2

x<-runif(J)
t2[t2==0 & x<=(1-cc-c2)*p*m/((1-(c2+cc)*p))]<-2
t2[t2==0 & x>(1-cc-c2)*p*m/((1-(c2+cc)*p))& x<=(1-cc-c2)*p/((1-(c2+cc)*p))]<-1
x<-runif(J)
t1[t1==0 & x<=(1-cc-c2)*p*m/((1-(c2+cc)*p))]<-2
t1[t1==0 & x>(1-cc-c2)*p*m/((1-(c2+cc)*p))& x<=(1-cc-c2)*p/((1-(c2+cc)*p))]<-1


}


if (c3==0 )
{
x<-runif(J)
t4[t4==0 & x<=(1-cc)*p*m/(1-cc*p)]<-2
t4[t4==0 & x>(1-cc)*p*m/(1-cc*p)& x<=(1-cc)*p/(1-cc*p)]<-1
x<-runif(J)
t3[t3==0 & x<=(1-cc)*p*m/(1-cc*p)]<-2
t3[t3==0 & x>(1-cc)*p*m/(1-cc*p)& x<=(1-cc)*p/(1-cc*p)]<-1
}


else 
{

x<-runif(J)
t4[t4==0 & x<=c3*p*m/(1-cc*p)]<-2
t4[t4==0 & x>c3*p*m/(1-cc*p)& x<=c3*p/(1-cc*p)]<-1
t3<-t4

x<-runif(J)
t4[t4==0 & x<=(1-cc-c3)*p*m/((1-(c3+cc)*p))]<-2
t4[t4==0 & x>(1-cc-c3)*p*m/((1-(c3+cc)*p))& x<=(1-cc-c3)*p/((1-(c3+cc)*p))]<-1
x<-runif(J)
t3[t3==0 & x<=(1-cc-c3)*p*m/((1-(c3+cc)*p))]<-2
t3[t3==0 & x>(1-cc-c3)*p*m/((1-(c3+cc)*p))& x<=(1-cc-c3)*p/((1-(c3+cc)*p))]<-1


}

if (rettype=="matrix") return(cbind(t1,t2,t3,t4))
else return(paste(t1,t2,t3,t4,sep=""))
}




####probabilies  for topology 1

prob4.full.top1<-function(pi,cc,c2,c3,m)
{cbind(

#1111   mean(t %in% c("1111"))*2 
cc*pi+c2*pi*c3*pi * (m^2+(1-m)^2)/ ((1-cc*pi) )+
c2*pi*(1-cc-c3)^2*pi^2 * (m^3+(1-m)^3)/ ((1-cc*pi)*(1-(cc+c3)*pi) )+
      c3*pi*(1-cc-c2)^2*pi^2 * (m^3+(1-m)^3)/ ((1-cc*pi)*(1-(cc+c2)*pi) )+
(1-cc-c2)^2*pi^2*(1-cc-c3)^2*pi^2 * (m^4+(1-m)^4)/ ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c3)*pi))
,
#######
#1112,1121  mean(t %in% c("1112","1121"))*2 

c2*pi*(1-cc-c3)^2*pi^2 * 2*(m^2*(1-m)+m*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c3)*pi) )+
(1-cc-c2)^2*pi^2 *(1-cc-c3)^2*pi^2 * 2*(m^3*(1-m)+m*(1-m)^3)/ ((1-cc*pi)*(1-(cc+c3)*pi) *(1-(cc+c2)*pi))
,

#1211,2111  mean(t %in% c("1211","2111")) *2
c3*pi*(1-cc-c2)^2*pi^2 * 2*(m^2*(1-m)+m*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2)*pi) )+
(1-cc-c2)^2*pi^2 *(1-cc-c3)^2*pi^2 * 2*(m^3*(1-m)+m*(1-m)^3)/ ((1-cc*pi)*(1-(cc+c3)*pi) *(1-(cc+c2)*pi))
,

#######
#1122  mean(t %in% c("1122")) *2
c2*pi*c3*pi*2*(m*(1-m))/ ((1-cc*pi))+
c3*pi*(1-cc-c2)^2*pi^2 * (m^2*(1-m)+m*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2)*pi) )+
c2*pi*(1-cc-c3)^2*pi^2 * (m^2*(1-m)+m*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c3)*pi) )+
(1-cc-c2)^2*pi^2 *(1-cc-c3)^2*pi^2 * (m^2*(1-m)^2+m^2*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c3)*pi) *(1-(cc+c2)*pi))
,


#1212,1221   mean(t %in% c("1212","1221")) *2
(1-cc-c2)^2*pi^2 *(1-cc-c3)^2*pi^2 * 2*(m^2*(1-m)^2+m^2*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c3)*pi) *(1-(cc+c2)*pi))
,



#######
#1110,1101   mean(t %in% c("1110","1101"))*2
c2*pi*(1-pi)*(1-cc-c3)*pi * 2*(m^2+(1-m)^2)/ ((1-cc*pi)*(1-(cc+c3)*pi) )+
(1-cc-c2)^2*pi^2 *(1-pi)*(1-cc-c3)*pi * 2*(m^3+(1-m)^3)/ ((1-cc*pi)*(1-(cc+c3)*pi)*(1-(cc+c2)*pi) )
,

#0111,1011   mean(t %in% c("0111","1011")) *2
c3*pi*(1-pi)*(1-cc-c2)*pi *2* (m^2+(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2)*pi) )+
(1-cc-c3)^2*pi^2 *(1-pi)*(1-cc-c2)*pi * 2*(m^3+(1-m)^3)/ ((1-cc*pi)*(1-(cc+c3)*pi)*(1-(cc+c2)*pi) )
,


#######

#1012,1021,0112,0121 mean(t %in% c("1012","1021","0112","0121")) *2
(1-cc-c3)^2*pi^2 *(1-pi)*(1-cc-c2)*pi * 4*(m^2*(1-m)+m*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c3)*pi)*(1-(cc+c2)*pi) )
,


#1201,2101,1210,2110 mean(t %in% c("1201","2101","1210","2110"))  *2
(1-cc-c2)^2*pi^2 *(1-pi)*(1-cc-c3)*pi * 4*(m^2*(1-m)+m*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c3)*pi)*(1-(cc+c2)*pi) )
,

#0122,1022 mean(t %in% c("0122","1022")) *2
c3*pi*(1-pi)*(1-cc-c2)*pi * 4*(m*(1-m))/ ((1-cc*pi)*(1-(cc+c2)*pi) )+
(1-cc-c3)^2*pi^2 *(1-pi)*(1-cc-c2)*pi * 2*(m^2*(1-m)+m*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c3)*pi)*(1-(cc+c2)*pi) )
,

#2210, 2201   mean(t %in% c("2210", "2201")) *2
c2*pi*(1-pi)*(1-cc-c3)*pi * 4*(m*(1-m))/ ((1-cc*pi)*(1-(cc+c3)*pi) )+
(1-cc-c2)^2*pi^2 *(1-pi)*(1-cc-c3)*pi * 2*(m^2*(1-m)+m*(1-m)^2)/ ((1-cc*pi)*(1-(cc+c3)*pi)*(1-(cc+c2)*pi) )
,

#######

#0011   mean(t %in% c("0011")) *2
(1-pi)^2*c3*pi / ((1-cc*pi)*(1-(cc+c2)*pi) )+
(1-pi)^2*(1-cc-c3)^2*pi^2 *(m^2+(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c3)*pi) )
,

#1100 mean(t %in% c("1100")) *2 
(1-pi)^2*c2*pi / ((1-cc*pi)*(1-(cc+c3)*pi) )+
(1-pi)^2*(1-cc-c2)^2*pi^2 *(m^2+(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c3)*pi) )
,



#1010,1001,0110,0101   mean(t %in% c("1010","1001","0110","0101")) *2 
(1-pi)^2*(1-cc-c2)*pi*(1-cc-c3)*pi *4*(m^2+(1-m)^2)/ ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c3)*pi) )
,

#######

#0012   mean(t %in% c("0012")) *2 
(1-pi)^2*(1-cc-c3)^2*pi^2 *2*(m*(1-m))/ ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c3)*pi) )
,

#1200  mean(t %in% c("1200")) *2 
(1-pi)^2*(1-cc-c2)^2*pi^2 *2*(m*(1-m))/ ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c3)*pi) )
,

#0102,0120,1002,1020  mean(t %in% c("0102","0120","1002","1020")) *2
(1-pi)^2*(1-cc-c2)*pi*(1-cc-c3)*pi *2*4*(m*(1-m))/ ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c3)*pi) )
,

#######

#0001,0010  mean(t %in% c("0001","0010"))*2  

(1-pi)^3*(1-cc-c3)*pi*2 / ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c3)*pi) )
,


#1000,0100 mean(t %in% c("1000","0100" ))*2 
(1-pi)^3*(1-cc-c2)*pi*2 / ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c3)*pi) )
,

#######

#0000  mean(t %in% c("0000")) 
(1-pi)^4/ ((1-cc*pi)*(1-(cc+c2)*pi)*(1-(cc+c3)*pi) )

)
}



categ.full.top1<-list(
c("1111","2222"   ),
c("1112","1121","2221","2212"),
c("1211","2111","2122","1222") 
,
#######
c("1122","2211") 
,
c("1212","1221","2121","2112")  
,
#######
c("1110","1101","2220","2202")  
,
c("0111","1011","0222","2022")  
,
#######
c("1012","1021","0112","0121","2021","2012","0221","0212") 
,
c("1201","2101","1210","2110","2102","1202","2120","1220") 
,
c("0122","1022","0211","2011")
,
c("2210"," 2201","1120","1102")
,
#######

c("0011","0022")  
,
c("1100","2200")
,
c("1010","1001","0110","0101","2020","2002","0220","0202")
,
#######
c("0012","0021") 
,
c("1200","2100")
,
c("0102","0120","1002","1020","0201","0210","2001","2010")
,
#######
c("0001","0010","0002","0020")
,
c("1000","0100","2000","0200")
,
#######
c("0000")
)


outcome4.full.top1<-function(t)
{out<-NULL
for (i in 1:length(t))
out<-rbind(out,grepl(t[i],categ.full.top1))
out
}


lik4.full.top1<-function(tum,cc,c2,c3,pi,m=0.5)
{if (cc<0 | c2<0 |c3<0 | cc+c2>=1 | cc+c3>=1) return(-10000)
else {
t<-apply(tum,1,paste,collapse="")
probs<-prob4.full.top1(pi,cc,c2,c3,m)
ev<-outcome4.full.top1(t)
max(sum(ev*log(probs)),-10000)
}}

lik4.full.top1fun<-function(cvec,tum,pi,m=0.5){lik4.full.top1(tum,cvec[1],cvec[2],cvec[3],pi,m=0.5)}
lik4.full.top2fun<-function(cvec,tum,pi,m=0.5){lik4.full.top2(tum,cvec[1],cvec[2],cvec[3],pi,m=0.5)}

   

permutations4<-NULL
for (i1 in c(1:4))
for (i2 in c(1:4))
for (i3 in c(1:4))
for (i4 in c(1:4))
if (length(unique(c(i1,i2,i3,i4)))==4) permutations4<-rbind(permutations4,c(i1,i2,i3,i4))

permutations4.top1<-permutations4[c(1,3,5),]
permutations4.top2<-permutations4[seq(1,24,2),]



#creating reference distribution for the test of 4 tumors
LR4.ref<-function(pi,m=0.5,nsim=100,error=FALSE)
{

J<-length(pi)
likcombined<-rep(NA,nsim)
for (sim in 1:nsim)
{
 cat(".")
tstar<-tumors4top1(pi,cc=0,c2=0,c3=0,m=m)
maxlik<-nulllik<-rep(NA,15)

for (i in 1:3)
{
o<-optim(c(0,0,0),fn=lik4.full.top1fun,tum=tstar[, permutations4.top1[i,]],pi=pi,m=m,control=list(fnscale=-1))
maxlik[i]<-o$value
nulllik[i]<-lik4.full.top1fun(c(0,0,0),tum=tstar[, permutations4.top1[i,]],pi=pi,m=m)
}

for (i in 1:12)
{
o<-optim(c(0,0,0),fn=lik4.full.top2fun,tum=tstar[, permutations4.top2[i,]],pi=pi,m=m,control=list(fnscale=-1))
maxlik[i+3]<-o$value
nulllik[i+3]<-lik4.full.top2fun(c(0,0,0),tum=tstar[, permutations4.top2[i,]],pi=pi,m=m)
}

likcombined[sim]<-max(maxlik-nulllik)
}
likcombined
}



LR4.test<-function(tum,pi,m=0.5,ref=NULL,nref=100)
{


maxlik<-nulllik<-rep(NA,15)
mx<--10000
for (i in 1:3)
{
o<-optim(c(0,0,0),fn=lik4.full.top1fun,tum=tum[, permutations4.top1[i,]],pi=pi,m=m,control=list(fnscale=-1))
maxlik[i]<-o$value
nulllik[i]<-lik4.full.top1fun(c(0,0,0),tum=tum[, permutations4.top1[i,]],pi=pi,m=m)
if (maxlik[i]-nulllik[i]>mx) {mx<-mx1<-maxlik[i]-nulllik[i]; o1perm<-permutations4.top1[i,]; o1<-o}

}


mx<--10000
for (i in 1:12)
{
o<-optim(c(0,0,0),fn=lik4.full.top2fun,tum=tum[, permutations4.top2[i,]],pi=pi,m=m,control=list(fnscale=-1))
maxlik[i+3]<-o$value
nulllik[i+3]<-lik4.full.top2fun(c(0,0,0),tum=tum[, permutations4.top2[i,]],pi=pi,m=m)
if (maxlik[i+3]-nulllik[i+3]>mx) {mx<-mx2<-maxlik[i+3]-nulllik[i+3]; o2perm<-permutations4.top2[i,]; o2<-o}

}

LR<-max(maxlik-nulllik)

if (is.null(ref)) ref<-LR4.ref(pi,m,nref)
if (mx2>mx1)   out<-list(LR, mean(LR<=ref),paste("c=",o2$par[1],", c_",o2perm[2],o2perm[3],o2perm[4],"=",o2$par[2],", c_",o2perm[3],o2perm[4],"=",o2$par[3],sep=""), paste("topology 2, tumor order:",paste(o2perm, collapse=" "))) 
else out<-list(LR, mean(LR<=ref),paste("c=",o1$par[1],", c_",o1perm[1],o1perm[2],"=",o1$par[2],", c_",o1perm[3],o1perm[4],"=",o1$par[3],sep=""), paste("topology 1, tumor order:",paste(o1perm, collapse=" "))) #paste(o1perm[1],", ",o1perm[2],"; ",o1perm[3],sep=""))

 names(out)<-c("logLR","p-value","C.hats","topology")
 
 top<-cbind(c(rep(1,3),rep(2,12)), rbind(permutations4.top1,permutations4.top2),maxlik-nulllik )
 top<-as.data.frame(top)
 names(top)<-c("topology","order1","order2","order3", "order4", "logLR")
     #(c(LR, mean(LR<=ref),mx1,o1$par,o1perm,mx2,o2$par,o2perm,maxlik,nulllik))
out<-list(out,top)
names(out)<-c("results","logLR_for_all_topologies")
out
}







##### 3 tumors

 permutations3<-NULL
for (i1 in c(1:3))
for (i2 in c(1:3))
for (i3 in c(1:3))
if (length(unique(c(i1,i2,i3)))==3) permutations3<-rbind(permutations3,c(i1,i2,i3))
permutations3<-permutations3[c(1,2,6),] 

categ3.full<-list(
c("000" ),
c("001","002"),
c("100","200","010","020") 
,
c("110","220") 
,
c("120","210")  
,
c("101","011","202","022")  
,
c("102","012","201","021")  
,
c("112","221") 
,
c("121","211","212","122") 
,
c("111","222")
)


outcome3.full<-function(t)
{out<-NULL
for (i in 1:length(t))
out<-rbind(out,grepl(t[i],categ3.full))
out
}

prob3.full<-function(pi,cc,c2,m)
{cbind(   ((1-pi)^3)/((1-cc*pi)*(1-(cc+c2)*pi)),
((1-pi)^2* (1-cc)*pi)/((1-cc*pi)*(1-(cc+c2)*pi)),
2*((1-pi)^2* (1-cc-c2)*pi)/((1-cc*pi)*(1-(cc+c2)*pi)),
((1-pi)*c2*pi/(1-cc*pi))  +   ((1-pi)* (1-cc-c2)^2*pi^2*(m^2+(1-m)^2))/((1-cc*pi)*(1-(cc+c2)*pi)),
(1-pi)* (1-cc-c2)^2*pi^2*(2*m*(1-m))/((1-cc*pi)*(1-(cc+c2)*pi)),
(1-cc)*pi*(1-pi)* (1-cc-c2)*pi*2*(m^2+(1-m)^2)/((1-cc*pi)*(1-(cc+c2)*pi)),
(1-cc)*pi*(1-pi)* (1-cc-c2)*pi*2*(2*m*(1-m))/((1-cc*pi)*(1-(cc+c2)*pi)),
((1-cc)*pi*c2*pi*2*m*(1-m)/(1-cc*pi))  +   ((1-cc)*pi* (1-cc-c2)^2*pi^2*(m^2*(1-m)+m*(1-m)^2))/((1-cc*pi)*(1-(cc+c2)*pi)),

  ((1-cc)*pi* (1-cc-c2)^2*pi^2*2*(m^2*(1-m)+m*(1-m)^2))/((1-cc*pi)*(1-(cc+c2)*pi)),
  
cc*pi+((1-cc)*pi*c2*pi*(m^2+(1-m)^2)/(1-cc*pi))  + 
 ((1-cc)*pi* (1-cc-c2)^2*pi^2*(m^3+(1-m)^3))/((1-cc*pi)*(1-(cc+c2)*pi))
)
}



lik3.full<-function(tum,cc,c2,pi,m=0.5)
{if (cc<0 | c2<0  | cc+c2>=1 ) return(-10000)
else {
t<-apply(tum,1,paste,collapse="")
probs<-prob3.full(pi,cc,c2,m)
ev<-outcome3.full(t)
max(sum(ev*log(probs)),-10000)
}}

 tumors3<-
function(p,cc=0,c2=0,m=0.5)
{
J<-length(p)
if (cc==0) t0<-rep(0,J) else 
{
t0<-runif(J)
t0[t0<=p*cc]<-1
t0[t0<1]<-0
t0[t0==1 & runif(J)<=m]<-2
}
t1<-t0
t2<-t0
t3<-t0 

if (c2==0 )
{
x<-runif(J)
t2[t2==0 & x<=(1-cc)*p*m/(1-cc*p)]<-2
t2[t2==0 & x>(1-cc)*p*m/(1-cc*p)& x<=(1-cc)*p/(1-cc*p)]<-1
x<-runif(J)
t1[t1==0 & x<=(1-cc)*p*m/(1-cc*p)]<-2
t1[t1==0 & x>(1-cc)*p*m/(1-cc*p)& x<=(1-cc)*p/(1-cc*p)]<-1
}


else 
{

x<-runif(J)
t2[t2==0 & x<=c2*p*m/(1-cc*p)]<-2
t2[t2==0 & x>c2*p*m/(1-cc*p)& x<=c2*p/(1-cc*p)]<-1
t1<-t2

x<-runif(J)
t2[t2==0 & x<=(1-cc-c2)*p*m/((1-(c2+cc)*p))]<-2
t2[t2==0 & x>(1-cc-c2)*p*m/((1-(c2+cc)*p))& x<=(1-cc-c2)*p/((1-(c2+cc)*p))]<-1
x<-runif(J)
t1[t1==0 & x<=(1-cc-c2)*p*m/((1-(c2+cc)*p))]<-2
t1[t1==0 & x>(1-cc-c2)*p*m/((1-(c2+cc)*p))& x<=(1-cc-c2)*p/((1-(c2+cc)*p))]<-1


}

x<-runif(J)
t3[t3==0 & x<=(1-cc)*p*m/(1-cc*p)]<-2
t3[t3==0 & x>(1-cc)*p*m/(1-cc*p)& x<=(1-cc)*p/(1-cc*p)]<-1



return(cbind(t1,t2,t3))
}



LR3.ref<-function(pi,m=0.5,nsim=100)  #assume 12 more clonal
{ref2<-NULL
for (sim in 1:nsim)
{  cat(".")
tstar<-tumors3(pi,cc=0,c2=0,m=m)
maxlik<-nulllik<-rep(NA,3)
mx<--10000
for (i in 1:3)
{
o<-optim(c(0,0),fn=lik3.full.fun,tum=tstar[, permutations3[i,]],pi=pi,m=0.5,control=list(fnscale=-1))
maxlik[i]<-o$value
nulllik[i]<-lik3.full.fun(c(0,0),tum=tstar[, permutations3[i,]],pi=pi,m=m)
}

ref2<-c(ref2,max(maxlik-nulllik))
}
ref2
}



LR3.test<-function(tum,pi,m=0.5,ref=NULL,nref=100)
{
maxlik<-nulllik<-rep(NA,3)
mx<--10000
for (i in 1:3)
{
o<-optim(c(0,0),fn=lik3.full.fun,tum=tum[, permutations3[i,]],pi=pi,m=0.5,control=list(fnscale=-1))
maxlik[i]<-o$value
nulllik[i]<-lik3.full.fun(c(0,0),tum=tum[, permutations3[i,]],pi=pi,m=m)

if (maxlik[i]-nulllik[i]>mx) {mx<-maxlik[i]-nulllik[i]; o1perm<-permutations3[i,]; o1<-o}

}
LR<-max(maxlik-nulllik)

if (is.null(ref)) ref<-LR3.ref(pi,m,nref)
out<-list(LR, mean(LR<=ref),paste("c=",o1$par[1],", c_",o1perm[1],o1perm[2],"=",o1$par[2],sep=""), paste("topology 1, tumor order:",paste(o1perm, collapse=" "))) #paste(o1perm[1],", ",o1perm[2],"; ",o1perm[3],sep=""))

 names(out)<-c("logLR","p-value","C.hats","topology")
 
 top<-cbind(permutations3,maxlik-nulllik )
 top<-as.data.frame(top)
 names(top)<-c("order1","order2","order3","logLR")

out<-list(out,top)
names(out)<-c("results","logLR_for_all_topologies")
out
}

lik3.full.fun<-function(cvec,tum,pi,m=0.5){lik3.full(tum,cvec[1],cvec[2],pi,m=0.5)}








if (ncol( LOHtable)-1!=length(ptlist))  stop("Unequal number of tumors and patient labels. First column of LOHtable has to be the marker name.")

if (!all(as.matrix(  LOHtable[,-1]) %in% c(noloh,loh1,loh2,NA)))
  {print(table(   as.matrix(  LOHtable[,-1])))
  stop("Unrecognized symbols in  LOHtable: it should include only 
  symbols noloh for non-informative markers and loh1 or loh2 for LOH\n") 
   }
   
   
   if (!is.null(refLOHtable) & !is.null(pfreq))
  warning("Both refLOHtable and pfreq are provided - will use pfreq only")
if (!is.null(refLOHtable) & is.null(pfreq))
  if (!all(as.matrix(  refLOHtable[,-1]) %in% c(noloh,loh1,loh2,NA)))
  {print(table(   as.matrix(  refLOHtable[,-1])))
  stop("Unrecognized symbols in  refLOHtable: it should include only 
  symbols noloh for non-informative markers and loh1 or loh2 for LOH\n")
  }

if (!is.null(refLOHtable) & is.null(pfreq))
    {if (!all(refLOHtable[,1]==LOHtable[,1]) | nrow(LOHtable)!=nrow(refLOHtable))
       stop("refLOHtable have to have the same number of markers (rows)
        as LOHtable and the markers in the first column should be in the same order\n")
    pfreq<-apply(refLOHtable[,-1]!=noloh,1,mean,na.rm=TRUE)
    }
if (is.null(refLOHtable) & is.null(pfreq))
    pfreq<-apply(LOHtable[,-1]!=noloh,1,mean,na.rm=TRUE)
if (!is.null(pfreq) & length(pfreq)!=nrow(LOHtable))
 stop("pfreq must have as many frequencies as there are markers (rows) in LOHtable")
if (any(pfreq<0.05)) 
 {warning("Values of pfreq below 0.05 are substituted with 0.05") 
 pfreq[pfreq<0.05]<-0.05
 }
  if (any(pfreq>0.95)) 
 {warning("Values of pfreq above 0.95 are substituted with 0.95") 
 pfreq[pfreq>0.95]<-0.95
 }

 samnms<-names(LOHtable)[-1]
  cat("Calculating reference distribution, might take time...\n")  
 if (any(table(ptlist)==3)) ref3<-LR3.ref(pfreq,m,Nsim)
  if (any(table(ptlist)==4)) ref4<-LR4.ref(pfreq,m,Nsim)
  
out<-list()
nms<-unique(ptlist)
cat("Testing clonality for patient ")    
for (i in 1:length(nms))
{w<-which(ptlist==nms[i])
if (length(w)<3 | length(w)>4) cat (paste("\n skipping patient ",unique(ptlist)[i]," who doesn't have 3 or 4 tumors, ",sep=""))
else  
{cat (paste(unique(ptlist)[i],", ",sep=""))
 if (length(w)==3) out[[i]]<-LR3.test(LOHtable[,1+w],pfreq,m=m,ref=ref3)
 if (length(w)==4) out[[i]]<-LR4.test(LOHtable[,1+w],pfreq,m=m,ref=ref4)
 }
 }
 names(out)<-nms

out
}

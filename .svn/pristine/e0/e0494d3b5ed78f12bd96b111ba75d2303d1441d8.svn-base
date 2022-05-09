grantLR <-
function(tum1,tum2,p,pi,pc,cvalue=0.5,rescale=TRUE,prnfile=NULL,nm=NULL,chrlist)
{
######## calculates likelihood ratio as in the paper
######## input is G/L/N status for each chromosome and both tumors, marginal frequencies,density estimates###for individual comparisons
unlogit<- function(x){exp(x)/(1+exp(x))}
logit<- function(x){log(x/(1-x))}


pi[is.na(pi)]<-1
pc[is.na(pc)]<-1
pi<-as.numeric(pi)
pc<-as.numeric(pc)

tum1<-factor(tum1,levels=c("Gain","Loss","Normal"))
tum2<-factor(tum2,levels=c("Gain","Loss","Normal"))

nchr<-length(tum1)

PG<-p[,1]
PL<-p[,2]
PN<-p[,3]

if (rescale) #rescales the frequencies to conform to number of gains and losses at particular pair of tumors
{
currfr<-(table(tum1)+table(tum2))/(2*nchr)
w<-currfr==0
if (any(w)) currfr<-(currfr+0.05)/sum(currfr+0.05)

if (currfr[1]!=0.5) PGp<-unlogit(logit(PG)*logit( currfr[1])/( mean(logit(PG)))) else PGp<-PG* currfr[1]/ mean(PG)
if (currfr[2]!=0.5) PLp<-unlogit(logit(PL)*logit(currfr[2])/( mean(logit(PL)))) else PLp<-PL* currfr[2]/ mean(PL)


w<-PLp+PGp
PLp[w>1]<-PLp[w>1]/(1.05*w[w>1])
PGp[w>1]<-PGp[w>1]/(1.05*w[w>1])
PNp<-1-PGp-PLp
PG<-PGp
PL<-PLp
PN<-pmin(1,PNp)

}


if (any(round(PG+PL+PN,8)!=1 | round(PG,8)<=0 |round(PL,8)<=0 |round(PN,8)<=0))
 stop("trouble with p's")

pc[tum1!=tum2 | tum1=="Normal" | tum2=="Normal"]<-1
pi[tum1!=tum2 | tum1=="Normal" | tum2=="Normal"]<-1


Rgg<-(tum1==tum2 & tum1=="Gain")
Rll<-(tum1==tum2 & tum1=="Loss")
Rnn<-(tum1==tum2 & tum1=="Normal")
Rgl<-(tum1=="Loss" & tum2=="Gain") | (tum2=="Loss" & tum1=="Gain")
Rgn<-(tum1=="Normal" & tum2=="Gain") | (tum2=="Normal" & tum1=="Gain")
Rln<-(tum1=="Loss" & tum2=="Normal") | (tum2=="Loss" & tum1=="Normal")


cc<-0
LI<-(cc*PG+(1-cc)^2*PG^2/(1-cc*PG-cc*PL))^Rgg  * 
 (cc*PL+(1-cc)^2*PL^2/(1-cc*PG-cc*PL))^Rll  * 
 (2*(1-cc)^2*PG*PL/(1-cc*PG-cc*PL))^Rgl  *
  (2*(1-cc)*PG*PN/(1-cc*PG-cc*PL))^Rgn *  
  (2*(1-cc)*PL*PN/(1-cc*PG-cc*PL))^Rln * (PN^2/(1-cc*PG-cc*PL))^Rnn


cc<-cvalue
LC<-(cc*PG+(1-cc)^2*PG^2/(1-cc*PG-cc*PL))^Rgg  * 
 (cc*PL+(1-cc)^2*PL^2/(1-cc*PG-cc*PL))^Rll  * 
 (2*(1-cc)^2*PG*PL/(1-cc*PG-cc*PL))^Rgl  *
(2*(1-cc)*PG*PN/(1-cc*PG-cc*PL))^Rgn *  (2*(1-cc)*PL*PN/(1-cc*PG-cc*PL))^Rln * 
(PN^2/(1-cc*PG-cc*PL))^Rnn

Bg<-cc*PG/(cc*PG+(1-cc)^2*PG^2/(1-cc*PG-cc*PL))
Bl<-cc*PL/(cc*PL+(1-cc)^2*PL^2/(1-cc*PG-cc*PL))

Bg[tum1!="Gain" |tum2!="Gain" ]<-0.5
Bl[tum1!="Loss" |tum2!="Loss" ]<-0.5


BBg<-Bg*pc+(1-Bg)*pi
BBl<-Bl*pc+(1-Bl)*pi
BBg[tum1!="Gain" |tum2!="Gain" ]<-1
BBl[tum1!="Loss" |tum2!="Loss" ]<-1

a<-(BBg*BBl/pi)[pi!=1]
names(a)<-chrlist[pi!=1]

if (!is.null(prnfile)) 
{write.table(cbind(nm,chrlist,PG,PL,PN,tum1,tum2,LI, LC, BBg*BBl/pi),
file=prnfile,append=TRUE)
}
return(c(prod(LC/LI),prod(LC*BBg*BBl/(LI*pi)),a))

}


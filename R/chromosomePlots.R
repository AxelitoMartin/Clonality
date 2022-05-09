chromosomePlots <-
function(data.seg1,ptlist,ptname,nmad)
{ colorS<-function(x)
{
x[x=="Normal"]<-"grey"
x[x=="Gain"]<-"blue"
x[x=="Loss"]<-"red"
return(x)
}

samnms<-names(data.seg1$data)[-c(1,2)]
chrlist<-unique(data.seg1$data$chrom)
data.seg1GL<-GL(data.seg1,nmad)
gainthres<-data.seg1GL[[2]]
lossthres<-data.seg1GL[[3]]

datseg2<-subset(data.seg1,sample=samnms[ptlist==ptname])

ns2<-dim(datseg2$data)[2]-2
samnms2<-names(datseg2$data)[-c(1,2)]


upper<-gainthres[ptlist==ptname]
lower<-lossthres[ptlist==ptname]
   int.dev <- dev.interactive()
     parask <- par("ask")
    if (int.dev & !parask )
        par(ask = TRUE)

par(mfrow=c(ns2,1))
for (chr in chrlist)
{
for (p1 in c(1:ns2))
{
a<-subset(datseg2,sample=samnms2[p1],chrom=chr)
plot(a$data[,3],pch=".",cex=2,col="grey",ylab="log-ratio",main=samnms2[p1],ylim=c(-max(abs(range(a$data[,3],na.rm=T))),max(abs(range(a$data[,3],na.rm=T)))) )
abline(h=median( subset(datseg2,sample=samnms2[p1])$data[,3],na.rm=T))
abline(h=c(upper[p1],lower[p1]),lty=2)
for (k in 1:nrow(a$output))
lines(which(!is.na(a$data[,3]))[max(1,(cumsum(a$output$num.mark)[k-1]+1),na.rm=T):(cumsum(a$output$num.mark)[k])],rep(a$output$seg.mean[k],a$output$num.mark[k]),
      col=colorS(rep(a$output[k,7],a$output$num.mark[k])),lwd=4)
      
}
mtext(paste("Patient",ptname,"-",chr), side = 1, line = 2)

}
}


histogramPlot <-
function(ptLRvec,refLRvec)
{
par(lwd=2)
a<-hist(log(refLRvec),xlim=c(min(c(log(refLRvec),log(ptLRvec)),na.rm=TRUE),
max(c(log(refLRvec),log(ptLRvec)) ,na.rm=TRUE)),breaks=30,
main=paste("Reference distribution of logLR (black), tested pairs (red)"),xlab="")
b<-hist( log(ptLRvec),plot=FALSE)
mult<-max(1,round(max(a$counts)/max(b$counts)))
a<-hist(rep(log(ptLRvec),mult),border="red",col="red",lwd=10,add=TRUE,density=5,
breaks=length(a$breaks))
axis(side=4,at=sort(unique(a$counts)),labels=as.character(sort(unique(a$counts))/mult),
adj=0,col="red",col.axis="red")
}


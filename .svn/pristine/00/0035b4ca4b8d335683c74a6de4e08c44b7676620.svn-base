prettyplot <-
function(datseg,path,nm,lab.general="",t1lab="Tumor 1",t2lab="Tumor 2",
many=TRUE,gains=NULL,losses=NULL)
{#by default several plots are in one pdf file, pdf and dev.off outside
#xlab=nm; mtext on top = lab.general
 

########### color scheme for gain/loss/normal
colorS<-function(x)
{
x[x=="Normal"]<-"grey"
x[x=="Gain"]<-"blue"
x[x=="Loss"]<-"red"
return(x)
}


chrl<-sort(unique(datseg$data[,1]))
chrful<-as.numeric(substr(chrl,4,5))
samnms<-names(datseg$data)[-c(1:2)]
ns<-dim(datseg$data)[2]-2
if (!many) pdf(paste(path,nm,".pdf",sep=""),height=8,width=11)
if (lab.general=="") lab.general<-nm


seg1<-subset(datseg, sample=samnms[1])
seg2<-subset(datseg, sample=samnms[2])



r1<-range(seg1$data[,3],na.rm=TRUE)
r2<-range(seg2$data[,3],na.rm=TRUE)
nn<-length(seg1$data[,3])/2
ntot<-length(seg1$data[,3])
bet<-r2[2]-r2[1]
cntr1<-bet+abs(r1[1])
cntr2<-bet-abs(r2[2])

plot(0,0,ylim=c(0,r1[2]-r1[1]+r2[2]-r2[1]),xaxt="n",xlab=paste(nm),
    ylab="LogRatio",yaxt="n",xlim=c(1,length(seg1$data[,3])),col=0)
axis(at=c(cntr1+0.8*r1[1],cntr1,cntr1+r1[2],cntr2+r2[1],cntr2,cntr2+0.8*r2[2]),
    labels=round(c(0.8*r1[1],0,r1[2],r2[1],0,0.8*r2[2]) ,2),side=2,cex.axis=0.6)
mtext(lab.general, side = 3, line = 0,cex=1)
lines(rep(1,2),c(0,r1[2]-r1[1]+r2[2]-r2[1]+0.1),col="grey")


topgl<-abs(bet-0.9*max(abs(min(cntr1+seg1$output$seg.mean)),
      max(seg2$output$seg.mean+cntr2)))


locs<-cumsum(table(seg1$data[,1]))
lens<-table(seg1$data[,1])/2
for (chr in c(1:22))
{
if (sum(chrful==chr)==2)
    {lines(rep(locs[which(chrful==chr)[2]],2),c(0,r1[2]-r1[1]+r2[2]-r2[1]+0.1),
        col="grey")
     lines(rep(locs[which(chrful==chr)[1]],2),c(0,r1[2]-r1[1]+r2[2]-r2[1]+0.1),
        col="grey",lty=3)
    text(locs[which(chrful==chr)[1]],rep(0.01-0.01*(chr %%2==0),length(chrl)),
        chr,cex=0.7,pos=1)
    }
else
    { lines(rep(locs[which(chrful==chr)[1]],2),c(0,r1[2]-r1[1]+r2[2]-r2[1]+0.1),
       col="grey")
    text(locs[which(chrful==chr)[1]]-table(seg1$data[,1])[which(chrful==chr)[1]]/2,
        rep(0.01-0.01*(chr %%2==0),length(chrl)),chr,cex=0.7,pos=1)
    }
}


ppl<-cntr1

points(ppl+seg1$data[,3],pch=".",col="darkgrey",cex=2)


points(c(1:ntot)[!is.na(seg1$data[,3])],
      rep(ppl+seg1$output$seg.mean,seg1$output$num.mark),
      col=colorS(rep(seg1$output[,7],seg1$output$num.mark)),pch=".",cex=2)
#abline(h=ppl,lwd=0.5)
abline(h=ppl+median(seg1$data[,3],na.rm=TRUE),lwd=0.5,lty=2)
seg1<-seg1$output

text(nn,r1[2]-r1[1]+r2[2]-r2[1]-(r1[2]-r1[1]+r2[2]-r2[1])/40,t1lab,cex=0.7,pos=2)
ppl<-cntr2

points(ppl+seg2$data[,3],pch=".",col="darkgrey",cex=2)

points(c(1:ntot)[!is.na(seg2$data[,3])],  
      rep(ppl+seg2$output$seg.mean,seg2$output$num.mark),
      col=colorS(rep(seg2$output[,7],seg2$output$num.mark)),pch=".",cex=2)
#abline(h=ppl,lwd=0.5)
abline(h=ppl+median(seg2$data[,3],na.rm=TRUE),lwd=0.5,lty=2)
seg2<-seg2$output

text(nn,bet-(r1[2]-r1[1]+r2[2]-r2[1])/40,t2lab,cex=0.7,pos=2)
abline(h=bet)
if (!is.null(gains[1]))
{arrows(c(1:sum(seg1$num.mark)),rep(bet,length(gains)),
    c(1:sum(seg1$num.mark)),bet+gains*topgl/max(gains),length=0,col="red")
arrows(c(1:sum(seg1$num.mark)),rep(bet,length(gains)),
    c(1:sum(seg1$num.mark)),bet-losses*topgl/max(losses),length=0,col="blue")
}
if (!many) dev.off()
}


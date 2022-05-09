genomewidePlots <-
function(data.seg1,classall,ptlist,ptpair,ptLR,plot.as.in.analysis=TRUE)
{   

  nicen<-function(x)
{#rounding of LR
if (x<100) round(x,1) else formatC(x,2)
}
 if ((!is.character(data.seg1$data$chrom) | !all(nchar(data.seg1$data$chrom==6)) | ! all(substr(data.seg1$data$chrom,1,3) %in% c("Chr","chr")) | ! all(substr(data.seg1$data$chrom,6,6) %in% c("p","q"))))
 stop("Chromosome should have format 'chr' + two digits + 'p or q', like in 'chr01p','Chr11q' etc")


samnms<-names(data.seg1$data)[-c(1,2)]
chrlist<-unique(data.seg1$data$chrom)
nchr<-length(chrlist)

datseg<-subset(data.seg1,sample=ptpair)
dsamnms<-ptpair

if (plot.as.in.analysis)  #if true only segments that classify chromosomes are colored; if false the one step CBS results are plotted
{ 
tum1<-classall[,samnms==ptpair[1]]
tum2<-classall[,samnms==ptpair[2]]
for (ch in c(1:nchr))
{
 datseg$output[datseg$output$ID==dsamnms[1] & datseg$output$chrom==chrlist[ch] &
     datseg$output[,7]!=tum1[ch],7]<-"Normal"
 datseg$output[datseg$output$ID==dsamnms[2] & datseg$output$chrom==chrlist[ch] &
     datseg$output[,7]!=tum2[ch],7]<-"Normal"
qq<-datseg$output[datseg$output$ID==dsamnms[2] & datseg$output$chrom==chrlist[ch],]
nq<-nrow(qq)
if (nq>1 & (all(qq[,7]=="Gain") | all(qq[,7]=="Loss")))
{
w<-which(datseg$output$ID==dsamnms[2] & datseg$output$chrom==chrlist[ch])
ll<-qq[1,]
ll$num.mark<-sum(qq$num.mark)
ll$seg.mean<-sum(qq$num.mark*qq$seg.mean)/ sum(qq$num.mark)
datseg$output[w[1],]<-ll
datseg$output<-datseg$output[-c(w[2]:w[nq]),]
}

qq<-datseg$output[datseg$output$ID==dsamnms[1] & datseg$output$chrom==chrlist[ch],]
nq<-nrow(qq)
if (nq>1 & (all(qq[,7]=="Gain") | all(qq[,7]=="Loss")))
{
w<-which(datseg$output$ID==dsamnms[1] & datseg$output$chrom==chrlist[ch])
ll<-qq[1,]
ll$num.mark<-sum(qq$num.mark)
ll$seg.mean<-sum(qq$num.mark*qq$seg.mean)/ sum(qq$num.mark)
datseg$output[w[1],]<-ll
datseg$output<-datseg$output[-c(w[2]:w[nq]),]
}
}
}

line<-paste("Patient ",ptlist[samnms==ptpair[1]],sep="")
if (!is.null( ptLR))
{
numb<-as.numeric(ptLR[ptLR[,1]==dsamnms[1] & ptLR[,2]==dsamnms[2],4])
if (numb<=1) line<-paste(line,", odds in favor of independence = ",nicen(1/numb),sep="") else 
line<-paste(line,", odds in favor of clonality (metastasis) = ",nicen(numb),sep="")
}
prettyplot(datseg,path="",lab.general=line,nm="",t1lab=paste("Sample",dsamnms[1]),
      t2lab=paste("Sample",dsamnms[2]))

}


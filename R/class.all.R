class.all <-
function(seg1)
{
############ classifies all chomosome in 2 tumors based on one-step CBS pattern: by middle segment if there
## are 3 segments, by most outstanding segment if there are 2 segments.
chrlist<-sort(unique(seg1$output$chrom))
samnms<- names(seg1$data)[-c(1,2)]
ns<-length(samnms)
tum<-matrix("Normal",nrow=length(chrlist),ncol=ns)
for (chr in c(1:length(chrlist)))
for (pt in 1:ns)
{
maplocs<-seg1$data$maploc[seg1$data$chrom==chrlist[chr]]
ss1<-subset(seg1,chrom=chrlist[chr],sample=samnms[pt])
s1<-ss1$output
n1<-nrow(s1)
w1<-which(s1[,7]!="Normal")
if (n1==1) tum[chr,pt]<-s1[,7]
else if (length(w1)==0 | (n1==3 & s1[2,7]=="Normal")) tum[chr,pt]<-"Normal"
else
    {if (n1==2 & length(w1)==1) ind1<-w1
    if (n1==2 & length(w1)==2 & abs(s1$seg.mean[w1[1]])>=abs(s1$seg.mean[w1[2]])) 
          ind1<-w1[1]
    if (n1==2 & length(w1)==2 & abs(s1$seg.mean[w1[1]])<abs(s1$seg.mean[w1[2]]))
          ind1<-w1[2]
    if (n1==3) ind1<-2
    tum[chr,pt]<-s1[ind1,7]
    }
}
colnames(tum)<-samnms
rownames(tum)<- chrlist
return(tum)
}


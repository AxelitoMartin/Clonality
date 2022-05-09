LOHclonality <-
function(LOHtable,ptlist,refLOHtable=NULL, pfreq=NULL,noloh,loh1,loh2,method="both")
{ 
nm<-nrow( LOHtable)
if (length(ptlist)!=ncol( LOHtable)-1) stop("ptlist must have the same number of
 samples as LOHtable. First column in LOHtable must be the marker list.")
if (is.null(names(LOHtable)))   names(LOHtable)<-c("marker",ptlist)
if (!all(as.matrix(  LOHtable[,-1]) %in% c(noloh,loh1,loh2,NA)))
  {print(table(   as.matrix(  LOHtable[,-1])))
  stop("Unrecognized symbols in  LOHtable: it should include only 
  symbols noloh for non-informative markers and loh1 or loh2 for LOH\n") 
   }
if (method=="both" | method=="LR")
{   
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
}
 samnms<-names(LOHtable)[-1]
res<-NULL
cat("Testing clonality for patient ")    
for (i in unique(ptlist))
{w<-which(ptlist==i)
cat (paste(i,", ",sep=""))
ns<-length(w)
for (p1 in c(1:(ns-1)))
for (p2 in c((p1+1):ns))
{
v1<-LOHtable[,1+w[p1]]
v2<-LOHtable[,1+w[p2]]

notna<-!is.na(v1) & !is.na(v2)
v1<-v1[notna]
v2<-v2[notna]
    
rw<-c(sum(v2==v1 & v1!=noloh),sum(v2!=noloh & v1!=noloh),
sum(v1==noloh & v2!=noloh),sum(v2==noloh & v1!=noloh),sum(v1==noloh & v2==noloh))

rw2<-matrix(0,nrow=5,ncol=length(v1))
for (j in c(1:length(v1)))
{if (v1[j]==noloh & v2[j]==noloh) rw2[5,j]<-1
if (v1[j]==noloh & v2[j]!=noloh) rw2[4,j]<-1
if (v1[j]!=noloh & v2[j]==noloh) rw2[3,j]<-1
if (v1[j]!=noloh & v2[j]!=noloh) rw2[2,j]<-1
if (v1[j]!=noloh & v2[j]!=noloh & v2[j]==v1[j]) rw2[1,j]<-1
}
if (method=="both")
res<-rbind(res,c(samnms[w[c(p1,p2)]], rw,sum(notna) , 
CM.pvalue(rw), LRpv(rw2,pfreq[notna],1000)))
else    if (method=="CM")
res<-rbind(res,c(samnms[w[c(p1,p2)]], rw,sum(notna) , 
CM.pvalue(rw)))
else    if (method=="LR")
res<-rbind(res,c(samnms[w[c(p1,p2)]], rw,sum(notna) , LRpv(rw2,pfreq[notna],1000)))
}
}
cat(" Done \n")
res<-as.data.frame(res)
if (method=="both") 
names(res)<-c("Sample1","Sample2","a","e","f","g","h","Ntot","CMpvalue","LRpvalue")
else   if (method=="CM")
names(res)<-c("Sample1","Sample2","a","e","f","g","h","Ntot","CMpvalue")
else if (method=="LR") 
names(res)<-c("Sample1","Sample2","a","e","f","g","h","Ntot","LRpvalue")
res
}


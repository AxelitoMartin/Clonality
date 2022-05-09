complexts <-
function(s1,s2,maplocs)
{
###### function calculating test statistic t - closeness between two concordant segments of gain/loss
###### input - two one-step CBS segmentations
###### maplocs are used to take into account missing values

maplocs<-as.numeric(as.character(maplocs))
s1$loc.start<-as.numeric(as.character(s1$loc.start))
s1$loc.end<-as.numeric(as.character(s1$loc.end))
s2$loc.start<-as.numeric(as.character(s2$loc.start))
s2$loc.end<-as.numeric(as.character(s2$loc.end))

i1<-match(round(s1$loc.start,2),round(maplocs,2))
j1<-match(round(s1$loc.end,2),round(maplocs,2))

i2<-match(round(s2$loc.start,2),round(maplocs,2))
j2<-match(round(s2$loc.end,2),round(maplocs,2))

n1<-nrow(s1)
n2<-nrow(s2)
w1<-which(s1[,7]!="Normal")
w2<-which(s2[,7]!="Normal")

if (n1==1 | n2==1 | length(w1)==0 | length(w2)==0 | all(s1[,7]=="Gain") | 
  all(s1[,7]=="Loss") |  all(s2[,7]=="Gain") | all(s2[,7]=="Loss") ) return(NA)
else  
  if ((n1==3 & s1[2,7]=="Normal") | (n2==3 & s2[2,7]=="Normal")) return(NA)
else 
  if (n1==2 & n2==2 & 
  sign(s1$seg.mean[1]-s1$seg.mean[2])!=sign(s2$seg.mean[1]-s2$seg.mean[2])) 
      return(NA)
else
{
if (n1==2 & length(w1)==1) ind1<-w1
if (n1==2 & length(w1)==2 & abs(s1$seg.mean[w1[1]])>abs(s1$seg.mean[w1[2]])) 
      ind1<-w1[1]
if (n1==2 & length(w1)==2 & abs(s1$seg.mean[w1[1]])<abs(s1$seg.mean[w1[2]])) 
      ind1<-w1[2]
if (n1==2 & length(w1)==2 & abs(s1$seg.mean[w1[1]])==abs(s1$seg.mean[w1[2]]) &
       s1$num.mark[w1[1]]>=s1$num.mark[w1[2]]) ind1<-w1[1]
if (n1==2 & length(w1)==2 & abs(s1$seg.mean[w1[1]])==abs(s1$seg.mean[w1[2]]) & 
        s1$num.mark[w1[1]]<s1$num.mark[w1[2]]) ind1<-w1[2]
if (n1==3) ind1<-2

if (n2==2 & length(w2)==1) ind2<-w2
if (n2==2 & length(w2)==2 & abs(s2$seg.mean[w2[1]])>abs(s2$seg.mean[w2[2]]))
     ind2<-w2[1]
if (n2==2 & length(w2)==2 & abs(s2$seg.mean[w2[1]])<abs(s2$seg.mean[w2[2]])) 
    ind2<-w2[2]
if (n2==2 & length(w2)==2 & abs(s2$seg.mean[w2[1]])==abs(s2$seg.mean[w2[2]]) &
     s2$num.mark[w2[1]]>=s2$num.mark[w2[2]]) ind2<-w2[1]
if (n2==2 & length(w2)==2 & abs(s2$seg.mean[w2[1]])==abs(s2$seg.mean[w2[2]]) & 
    s2$num.mark[w2[1]]<s2$num.mark[w2[2]]) ind2<-w2[2]
if (n2==3) ind2<-2

if (n1==2 & length(w1)==2 & n2==2 & length(w2)==2 ) if ( all(s1[,7]==s2[,7])) 
    ind1<-ind2<-1
if (s1[ind1,7]!=s2[ind2,7] ) return(NA)   # discordnat pairs
else{
overl<-abs(i1[ind1]-i2[ind2])+abs(j1[ind1]-j2[ind2])
nono<-!(i1[ind1]>j2[ind2] | j1[ind1]<i2[ind2])  # is there overlap
b1<-ind1; b2<-ind2
return(c(overl,b1,b2,nono))
    }
}
}


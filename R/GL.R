GL <-
function(seg,nmad)
{
###### Assignment of Gain/Loss/Normal status to each segment based on MAD criteria. 
#Segment if Gain/Loss of its segment mean is above/below nmad MADS of noise level from the median
ns<-dim(seg$data)[2]-2
samnms<-names(seg$data)[-c(1,2)]
gainthres<-rep(NA,ns)
lossthres<-rep(NA,ns)
seg$output[,7]<- "Normal"
for (i in c(1:ns))
  {resid<-seg$data[!is.na(seg$data[,2+i]),2+i]-
  rep(seg$output$seg.mean[seg$output[,1]==samnms[i]],
  seg$output$num.mark[seg$output[,1]==samnms[i]])
  madd<-mad(resid,na.rm=TRUE)
  seg$output[seg$output[,1]==samnms[i]& seg$output$seg.mean>
    median(seg$data[,2+i],na.rm=TRUE)+nmad*madd,7]<-"Gain"
  seg$output[seg$output[,1]==samnms[i] & seg$output$seg.mean<
    median(seg$data[,2+i],na.rm=TRUE)-nmad*madd,7]<-"Loss"
  gainthres[i] <-median(seg$data[,2+i],na.rm=TRUE)+nmad*madd
  lossthres[i] <-median(seg$data[,2+i],na.rm=TRUE)-nmad*madd
  }
names(seg$output)[7]<-"state"
return(list(seg,gainthres,lossthres))
}


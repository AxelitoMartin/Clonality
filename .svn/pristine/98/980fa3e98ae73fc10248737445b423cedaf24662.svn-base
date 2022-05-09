calculateLR <-
function(data.seg1, classall, ptlist, pfreq, reference,
                        allpairs=TRUE, gainthres, lossthres, segmethod,
                        segpar) {
  npts<-length(unique(ptlist))
  samnms<-names(data.seg1$data)[-c(1,2)]
  chrlist<-unique(data.seg1$data$chrom)
  nchr<-length(chrlist)
  testset<-NULL
  if (reference) {
    for (i in (1:(npts-1)))
      for (j in ((i+1):npts)) {
        w1<-which(ptlist==unique(ptlist)[i])
        w2<-which(ptlist==unique(ptlist)[j])
        for (p1 in 1:length(w1))
          for (p2 in 1:length(w2))
            if (allpairs | (!allpairs & p1!=p2))   testset<-rbind(testset,c(w1[p1],w2[p2]))
      }
    if (nrow(testset)<20) warning("too few patients to estimate the reference distribution reliably!! Proceeding anyway...\n")
  } else {       
    for (i in unique(ptlist)) {
      w<-which(ptlist==i) 
      ns<- length(w)
      if (ns>1) {
        for (p1 in c(1:(ns-1)))
          for (p2 in c((p1+1):ns))
            testset<-rbind(testset,c(w[p1],w[p2]))
      }
    }
  }
  if (!reference) cat("Calculating LR")
  else cat("Calculating reference LR: %completed ")
  ncomp<-nrow(testset)
  iLRs<-NULL
  for (i in (1:ncomp)) {
    if (reference & (i %in% round(c(1:10)*ncomp/10))) cat(paste(round(100*i/ncomp),", ",sep=""))
    if (!reference) cat(".")
    x<-rep(NA,nchr)
    y<-rep(NA,nchr)
    for (chr in c(1:nchr)) {
      b<-indiv.test(subset(data.seg1,chrom=chrlist[chr],sample=samnms[testset[i,1]]),
 subset(data.seg1,chrom=chrlist[chr],sample=samnms[testset[i,2]]),func,
 gainthres[testset[i,]],lossthres[testset[i,]],segmethod=segmethod,segpar=segpar)
      if  (!is.na(b[1]) )
{x[chr]<-b[[3]]
y[chr]<-b[[4]]
}
}
tum1<-classall[,testset[i,1]]
tum2<-classall[,testset[i,2]]
a<-grantLR(tum1,tum2,pfreq,x,y,cvalue=0.5,chrlist=chrlist)

pattern<-c(sum(tum1==tum2 & tum1!="Normal"),sum(tum1==tum2 & tum1=="Normal"),
sum(tum1!=tum2 & tum1!="Normal" & tum2!="Normal") )
pattern<-c(pattern,length(chrlist)-sum(pattern))
pattern<-c(pattern,paste(names(a[-c(1,2)]),round(a[-c(1,2)],2),collapse="; "))
iLRs<-rbind(iLRs,c(samnms[testset[i,]],a[c(1,2)],pattern))
}

iLRs<-as.data.frame(iLRs)  
names(iLRs)<-c("Sample1",	"Sample2",	"LR1"	,"LR2",	"GGorLL",	"NN",	"GL",	"GNorLN",	"IndividualComparisons")
iLRs[,4]<-as.numeric(as.character(iLRs[,4]))
iLRs[,3]<-as.numeric(as.character(iLRs[,3]))
cat("\n")
return(iLRs)
}


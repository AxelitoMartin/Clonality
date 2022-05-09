calc.freq <-
function(pfreq, refdata, classall, nmad, segmethod, segpar) {
  chrlist<-rownames(classall)
  if (!is.null(pfreq)) {
    if (ncol(pfreq)!=4) 
      stop("pfreq should have 4 columns: chromosome arm name, frequency of gain, loss, and normal.") 

    if (!is.character(pfreq[,1]) | !all(nchar(pfreq[,1]==6)) | 
        ! all(substr(pfreq[,1],1,3) %in% c("Chr","chr")) | 
        ! all(substr(pfreq[,1],6,6) %in% c("p","q")))
      stop("Chromosome arm in the first column of pfreq should have format 'chr' + two digits + 'p or q', like in 'chr01p','Chr11q' etc")

    pfreq<-pfreq[match(chrlist,as.character(pfreq[,1])),2:4]
    pfreq<-as.data.frame(pfreq)
    pfreq[,1]<-as.numeric(as.character(pfreq[,1]))
    pfreq[,2]<-as.numeric(as.character(pfreq[,2]))
    pfreq[,3]<-as.numeric(as.character(pfreq[,3]))
    if (!is.null(refdata)) cat("Refdata is not used since pfreq is specified")
    return(pfreq)
  }

  if (is.null(pfreq)) {
    if (!is.null(refdata)) {
      cat("Resolution of the reference cohort arrays should be similar to the resolution of the averaged data. Please make sure. Proceeding...\n")
      refseg1<-segment1(refdata[refdata$chrom %in% chrlist,],segmethod=segmethod,segpar=segpar)
      classall<-class.all(GL(refseg1,nmad)[[1]])
    }

    nppl<-ncol(classall)
    if (nppl<10) 
      warning("too few patients to estimate the frequencies reliably!! Proceeding anyway...\n")

    pg<-apply(classall=="Gain",1,mean,na.rm=TRUE)
    pl<-apply(classall=="Loss",1,mean,na.rm=TRUE)

    pg[pg<0.05]<-0.05
    pl[pl<0.05]<-0.05
    pg[pg>0.9]<-0.9
    pl[pl>0.9]<-0.9
    
       w<-which(pg+pl>=0.95)    
   diff<-pg[w]+pl[w]-0.95
   pg[w]<-pg[w]- diff/2
    pl[w]<-pl[w]- diff/2
    
    pfreq<-cbind(pg,pl,1-pg-pl)

    if (any(pfreq<0)) {stop("negative frequencies")}
    return(pfreq)
  }
}


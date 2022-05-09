clonality.analysis <-
function(data, ptlist, pfreq=NULL, refdata=NULL,
                               nmad=1.25, reference=TRUE, allpairs=TRUE,
                               segmethod="oneseg", segpar=NULL) {
  if (!inherits(data, "CNA")) 
    stop("First arg must be a copy number array(CNA) object\n")
  if (segmethod=="oneseg") {
    if (missing(segpar)) {
      segpar$alpha <- 0.01
      segpar$nperm <- 2000
      segpar$sbdry <- get.oneseg.bdry(alpha=0.01, nperm=2000)
    } else {
      if (exists("alpha", segpar) & exists("nperm", segpar)) {
        segpar$sbdry <- get.oneseg.bdry(alpha=segpar$alpha, nperm=segpar$nperm)
      } else {
        stop("segpar missing alpha and/or nperm")
      }
    }
  } else {
    warning(" you have given an alternate segmentation function\n ensure that this function and segpar meet the requirements\n")
  }
  npts <- length(unique(ptlist))
  if (!(any(table(ptlist)>=2)))
    stop("No pairs of tumors from the same patient. Check  ptlist")
  if (npts==1 & reference) {
    warning("Just one patient - can't evaluate reference distribution\n")
    reference <- FALSE
  }
  data <- data[data$chrom!="chrX" & data$chrom!="chrY" & 
               data$chrom!="ChrX" & data$chrom!="ChrY" &
               data$chrom!=23 & data$chrom!=24,]
  if (all(is.numeric(data$chrom)))
    warning("Chromosomes should be split into p and q arms to increase power")
  if (any(table(ptlist)<2)) 
    warning("Warning: some patients have only 1 tumor")
  if (nrow(data)>=15000) 
    warning("Averaging is highly recommended; use ave.adj.probes() function first")
  
  samnms <- names(data)[-c(1,2)]
  chrlist <- unique(data$chrom)
  if (any(table(data$chrom)<5)) {
    cat("Removing the following chromosomes since they have fewer than 5 markers\n")
    cat(paste(names(table(data$chrom))[table(data$chrom)<5],"\n"))
  }
  data <- data[!(data$chrom %in% chrlist[table(data$chrom)<5]),]
  chrlist <- unique(data$chrom)
  nchr <- length(chrlist)
  data.seg1 <- segment1(data,segmethod=segmethod, segpar=segpar)
  data.seg1GL <- GL(data.seg1,nmad)
  data.seg1 <- data.seg1GL[[1]]
  classall <- class.all(data.seg1)
  pfreq <- calc.freq(pfreq, refdata, classall, nmad, segmethod=segmethod, segpar=segpar)
  ptLR <- calculateLR(data.seg1,classall,ptlist,pfreq,reference=FALSE,
                      allpairs=allpairs,gainthres=data.seg1GL[[2]],
                      lossthres=data.seg1GL[[3]],segmethod=segmethod,
                      segpar=segpar)
  if (reference) {
    refLR <- calculateLR(data.seg1,classall,ptlist,pfreq,reference=TRUE,
                         allpairs=allpairs,gainthres=data.seg1GL[[2]],
                         lossthres=data.seg1GL[[3]],segmethod=segmethod,
                         segpar=segpar)
    LR2pvalue<-NULL
    for (i in 1:nrow(ptLR)) 
      LR2pvalue<-c(LR2pvalue,mean(ptLR[i,4]<=refLR[,4],na.rm=TRUE))
    ptLR<-cbind(ptLR, LR2pvalue)
  }

  if (reference) return(list("LR"=ptLR,"OneStepSeg"=data.seg1,"ChromClass"=classall,"refLR"=refLR))
  else return(list("LR"=ptLR,"OneStepSeg"=data.seg1,"ChromClass"=classall))
}


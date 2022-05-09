segment1 <-
function(xcna, segmethod, segpar) {
  ######One step CBS - finds at most one (most prominent) copy number change 
  #on each  chromosome arm of every sample
  if (!inherits(xcna, "CNA"))
    stop("First arg must be a copy number array object")
  
  segres <- list()
  segres$data <- xcna
  outp<-NULL
  for (sam in c(1:(ncol(xcna)-2))) {
    for (chr in unique(xcna$chrom)) {
      if (!all(is.na(xcna[xcna$chrom==chr,2+sam]))) {
        mapl <- xcna$maploc[xcna$chrom==chr]
        x <- xcna[xcna$chrom==chr,2+sam]
        if (any(is.na(x))) {mapl<-mapl[!is.na(x)];  x<-x[!is.na(x)]}
        n <- length(mapl)

        segdata <- c(list(x=x), segpar)
        seg <- do.call(segmethod, segdata)
        if (seg[1]==0) 
          outp<-rbind(outp,c(names(xcna)[2+sam],chr,mapl[1],mapl[n],n,mean(x,na.rm=TRUE)))
        if (seg[1]==1) {
          outp<-rbind(outp,c(names(xcna)[2+sam],chr,mapl[1],
                             mapl[seg[3]],seg[3],mean(x[1:seg[3]],na.rm=TRUE)))
          outp<-rbind(outp,c(names(xcna)[2+sam],chr,mapl[seg[3]+1],
                             mapl[n],n-seg[3],mean(x[(seg[3]+1):n],na.rm=TRUE)))
        }
        if (seg[1]==2) {
          outp<-rbind(outp,c(names(xcna)[2+sam],chr,mapl[1],mapl[seg[2]-1],
                             seg[2]-1,mean(x[c(1:(seg[2]-1),(seg[3]+1):n)],na.rm=TRUE)))
          outp<-rbind(outp,c(names(xcna)[2+sam],chr,mapl[seg[2]],mapl[seg[3]],
                             seg[3]-seg[2]+1,mean(x[seg[2]:seg[3]],na.rm=TRUE)))
          outp<-rbind(outp,c(names(xcna)[2+sam],chr,mapl[seg[3]+1],mapl[n],
                             n-seg[3],mean(x[c(1:(seg[2]-1),(seg[3]+1):n)],na.rm=TRUE)))
        }
      }
    }
  }
  outp<-as.data.frame(outp)
  names(outp)<-c( "ID",  "chrom", "loc.start", "loc.end" ,"num.mark", "seg.mean")
  outp$seg.mean<-round(as.numeric(as.character(outp$seg.mean)),5)
  outp$num.mark<-as.numeric(as.character(outp$num.mark))
  outp$loc.start<-as.numeric(as.character(outp$loc.start))
  outp$loc.end<-as.numeric(as.character(outp$loc.end))
  segres$output <- outp
  class(segres) <- "DNAcopy"
  segres
}


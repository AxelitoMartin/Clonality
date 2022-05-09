ave.adj.probes <-
function(data, K) {
  all.ave <- NULL
  chrs <- NULL
  i <- 1
  while (i < nrow(data)) {
    n <- sum(data$chrom[i:(i+K-1)]==data$chrom[i])
    if (!is.na(n) & n==K) {
      chrs <- c(chrs,data$chrom[i+round(K/2)-1])
      all.ave <- rbind(all.ave,(apply(data[i:(i+K-1),-1],2,mean,na.rm=TRUE)))
      i <- i+K
    }
    else if (!is.na(n)) i <- i+n
    else i <- nrow(data)+1
  }
  data <- CNA(as.matrix(all.ave[,-1]),chrom=chrs,maploc=as.numeric(all.ave[,1]),sampleid=names(data)[-c(1,2)])
  cat(paste("Total number of markers after averaging is",nrow(data),"\n"))
  return(data)
}


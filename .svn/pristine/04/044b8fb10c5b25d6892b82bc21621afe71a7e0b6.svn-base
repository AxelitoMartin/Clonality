\name{splitChromosomes}
\alias{splitChromosomes}
\title{
Chromosome splitting
}
\description{
Divides the chromosomes into p and q arms.   }
\usage{
splitChromosomes(chrom,maploc)
}

\arguments{
  \item{chrom}{
Vector of chromosomes. They should be numeric 1 to 22.
}
  \item{maploc}{
Vector of genomic locations. They should be in Kilobases. 
}
}
\details{
The function returns the vector of chromosome arms labeled "chr01p", "chr01q", etc. The split into arms is accomplished using the following  centers (in Kb) for  chromosomes 1 through 22: (122356.96,  93189.90,  92037.54 , 50854.87  ,47941.40,  60438.12 , 59558.27,  45458.05 , 48607.50,  40434.94 , 52950.78,  35445.46 ,      16934.00,  16570.00,  16760.00 , 36043.30 , 22237.13,  16082.90  ,      28423.62 , 27150.40,  11760.00, 12830.00 ). 

}

\examples{
#simulated data
 
set.seed(100)
chrom<-rep(c(1:22),each=100)
maploc<- runif(2200)* 200000
chromarm<-splitChromosomes(chrom,maploc)
 
}
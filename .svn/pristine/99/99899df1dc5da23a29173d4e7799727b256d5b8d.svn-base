 \name{ECMtesting}
\alias{ECMtesting}

\title{
Clonality testing of >=3 tumors using Extended Concordant Mutations (ECM) test based on LOH (Loss of Heterozygosity) profiles
}
\description{
Function to test clonality of three and more tumors from the same patient based on their LOH profiles. This function implements Extended Concordant Mutations for all possible subsets of tumors from the same patient and minP multiplicity adjustment using simulated tumors.
}
\usage{
ECMtesting(LOHtable,ptlist,noloh,loh1,loh2,Nsim=100)
}

\arguments{
  \item{LOHtable}{
Matrix of LOH  calls. Each row corresponds to a marker. First column contains the names of the markers. Each other column represents a sample and contains LOH calls. 
}
  \item{ptlist}{
Vector of the patient IDs in the order the samples appear in the data. For example, if the first three tumors (columns 2, 3, 4 of data) belong to patient A, and the following two (columns 5, 6 of data) belong to patient B, then ptlist=c('ptA', 'ptA', 'ptA', 'ptB', 'ptB').
}
  \item{noloh}{
The string or a number that denotes absence of LOH.
}
  \item{loh1}{
The string or a number that denotes presence of LOH.
}
  \item{loh2}{
The string or a number that denotes presence of LOH that is discordant from loh1.
}
  \item{Nsim}{
Number of simulations used to calculate minP adjusted p-values}
}
\details{
Extended Concordant Mutations test is done for every subset of tumors. It uses number of concordant mutations in all tumors of the subset as a test statistic, and its reference distribution is calculated assuming fixed  counts of LOH per tumor and equal probability of maternal and paternal alleles being affected. Note that ECM test for 2 tumors and original CM test will give slightly different p-values since continuity correction is done in ECM test.



}
\value{
The function returns a list with number of elements equal to the number of patients. Each element is a matrix with two rows: ECM p-values for all possible subsets of tumors from this patient, and minP adjusted p-values. The tumors are denoted 1,2,3,... in the order they appear in LOHtable. Any tumor subsets with minP adjusted p-value <=0.05 can be considered clonal.

}
\references{ Ostrovnaya, I. "Testing clonality of three and more tumors using their loss of heterozygosity profiles", Statistical Applications in Genetics and Molecular Biology, 2012
}
\examples{
set.seed(25)
LOHtable<-cbind(1:15,matrix(sample(c(0,1,2),15*12,replace=TRUE),ncol=12))
ECMtesting(LOHtable,rep(1:3,each=4),noloh=0,loh1=1,loh2=2,Nsim=100)
  }

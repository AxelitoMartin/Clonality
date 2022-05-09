 \name{LRtesting3or4tumors}
\alias{LRtesting3or4tumors}

\title{
Clonality testing of 3 or 4 tumors using Likelihood model based on LOH (Loss of Heterozygosity) profiles
}
\description{
Function to test clonality of 3 or 4 tumors from the same patient based on their LOH profiles.
}
\usage{
LRtesting3or4tumors(LOHtable,ptlist,refLOHtable=NULL, pfreq=NULL,noloh,loh1,loh2,Nsim=100,m=0.5)
}

\arguments{
  \item{LOHtable}{
Matrix of LOH  calls. Each row corresponds to a marker. First column contains the names of the markers. Each other column represents a sample and contains LOH calls. 
}
  \item{ptlist}{
Vector of the patient IDs in the order the samples appear in the data. For example, if the first three tumors (columns 2, 3, 4 of data) belong to patient A, and the following two (columns 5, 6 of data) belong to patient B, then ptlist=c('ptA', 'ptA', 'ptA', 'ptB', 'ptB').
}
  \item{refLOHtable}{
Matrix of LOH  calls that should be used to calculate the LOH frequencies used in Likelihood Ratio calculation. The structure is similar to LOHtable.  If refLOHtable is not specified,  frequencies are calculated from LOHtable.
}
  \item{pfreq}{
Vector of LOH frequencies known from the literature. Should be in the same order as the markers in LOHtable.    If pfreq is not specified,  frequencies are calcualted from LOHtable.
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
Number of simulations used to calculate minP adjusted p-values
}

  \item{m}{
Probability that a favored allele is  affected given that LOH has occurred. Must be a number above 0.5 (equal probability of maternal and paternal allelic loss)  
}
}
\details{
Likelihood ratio test for 3 and 4 tumors. For 3 tumors there are 3 possible tumor orderings, and for 4 tumors there are 2 topologies with 3 and 12 orderings each. The test calculates maximum likelihood ratio across all possible orderings, and the p-value is calculated using simulated reference distribution.



}
\value{
The function returns a list with number of elements equal to the number of patients. Each element is list with two elements. First contains log maximum likelihood ratio value, p-value, and estimates of parameters c, the topology and tumor ordering  that have maximum likelihood ratio.  If p-value is significant, then the null hypothesis that all tumors are independent can be rejected. 
The second element has a matrix with all possible topologies and tumor orderings and their corresponding log likelihood ratios. 

}
\references{ Ostrovnaya, I. "Testing clonality of three and more tumors using their loss of heterozygosity profiles", Statistical Applications in Genetics and Molecular Biology, 2012
}
\examples{
set.seed(25)
LOHtable<-cbind(1:15,matrix(sample(c(0,1,2),15*12,replace=TRUE),ncol=12))
q<-LRtesting3or4tumors(LOHtable,rep(1:4,each=3),refLOHtable=NULL, pfreq=NULL,noloh=0,loh1=1,loh2=2,Nsim=100,m=0.5)
 
   }

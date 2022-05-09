 \name{LOHclonality}
\alias{LOHclonality}

\title{
Clonality testing using LOH (Loss of Heterozygosity) profiles
}
\description{
Function to test clonality of two tumors from the same patient based on their LOH profiles. This function implements Concordant Mutations and Likelihood Ratio tests.
}
\usage{
LOHclonality(LOHtable, ptlist, refLOHtable = NULL, pfreq = NULL, noloh, loh1, loh2,method="both")
}

\arguments{
  \item{LOHtable}{
Matrix of LOH  calls. Each row corresponds to a marker. First column contains the names of the markers. Each other column represents a sample and contains LOH calls. 
}
  \item{ptlist}{
Vector of the patient IDs in the order the samples appear in the data. For example, if the first three tumors (columns 3, 4, 5 of data) belong to patient A, and the following two (columns 6, 7 of data) belong to patient B, then ptlist=c('ptA', 'ptA', 'ptA', 'ptB', 'ptB').
}
  \item{refLOHtable}{
Matrix of LOH  calls that should be used to calculate the LOH frequencies used in Likelihood Ration calculation. The structure is similar to LOHtable.  If refLOHtable is not specified,  frequencies are calculated from LOHtable.
}
  \item{pfreq}{
Vector of LOH frequencies known from the literature. Should be in the same order as the markers in LOHtable.    If pfreq is not specified,  frequencies are calculated from LOHtable.
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
  \item{method}{
Takes values "CM", "LR" or "both" if only Concordant Mutations test, or only Likelihood Ratio test, or both should be performed. Default value is "both".
}
}
\details{
Function tests clonality of LOH profiles of tumors from the same patient using two tests. Concordant Mutations test has number of markers with concordant LOH as its test statistic. Its theoretical reference distribution under independence is calculated assuming that the maternal and paternal alleles are equally likely to be lost and that the frequencies of LOH are about the same across different markers. 

 Likelihood Ratio test uses pre-specified frequencies of LOH to compute Likelihood Ratio statistic. Its reference distribution is obtained by simulating tumors with the given LOH probabilities, and probability of maternal/paternal mutation estimated from the data. If LOH frequencies are not specified then they are estimated from the data.

}
\value{
The function returns a data frame where each row corresponds to the pair of samples that are compared. 
Columns are:
Sample1 - name of sample 1;
Sample2 - name of sample 2;
a - number of markers with concordant LOH in both tumors (test statistic for Concordant Mutations test);
e -  number of markers with  LOH in both tumors, concordant or discordant;
f - number of markers with LOH in the first tumor and Normal in the 2nd tumor;
g - number of markers with LOH in the second tumor and Normal in the first tumor;
h - number of markers that are Normal in both tumors;
Ntot - total number of informative markers for both tumors;
CMpvalue - p-value for Concordant Mutations test;
LRpvalue - p-value for Likelihood Ratio test.

}
\references{

Begg CB, Eng KH, Hummer AJ. Statistical tests for clonality. Biometrics 2007; 63:522-530

Ostrovnaya I, Seshan VE, Begg CB. Comparison of properties of tests for assessing tumor clonality. Biometrics 2008; 68:1018-1022.
}
\examples{
set.seed(25)
LOHtable<-cbind(1:20,matrix(sample(c(0,1,2),20*20,replace=TRUE),20))
LOHclonality(LOHtable,rep(1:10,each=2),pfreq=NULL,noloh=0,loh1=1,loh2=2)
  }

\name{genomewidePlots}
\alias{genomewidePlots}
\title{
Plot of the genomewide copy number profiles of a pair of tumors.
}
\description{
Plot contains genomewide profiles from a pair of tumors. It uses the output from the function clonality.analysis().
}
\usage{
genomewidePlots(data.seg1, classall, ptlist, ptpair, ptLR, plot.as.in.analysis = TRUE)
}

\arguments{
  \item{data.seg1}{
Output of one-step segmentation - output OneStepSeg of clonality.analysis().  The chromosomes should be in the format "chr01p", "chr01q" etc.
}
  \item{classall}{
Classifications of the chromosomes  - output ChromClass of clonality.analysis()
}

  \item{ptlist}{
Vector of the patient IDs in the order of the samples appearing in the data. 
}
  \item{ptpair}{
Two sample names for which the plot is desired}
  \item{ptLR}{
Matrix with the likelihood ratios -  output LR of clonality.analysis()
}
  \item{plot.as.in.analysis}{
If TRUE then the gain/loss  patterns  will be highlighted in accordance with the chromosome classification.  For example, if there are three segments in a chromosome, then the middle one determines the chromosome status. If it is normal, no color will be plotted in the chromosome even if the 1st and 3rd segments are gains or losses. Another example: if there are 2 or 3 different segments of gains, they will be combined and only one segment will be plotted. If plot.as.in.analysis is equal to FALSE, the original one-step CBS segmentation will be plotted.
}
}
\details{
Function produces genomewide plots of a pair of tumors.  The log-ratios are plotted in grey in the order of their genomic locations, gains are plotted in blue, and losses are plotted in red.
}


\examples{
 # See example as in clonality.analysis()
   
}

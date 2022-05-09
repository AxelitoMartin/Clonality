
\name{chromosomePlots}
\alias{chromosomePlots}
\title{
Per-chromosome plots of the copy number arrays  from a particular patient
}
\description{
The function produces a sequence of plots for each chromosome with one-step segmented data  of all samples of a particular patient.
}
\usage{
chromosomePlots(data.seg1, ptlist, ptname,nmad)
}

\arguments{
  \item{data.seg1}{
Output of one-step segmentation - output OneStepSeg of clonality.analysis(). 
}
  \item{ptlist}{
Vector of the patient IDs in the order of the samples appearing in the data.  For example, if the first three tumors  belong to patient A, and the following two belong to patient B, then ptlist=c('ptA', 'ptA', 'ptA', 'ptB', 'ptB'). 
}
  \item{ptname}{
Name of the patient from ptlist for which the data should be plotted}
  \item{nmad}{
Number of MADs (median absolute deviations) that is used for Gain/Loss calls. Used to mark the Gain/Loss threshold on the plots.
}
}
\details{ 
The function produces a sequence of plots  for each chromosome with one-step segmented data of all samples of a particular patient. The dotted horizontal lines denote the gain and loss thresholds.
}

\examples{
 # See example as in clonality.analysis()
   
}

\name{histogramPlot}
\alias{histogramPlot}

\title{
Histrograms of Log-Likelihood Ratios
}
\description{
Function produces the histograms of the within-patient and between-patient log-Likelihood Ratios.
}
\usage{
histogramPlot(ptLRvec, refLRvec)
}

\arguments{
  \item{ptLRvec}{
Vector with the within-patient likelihood ratios -  output LR of clonality.analysis()
}
  \item{refLRvec}{
Vector with the between-patient likelihood ratios -  output refLR of clonality.analysis()
}
}
\details{
Functions plots two overlapping histograms:  within-patient log-likelihood ratios are in red and between-patient log-likelihood ratios (reference distribution under the hypothesis of independence) are in black.
}

\examples{
 # See example as in clonality.analysis()


}
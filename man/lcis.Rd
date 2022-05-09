\name{lcis}
\alias{lcis}

\title{
  Breast cancer data
}
\description{
Data from the LCIS study, with mutations listed for all pairs of LCIS-invasive tumors in a series of breast cancer cases.

}
\usage{
  data(lcis)
}

\details{
  This is exome sequencing data from study of Lobular Carcinoma in Situ (LCIS) and  
  Invaisve lobular carcinomas  (ILC) or Invasive Ductal Carcinomas (IDC) in the same patients.
  Each row represents single mutation in a patient. This structure is similar to MAF (mutation annotation format) files.

The columns are "PatientID, Lesion,Hugo_Symbol, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2 (Alternative allele), Tumor_Sample_Barcode (sample ID). Lesion contains Lesion type, which takes values DCIS (ductal carcinoma in situ), IDC (invasive ductal carcinoma), IDC2 (second profiled IDC), ILC (invasive lobular carcinoma), LCIS and LCIS1,2, etc.
}

\references{
  Begg CB, Ostrovnaya I, Carniello JV, Sakr RA, Giri D, Towers R, Schizas M, De Brot M, Andrade VP, Mauguen A, Seshan VE, King TA. "Clonal relationships between lobular carcinoma in situ and other breast malignancies.", Breast Cancer Res. 2016 Jun 23;18(1):66. doi: 10.1186/s13058-016-0727-z.

}


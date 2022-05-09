\name{get.mutation.frequencies}
\alias{get.mutation.frequencies}

\title{Estimating mutation frequencies based on the TCGA data or a submitted reference mutations file}
\description{
Function to estimate the mutation frequencies in the specific cancer subtype using TCGA data or using  a submitted reference mutations file
}
\usage{
get.mutation.frequencies (xmut.ids, tcga.cancer.type=NULL, reference.data=NULL,combine.with.TCGA=FALSE  )
}

\arguments{
  \item{xmut.ids}{
vector of mutation IDs for which frequencies are needed. Usually these are mutation IDs observed in one or two tumors from the same patient that will be tested for clonality. xmut.ids should have the following format: {Chromosome Location RefAllele AltAllele}, each entry separated by space,  where chromosome is a number 1-22 or X or Y; location is genomic location in  GRCh37 build; RefAllele is a reference allele and AltAllele is Alternative allele/Tumor_Seq_Allele2.  For example '10 100003849 G A' is the mutation at chromosome 10, genomic location 100003849, where reference allele G is substituted with A, or  '10 100011448 - CCGCTGCAAT' is the insertion of 'CCGCTGCAAT' at chromosome 10, location 100011448.  The ref and alt alleles follow standard TCGA maf file notations.
}
  \item{tcga.cancer.type}{
String denoting one of 33 TCGA cancer types:  ACC,  BLCA, BRCA, CESC, CHOL, COAD, DLBC, ESCA, GBM,  HNSC, KICH, KIRC, KIRP, LAML, LGG,  LIHC, LUAD, LUSC, MESO, OV,   PAAD, PCPG, PRAD, READ, SARC, SKCM, STAD, TGCT, THCA, THYM, UCEC, UCS,  or UVM. See https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations for details. If reference.data is supplied, cancer.type can be NULL

}
  \item{reference.data}{
is a matrix that contains the mutations in external dataset from which the frequencies can be estimated. Maf files (Mutation Annotation Format) can be used here.  Each line represents a mutation, and the matrix should include the following columns: 'PatientID','Tumor_Sample_Barcode' (used as a sample ID), and Chromosome, Start_Position, Reference_Allele, and Tumor_Seq_Allele2 that are used to create mutation IDs.

If the frequencies are estimated for a specific patient with mutations 'xmut.ids', then this matrix should not contain the data from this patient. In this case the denominator for the frequency will be (1+n) where n is number of patients in a reference. Alternatively, if the list of mutations in 'xmut.ids' is the same as list of all mutations in reference.data then the denominator will be n. If the reference dataset  contains multiple samples from  the same patient, union of them will be taken so that clonal mutations are not counted twice -  denominator is the number of patients and not the number of samples. 

Another word of caution that if you are planning to compare pairs of tumors for clonality profiled by exome sequencing, the reference dataset should also be from exome sequencing (like in TCGA), to avoid counting mutations in the genes that were not sequenced by targeted panel for example as non-observed/never mutated.

The possible reasons to use reference.data include: 1) the dataset you analyze is large and combining it with tcga data will improve the frequency; 2) you analyze tumors from cancer type that is not among 33 pancancer TCGA sites; or 3) you want to use specific subtype of cancers, for example Luminal A breast tumors, for which overall TCGA breast cancer frequencies are not applicable.
   

}
  \item{combine.with.TCGA}{
TRUE if you want to combine  mutations in reference dataset with TCGA for frequency calculation. In this case you need to both specify 'cancer.type' and 'reference.data'. This option makes sense when the dataset you assembled for clonality study is large and therefore will improve frequency calculation.  Default is FALSE.
   

}
}
\details{
For each mutation in 'xmut.ids' we calculate mutation frequency in the following way. 

If tcga.cancer.type is chosen and reference.data is NULL, we calculate X1, number of patients that have that mutation,  among n1 TCGA samples of specific subtype. Mutation frequency then is (1+X1)/(1+n1), where 1 is added to denominator and nominator to incorporate data in a current patient.

If reference.data is specified and tcga.cancer.type is not chosen, we calculate X2 - number of patients with a particular  mutation among n2 patients (not samples) in reference.data. Then there are two possibilities. If the list of mutations in 'xmut.ids' is the same as list of all mutations in reference.data, we assume the reference data is the dataset for which clonality is under question, and the frequency is defined as X2/n2. Otherwise if 'xmut.ids' is not the same as list of all mutations in reference.data then we assume it's an external data and like in TCGA the frequency then is (1+X2)/(1+n2).

If both tcga.cancer.type and reference.data are specified and  combine.with.TCGA is  TRUE, then the frequency is defined as either (1+X1+X2)/(1+n1+n2) or (X1+X2)/(n1+n2)  depending on whether the list of mutations in 'xmut.ids' is the same as list of all mutations in reference.data.
}
\value{
Vector of frequencies with names same as 'xmut.ids'
}
\references{

Ostrovnaya, Irina, Venkatraman E. Seshan, and Colin B. Begg. 2015. “USING SOMATIC MUTATION DATA TO TEST TUMORS FOR CLONAL RELATEDNESS.” The Annals of Applied Statistics 9 (3): 1533–48. https://doi.org/10.1214/15-AOAS836.

}
\author{
Irina Ostrovnaya \email{ostrovni@mskcc.org}
}

\examples{
data(lcis)

#Analysis of data from patient 53
mut.matrix<-create.mutation.matrix(lcis )
table(mut.matrix$TK53IDC2,mut.matrix$TK53LCIS2)

freq<-get.mutation.frequencies(rownames(mut.matrix),"BRCA") 

SNVtest(mut.matrix$TK53IDC2,mut.matrix$TK53LCIS2,freq)


}
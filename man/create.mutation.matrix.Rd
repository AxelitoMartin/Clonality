\name{create.mutation.matrix}
\alias{create.mutation.matrix}

\title{Formatting matrix of mutations }
\description{This functions reformats matrix of mutations for subsequent analysis

}
\usage{
create.mutation.matrix(maf,rem=FALSE)}

\arguments{
  \item{maf}{
maf is a mutations file in generic MAF (mutation annotation format) style: one row per mutation per sample. See this webpage for detailed description of the format: https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/ .
Object 'maf' should have the following columns: PatientID, Tumor_Sample_Barcode (sample ID), Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2 (reference allele. For compatibility with mutational frequencies that can be obtained using function "get.mutation.frequencies" chromosome should be a number 1-22 or X or Y; Start_Position is genomic location in  GRCh37 build; Reference_Allele is a reference allele and Tumor_Seq_Allele2 is Alternative allele.  Mutation IDs are created, e.g. '10 100003849 G A' is the mutation at chromosome 10, genomic location 100003849, where reference allele G is substituted with A, or  '10 100011448 - CCGCTGCAAT' is the insertion of 'CCGCTGCAAT' at chromosome 10, location 100011448.  The ref and alt alleles follow standard TCGA maf file notations.
}
  \item{rem}{
if TRUE, the mutationan matrix for random effect function will be prepared.

}
}
\details{
if rem=FALSE (default), binary mutational matrix will be created, where each row is a mutation, each column is a sample, and values are binary taking values 0 if there is no mutation, 1 if there is a mutation in this sample. If rem=TRUE matrix with each possible pair of samples from the same patient will be created. Each row represents mutation, and each column - pair of samples, where value 0  denotes that mutation is not 
observed, 1 if shared mutation is observed in both tumors, and 2 if it is a private mutation observed in only one tumor.

}
\value{
Data frame with matrix of mutations
}
\references{

Ostrovnaya, Irina, Venkatraman E. Seshan, and Colin B. Begg. 2015. “USING SOMATIC MUTATION DATA TO TEST TUMORS FOR CLONAL RELATEDNESS.” The Annals of Applied Statistics 9 (3): 1533–48. https://doi.org/10.1214/15-AOAS836.

}
\author{
Irina Ostrovnaya \email{ostrovni@mskcc.org}
}

\examples{
data(lcis)
#we want to analyze pair TK53IDC2.TK53LCIS2 that has only 1 shared mutation

mut.matrix<-create.mutation.matrix(lcis )
table(mut.matrix$TK53IDC2,mut.matrix$TK53LCIS2)

freq<-get.mutation.frequencies(rownames(mut.matrix),"BRCA") 
SNVtest(mut.matrix$TK53IDC2,mut.matrix$TK53LCIS2,freq)


}
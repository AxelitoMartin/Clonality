\name{freqdata}
\alias{freqdata}

\title{
 TCGA pancancer mutation frequencies
}
\description{
Data from the LCIS study, with mutations listed for all pairs of LCIS-invasive tumors in a series of breast cancer cases.

}
\usage{
  data(freqdata)
}

\details{
  Object 'freqdata' contains the frequencies of the mutations observed in the exome sequencing data in 33 cancer types profiled by TCGA. There are 33 columns for 33 cancer types abbreviated in TCGA as GBM, OV, LUAD, etc. (for the full list of abbreviations see https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations). The first column is the number of patients profiled in each cancer type. Each subsequent row is a mutation, where number of patients with this particular mutation in each cancer type is given. The mutation ID, contained in the row names of 'freqdata', is of the following format: {"Chromosome Location RefAllele AltAllele"}, each entry separated by space, where chromosome is a number 1-22 or "X" or "Y"; location is genomic location in  GRCh37 build; RefAllele is a reference allele and AltAllele is Alternative allele/Tumor_Seq_Allele2. For example "10 100003849 G A", is the mutation at chromosome 10, genomic location 100003849, where reference allele G is substituted with A, or "10 100011448 - CCGCTGCAAT" is the insertion of "CCGCTGCAAT" at chromosome 10, location 100011448. The ref and alt alleles follow standard TCGA maf file notations.
  
This data object was obtained using the following code:

#'  The file mc3.v0.2.8.PUBLIC.maf can be downloaded from https://gdc.cancer.gov/about-data/publications/mc3-2017
#' pancan<-read.table(  "mc3.v0.2.8.PUBLIC.maf",sep="\t",header=T,quote="") 
#' 
#'  The file clinical.tsv.csv can be downloaded from https://portal.gdc.cancer.gov/projects - > check TCGA and clinical on the left panel -> 
#' click on 11,315 Cases across 33 Projects on the upper right corner -> click on Clinical->TSV button below
#' pancaninfo<-read.table("clinical.tsv",sep="\t",header=T)
#' 
#' 
#' pancaninfo$project_id<-substr(pancaninfo$project_id,6,11)
#' pancan$Tumor_Sample_Barcode<-substr(pancan$Tumor_Sample_Barcode,1,12)
#' pancan$mut<-paste(pancan$Chromosome,pancan$Start_Position,pancan$Reference_Allele,pancan$Tumor_Seq_Allele2)
#' pancan$type<-pancaninfo$project_id[match(pancan$Tumor_Sample_Barcode,pancaninfo$submitter_id)]
#' pancan<-pancan[!duplicated(paste(pancan$Tumor_Sample_Barcode,pancan$mut)),]
#' pancan<-pancan[!is.na(pancan$type),]
#' 
#' pancan$mut<-as.factor(pancan$mut)
#' types<-unique(pancan$type)
#' ntypes<-NULL
#'freqdata<-NULL
#'for (i in types)
#'{w<-pancan$type==i
#'s<-length(unique(pancan$Tumor_Sample_Barcode[w]))
#' ntypes<-c(ntypes,s)
#' freqdata<-cbind(freqdata,table(pancan$mut[w]))
#' 
#' }
#' 
#' rownames(freqdata)<-names(table(pancan$mut[w]))
#' colnames(freqdata)<-types
#' freqdata<-rbind(ntypes,freqdata)


}
\references{
Ostrovnaya, Irina, Venkatraman E. Seshan, and Colin B. Begg. 2015. “USING SOMATIC MUTATION DATA TO TEST TUMORS FOR CLONAL RELATEDNESS.” The Annals of Applied Statistics 9 (3): 1533–48. https://doi.org/10.1214/15-AOAS836.

Ellrott K, Bailey MH, Saksena G, Covington KR, Kandoth C, Stewart C, Hess J, Ma S, Chiotti KE, McLellan M, Sofia HJ, Hutter C, Getz G, Wheeler D, Ding L; MC3 Working Group; Cancer Genome Atlas Research Network,  Scalable Open Science Approach for Mutation Calling of Tumor Exomes Using Multiple Genomic Pipelines. Cell Syst. 2018 Mar 28;6(3):271-281.e7. doi: 10.1016/j.cels.2018.03.002. PubMed PMID: 29596782.
}


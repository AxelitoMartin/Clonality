get.mutation.frequencies <- function(xmut.ids, tcga.cancer.type=NULL, reference.data=NULL,combine.with.TCGA=FALSE ) {
  library(data.table)
  # library(dplyr)
  # library(dtplyr)
  
  # data("d.gene")
  # data("GT_probs")
  # data("GT_agg_probs")
  
  # d.gene <- data.table::setDT(d.gene)
  # print(d.gene[1,])
  # make these internal data once done #
  # load("/Users/axel/Desktop/MSKCC Work/Clonality Project/temp_clonality/Clonality/data/TCGA_GT_probs.Rdata")
  # load("/Users/axel/Desktop/MSKCC Work/Clonality Project/temp_clonality/Clonality/data/TCGA_GT_agg_probs.Rdata")
  
  xmut.ids <- gsub("X",23,xmut.ids)
  plus1<-TRUE
  # if (! (all(substr(xmut.ids,1,1) %in% c(c(1:9),"X","Y") | substr(xmut.ids,1,2) %in% c(c(10:22)) )))
  #   stop("xmut.ids should be of the following format: {Chromosome Location RefAllele AltAllele}, each entry separated by space, \n where chromosome is a number 1-22 or X or Y; \n location is genomic location in  GRCh37 build; \n RefAllele is a reference allele and AltAllele is Alternative allele/Tumor_Seq_Allele2. \n For example '10 100003849 G A', is the mutation at chromosome 10, genomic location 100003849, where reference allele G is substituted with A, or \n '10 100011448 - CCGCTGCAAT' is the insertion of 'CCGCTGCAAT' at chromosome 10, location 100011448. \n The ref and alt alleles follow standard TCGA maf file notations.")
  
  if (is.null(tcga.cancer.type) & is.null(reference.data)) stop("You either need to specify 'tcga.cancer.type' to use TCGA based frequnecies or specify 'reference.data' to get reference dataset for frequency calculation")
  if (combine.with.TCGA & (is.null(tcga.cancer.type) | is.null(reference.data))) stop("If you choose combine.with.TCGA=TRUE, i.e. to combine reference dataset with TCGA data, then you  need to  specify both'tcga.cancer.type'  and 'reference.data' ")
  if (!combine.with.TCGA & (!is.null(tcga.cancer.type) & !is.null(reference.data))) warning("You choose combine.with.TCGA=FALSE, but specified both TCGA and reference data - only TCGA data will be used unless you choose combine.with.TCGA=TRUE")
  
  # match variants to their genes #
  matched_genes <- as.data.table(do.call('rbind',lapply(xmut.ids, function(x){
    # print(x)
    temp <- strsplit(x, " ")[[1]]
    chr <- as.numeric(temp[1])
    pos <- as.numeric(temp[2])
    # print("test")
    gene <- d.gene$Hugo_Symbol[d.gene$Chromosome == chr & pos > d.gene$Start & pos < d.gene$End]
    if(length(gene)==0)
      gene <- "MissingGeneName"
    return(c(x,gene[1],NA,NA))
  })))
  colnames(matched_genes) <- c("variant","gene","freq","aggregated_estimate")
  # get estimates for GT at that sites #
  site_spe_gt <- as.data.table(GT_probs[[tcga.cancer.type]])
  site_spe_gt <- site_spe_gt[site_spe_gt$warning ==0,]
  site_spe_agg_gt <- as.data.table(GT_agg_probs[[tcga.cancer.type]])
  # get estimates per genes found in the matching #
  genes <- unique(matched_genes$gene)
  
  for(i in 1:length(genes)){
    gene_temp <- genes[i]
    variants <- matched_genes$variant[matched_genes$gene == gene_temp]
    site_spe_gt_temp <- site_spe_gt[site_spe_gt$Gene == gene_temp,]
    if(nrow(site_spe_gt_temp) > 0){
      variants_temp <- ifelse(variants %in% site_spe_gt_temp$Variant, variants, "each_unseen")
      matched_genes[match(variants,matched_genes$variant) , 3] <- as.numeric(site_spe_gt_temp$GT[match(variants_temp,site_spe_gt_temp$Variant)])
      matched_genes[match(variants,matched_genes$variant) , 4] <- FALSE
    }
    else{
      variants_temp <- ifelse(variants %in% site_spe_agg_gt$Variant, variants, "each_unseen")
      matched_genes[match(variants,matched_genes$variant) , 3] <- as.numeric(site_spe_agg_gt$GT_agg[match(variants_temp,site_spe_agg_gt$Variant)])
      matched_genes[match(variants,matched_genes$variant) , 4] <- TRUE
    }
  }
  
  
  
  
  #   if (!is.null(reference.data)) #estimation using reference file
  #   { if (all(colnames(reference.data)!="Chromosome")) stop("reference.data should include column 'Chromosome' that will be used to create mutation IDs")
  #     if (all(colnames(reference.data)!="PatientID")) stop("reference.data should include column 'PatientID' - could be the same as sample IDs Tumor_Sample_Barcode")
  #     if (all(colnames(reference.data)!="Start_Position")) stop("reference.data should include column 'Start_Position' that will be used to create mutation IDs")
  #     if (all(colnames(reference.data)!="Reference_Allele")) stop("reference.data should include column 'Reference_Allele' that will be used to create mutation IDs")
  #     if (all(colnames(reference.data)!="Tumor_Seq_Allele2")) stop("reference.data should include column 'Tumor_Seq_Allele2' that will be used to create mutation IDs")
  #     if (all(colnames(reference.data)!="Tumor_Sample_Barcode")) stop("reference.data should include column 'Tumor_Sample_Barcode' that will be used in frequency calculation")
  #     
  #     reference.data$mut<-paste(reference.data$Chromosome,reference.data$Start_Position,reference.data$Reference_Allele,reference.data$Tumor_Seq_Allele2)
  #     if (any(duplicated(paste(reference.data$mut,reference.data$PatientID))))
  #     {
  #       reference.data<-reference.data[!duplicated(paste(reference.data$mut,reference.data$PatientID)),]
  #     }
  #     
  #     
  #     if (all(xmut.ids %in%  reference.data$mut) & all(reference.data$mut %in% xmut.ids )) 
  #     {#warning("All mutations in xmut.ids are observed in the reference dataset, and vice versa. See help file for more info. ")
  #       plus1<-FALSE
  #     }
  #     else
  #       
  #       n<-length(unique(reference.data$PatientID))  
  #     if (n<50 & is.null(tcga.cancer.type)) warning(paste("reference dataset has only",n,"patients - frequency estimation might be not very accurate")    )
  #     
  #     fref<-NULL
  #     for (i in xmut.ids) fref<-c(fref,length(unique(reference.data$PatientID[reference.data$mut==i])) )
  #     
  #     if (all(fref==0,na.rm=T))    warning("None of the mutations in 'xmut.ids' were seen in this reference dataset. Make sure  'xmut.ids' are in the correct format: check ?get.mutation.frequencies")
  #     if (  plus1==TRUE)
  #       freq<-(fref+1)/(n+1)
  #     else   freq<-fref/n
  #     
  #     
  #     
  #   }
  #   
  #   if (!is.null(tcga.cancer.type)) #estimation using TCGA
  #   {  data(freqdata)
  #     if (!(tcga.cancer.type %in% colnames(freqdata))) 
  #       stop("tcga.cancer.type should be one of 33 TCGA cancer types:  ACC  BLCA BRCA CESC CHOL COAD DLBC ESCA GBM  HNSC KICH KIRC KIRP LAML LGG  LIHC LUAD LUSC
  # MESO OV   PAAD PCPG PRAD READ SARC SKCM STAD TGCT THCA THYM UCEC UCS  UVM \n See https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations for details")
  #     f<-freqdata[match(xmut.ids,rownames(freqdata)), tcga.cancer.type]
  #     names(f)<-xmut.ids
  #     if (all(is.na(f)))    warning("None of the mutations were seen in TCGA.  Make sure  xmut.ids are in the correct format: check ?get.mutation.frequencies")
  #     if (all(f==0,na.rm=T))    warning("None of the mutations were seen in this TCGA subtype. Make sure  xmut.ids are in the correct format: check ?get.mutation.frequencies")
  #     f[is.na(f)]<-0
  #     
  #     
  #     freq<-(f+1)/(freqdata[1, tcga.cancer.type]+1)
  #     
  #   }
  #   
  #   
  #   
  #   if (combine.with.TCGA & !is.null(reference.data) & !is.null(tcga.cancer.type))
  #   {     if (  plus1==TRUE) freq<-(fref+f+1)/(freqdata[1, tcga.cancer.type]+n+1)
  #   else freq<-(fref+f)/(freqdata[1, tcga.cancer.type]+n)
  #   
  #   
  #   }  
  
  # freq <- as.numeric(matched_genes$freq)
  # names(freq)<-xmut.ids
  # freq
  matched_genes
}






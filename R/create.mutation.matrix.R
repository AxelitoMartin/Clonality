create.mutation.matrix<-function(maf,rem=FALSE)
{ if (all(colnames(maf)!="Chromosome")) stop("maf matrix should include column 'Chromosome' that will be used to create mutation IDs")
    if (all(colnames(maf)!="Start_Position")) stop("maf matrix should include column 'Start_Position' that will be used to create mutation IDs")
    if (all(colnames(maf)!="Reference_Allele")) stop("maf matrix should include column 'Reference_Allele' that will be used to create mutation IDs")
    if (all(colnames(maf)!="Tumor_Seq_Allele2")) stop("maf matrix should include column 'Tumor_Seq_Allele2' that will be used to create mutation IDs")
    if (all(colnames(maf)!="Tumor_Sample_Barcode")) stop("maf matrix should include column 'Tumor_Sample_Barcode' that will be used in frequency calculation")
   maf$mut<-paste(maf$Chromosome,maf$Start_Position,maf$Reference_Allele,maf$Tumor_Seq_Allele2)
   
   if (!rem)
    {mut.matrix<-matrix(0,ncol=length(unique(maf$Tumor_Sample_Barcode)),nrow=length(unique(maf$mut)))
    colnames(mut.matrix)<-sort(unique(maf$Tumor_Sample_Barcode))
    rownames(mut.matrix)<-unique(maf$mut)
    for (i in unique(maf$Tumor_Sample_Barcode))
      mut.matrix[maf$mut[maf$Tumor_Sample_Barcode==i],i]<-1
    return(as.data.frame(mut.matrix))
   }
   else 
  {
if (all(colnames(maf)!="PatientID")) stop("maf matrix should include column 'PatientID' to identify paired samples")
  
     mut.matrix<-matrix(0,ncol=0,nrow=length(unique(maf$mut)))
 
  pairnames<-NULL
    rownames(mut.matrix)<-unique(maf$mut)
   for (i in unique(maf$PatientID))
   {
   s<-unique(maf$Tumor_Sample_Barcode[maf$PatientID==i])
   ns<-length(s)  
   for (i1 in 1:(ns-1))
     for (i2 in (i1+1):ns)
     {x<-rep(0,length(unique(maf$mut)))
     x[unique(maf$mut) %in% maf$mut[maf$Tumor_Sample_Barcode==s[i1]] | unique(maf$mut) %in% maf$mut[maf$Tumor_Sample_Barcode==s[i2]]]<-2 
     x[unique(maf$mut) %in% maf$mut[maf$Tumor_Sample_Barcode==s[i1]] & unique(maf$mut) %in% maf$mut[maf$Tumor_Sample_Barcode==s[i2]]]<-1
     mut.matrix<-cbind(mut.matrix,x)
     pairnames<-c( pairnames,paste(s[i1],s[i2],sep="_"))
     }
   }
    colnames(mut.matrix)<-pairnames
   return(as.data.frame(mut.matrix))
    
  }
}
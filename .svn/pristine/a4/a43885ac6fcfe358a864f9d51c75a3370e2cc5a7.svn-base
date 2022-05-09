splitChromosomes <-
function(chrom,maploc)
 { 
 
#locations of the centromeres in Kb
centro<-c(122356.96,  93189.90,  92037.54 , 50854.87  ,47941.40,  60438.12 ,
        59558.27,  45458.05 , 48607.50,  40434.94 , 52950.78,  35445.46 , 
        16934.00,  16570.00,  16760.00 , 36043.30 , 22237.13,  16082.90  ,
        28423.62 , 27150.40,  11760.00, 12830.00 )

if (!all(is.numeric(chrom)) | any(chrom>22))
      stop("Chromosome should be a number between 1 and 22") 

maploc<-maploc[chrom<=22 ]
chrom<-chrom[chrom<=22 ]
kb=1
if (mean(maploc[chrom==1]<centro[1],na.rm=TRUE) < 0.01) 
warning("There are too few markers in one of the arms. The maplocs should be in Kb") 
if (mean(maploc[chrom==1]>centro[1],na.rm=TRUE) < 0.01) 
warning("There are too few markers in one of the arms. The maplocs should be in Kb") 
chr<-paste("chr",chrom,sep="")
chr[chrom<=9]<-paste("chr0",chrom[chrom<=9],sep="")

for (i in c(1:22))
		{w<-chrom==i
		chr[w & maploc<kb*centro[i]]<-
        paste(chr[w & maploc<kb*centro[i]],"p",sep="")
		chr[w & maploc>=kb*centro[i]]<-
        paste(chr[w & maploc>=kb*centro[i]],"q",sep="")
		}
chrom<-chr
print(table(chrom))
return(chrom)
   }


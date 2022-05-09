\name{ave.adj.probes}
\alias{ave.adj.probes}

\title{
Averaging of adjacent probes in copy number arrays
}
\description{
For each sample the log-ratios at each consecutive K number of probes are averaged.
}
\usage{
ave.adj.probes(data, K)
}

\arguments{
  \item{data}{
Copy Number Array object (output of function CNA() from the package DNAcopy). First column contains chromosomes, second column contains genomic locations. Each remaining column contains log-ratios from a particular tumor or sample.
}
  \item{K}{
Number of markers to be averaged. Should be selected so that the final resolution of the averaged data would be 5,000-10,000 markers.
}
}
\details{
Averages log-ratios in every K consecutive markers. The purpose of this step is to reduce the noise in the data, eliminate possible very small germline copy number variations, and get rid of a possible wave effect.
}
\value{
Returns CNA object of reduced resolution
}

          
\examples{
# Same example as in clonality.analysis()

set.seed(100)
chrom<-rep(c(1:22),each=100)
maploc<- runif(2200)* 200000
chromarm<-splitChromosomes(chrom,maploc)
 
 
#Simulate the dataset with 10 pairs of tumors with 22 chromosomes, 100 markers each
#Simulated log-ratios are equal to signal + noise
#Signal: each chromosome has 50% chance to be normal, 30% to be whole-arm loss/gain, and 20% to be partial arm loss/gain, where endpoints are drawn at random, loss/gain means are drawn from normal distribution
#There are no chromosomes with recurrent losses/gains
#Noise: drawn from normal distribution with mean 0, standard deviation 0.25
#First 9 patients have independent tumors, last patient has two tumors with identical signal, independent noise


set.seed(100)
chrom<-paste("chr",rep(c(1:22),each=100),"p",sep="")
chrom[nchar(chrom)==5]<-paste("chr0",substr(chrom[nchar(chrom)==5] ,4,5),sep="")
maploc<- rep(c(1:100),22)
data<-NULL
for (pt in 1:9)  #first 9 patients have independent tumors
{
tumor1<-tumor2<- NULL
mean1<- rnorm(22) 
mean2<- rnorm(22)
for (chr in 1:22)
{ 
  r<-runif(2) 
if (r[1]<=0.5) tumor1<-c(tumor1,rep(0,100))   
  else if   (r[1]>0.7)  tumor1<-c(tumor1,rep(mean1[chr],100))
  else  { i<-sort(sample(1:100,2))
        tumor1<-c(tumor1,mean1[chr]*c(rep(0,  i[1]),rep(1, i[2]-i[1]), rep(0,  100-i[2])))
        }
if (r[2]<=0.5) tumor2<-c(tumor2,rep(0,100))
  else if   (r[2]>0.7)  tumor2<-c(tumor2,rep(mean2[chr],100))
  else   {i<-sort(sample(1:100,2))
       tumor2<-c(tumor2,mean2[chr]*c(rep(0,  i[1]),rep(1, i[2]-i[1]), rep(0,  100-i[2])))
         }
}
data<-cbind(data,tumor1,tumor2)
}

#last patient has identical profiles
tumor1<- NULL
mean1<- rnorm(22) 
for (chr in 1:22)
{ 
  r<-runif(1) 
if (r<=0.4) tumor1<-c(tumor1,rep(0,100))   
  else if   (r>0.6)  tumor1<-c(tumor1,rep(mean1[chr],100))
  else  { i<-sort(sample(1:100,2))
        tumor1<-c(tumor1,mean1[chr]*c(rep(0,  i[1]),rep(1, i[2]-i[1]), rep(0,  100-i[2])))
        }

}
data<-cbind(data,tumor1,tumor1)

data<-data+matrix(rnorm( 44000,mean=0,sd=0.4) ,nrow=2200,ncol=20)
dataCNA<-CNA(data,chrom=chrom,maploc=maploc,sampleid=paste("pt",rep(1:10,each=2),rep(1:2,10)))
dim(dataCNA)
dataCNA2<-ave.adj.probes(dataCNA, 2)
dim(dataCNA2)
}

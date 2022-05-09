ECMtesting <-
function(LOHtable,ptlist,noloh,loh1,loh2,Nsim=100)
{ 
if (ncol( LOHtable)-1!=length(ptlist))  stop("Unequal number of tumors and patient labels. First column of LOHtable has to be the marker name.")

if (!all(as.matrix(  LOHtable[,-1]) %in% c(noloh,loh1,loh2,NA)))
  {print(table(   as.matrix(  LOHtable[,-1])))
  stop("Unrecognized symbols in  LOHtable: it should include only 
  symbols noloh for non-informative markers and loh1 or loh2 for LOH\n") 
   }
#minP ajustment for ECM p-values from all possible subsets of tumors
minpadjECM<-function(tum,Nsim=1000,noloh=0,loh1=1,loh2=2)
{


      
     gen.indep.tumors<-function(tum)
     {tumS<-apply(tum,2,sample)
      tumS[tumS==2]<-1
      tumS[tumS==1 & matrix(runif(ncol(tum)*nrow(tum)),ncol=ncol(tum))<=0.5] <-2
      tumS
      }

         
m<-apply(tum!=noloh, 2, sum)
J<-nrow(tum)
nt=ncol(tum)
pvs<-allECM(tum,noloh,loh1,loh2)
n<-length(pvs)
if (nt>2)
{
vec<-matrix(NA,ncol=n,nrow=Nsim)
for (simul in 1:Nsim)
{
tumS<- gen.indep.tumors(tum)
vec[simul,]<-as.numeric(allECM(tumS,noloh,loh1,loh2))
}
s<-order(pvs)
adjpvs<-(sum(apply(vec,1,min)<as.numeric(sort(pvs)[1]))+0.5*sum(apply(vec,1,min)==as.numeric(sort(pvs)[1])))/Nsim
for (i in 2:(n-1))
{adjpvs<-c(adjpvs,max(adjpvs,(sum(apply(vec[,s[i:n]],1,min)<as.numeric(sort(pvs)[i])) +  0.5*sum(apply(vec[,s[i:n]],1,min)==as.numeric(sort(pvs)[i])))/Nsim    ))
}
adjpvs<-c(adjpvs,max(adjpvs,(sum(vec[,s[n]]<as.numeric(sort(pvs)[n])) +0.5*sum(vec[,s[n]]==as.numeric(sort(pvs)[n])) )/Nsim  ))
out<-rbind(pvs,adjpvs[order(s)])
colnames(out)<-colnames(pvs)
}
else 
{out<-rbind(pvs,pvs)
colnames(out)<-"12"
}
rownames(out)<-c("p","adjusted.p")
out
}


#ECM test for all possible subsets of tumors
allECM<- function(tum,noloh=0,loh1=1,loh2=2)
{n<-ncol(tum)
all<-as.data.frame(0)
for (k in n:2)
{a<-NULL
#create all possible subsets of size k out if n tumors
co<-matrix(c(1:k),nrow=1)
colap<- apply(co,1,paste,collapse=" ")
while (nrow(co)<choose(n,k))
{s<-sort(sample(1:n,k))
if (all(colap!=paste(s,collapse=" ")))
      {co<-rbind(co,s)
        colap<-c(colap,paste(s,collapse=" ") )
      }
}
if (nrow(co)>1) co<-co[sort.list(colap),]

for (i in 1:nrow(co))
{
a<-c(a, ECM.pvalue(tum[,co[i,]],noloh,loh1,loh2 ))
}
a<-as.data.frame(t(a))
colnames(a)<- apply(co,1,paste,collapse="")
all<-cbind(all,a)
}
all[,-1]
}



#ECM test for specific set of tumors
ECM.pvalue<-function (tum,noloh=0,loh1=1,loh2=2)
{
tum<-tum[apply(!is.na(tum),1,all),]
n<-ncol(tum)

#specify distribution of e_{I}
peL<-list()
peL[[2]]<-
        function(astar,m,J,j){
         choose(m[1], astar) * choose(J - m[1], m[2] -astar)/choose(J, m[2])
        }
if (n>2) {
for (j in 3:n)
{
    peL[[j]]<-
        function(astar,m,J,j){
         e <- (astar:(min(m[-1])))
        sum( as.numeric(lapply(e,peL[[j-1]],m[-1],J,j-1)) *choose(m[1], astar) * 
        choose(J - m[1], e - astar)/choose(J, e))
        }

}
          }
        
#specify distribution of a_{I}        
paL<-list()
paL[[2]]<-
                function(astar,m,J,j){
                 e <- astar:m[1]
                sum((0.5)^e * choose(e, astar) * choose(m[1],
                            e) * choose(J - m[1], m[2] - e)/choose(J, m[2]))
                }
if (n>2) {
for (j in 3:n)
{
    paL[[j]]<-function(astar,m,J,j){
     e <- astar:m[1]
    sum((0.5^(j-1))^astar * (1-0.5^(j-1))^(e-astar) *choose(e, astar) * as.numeric(lapply(e,peL[[j]],m,J,j)))
    }
}
          }


a<-sum(apply(tum==loh1,1,all))+sum(apply(tum==loh2,1,all))  
J<-nrow(tum)
m<-apply(tum!=noloh, 2, sum)
m<-sort(m)
    tmp <- 0 
  for (astar in a:J) {
      if (astar ==a)  tmp <- tmp+0.5*paL[[n]](astar,m,J,n)
else tmp <- tmp+paL[[n]](astar,m,J,n)
    }

    return(tmp )
}


out<-list()
cat("Testing clonality for patient ")    
for (i in 1:length(unique(ptlist)))
{w<-which(ptlist==unique(ptlist)[i])
cat (paste(unique(ptlist)[i],", ",sep=""))


adjt<-minpadjECM(LOHtable[,w],Nsim=Nsim,noloh,loh1,loh2)
out[[i]]<-adjt
}
names(out)<- unique(ptlist)

out
}

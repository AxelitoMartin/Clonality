refvector <-
function(p,cc,m,nref=1000)
{

simvecmh<-NULL
for (i in c(1:nref))
{
tum<-tumors(p,cc,m)
simvecmh<-c(simvecmh,LRts(tum,p,0.8)[1])
}
return(simvecmh)
}


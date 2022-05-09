indiv.test <-
function(ss1, ss2, func, gainthres, lossthres, Nsim=100,
                       segmethod, segpar) {
###### tests whether two chromosomes have identical changes: returns distribution of the test statistic t# under hypotheses of clonality and independence, NA if discordant changes

  s1<-ss1$output
  s2<-ss2$output

  s1c<-cumsum(s1$num.mark)
  n1<-nrow(s1)
  s2c<-cumsum(s2$num.mark)
  n2<-nrow(s2)

  xx1<-ss1$data[,3]
  xx2<-ss2$data[,3]
  resid1<-xx1
  resid1[!is.na(xx1)]<-
    resid1[!is.na(xx1)]-rep(ss1$output$seg.mean,ss1$output$num.mark)
  resid2<-xx2
  resid2[!is.na(xx2)]<-
    resid2[!is.na(xx2)]-rep(ss2$output$seg.mean,ss2$output$num.mark)
  chrlen<-length(xx1)

  ts<-func(ss1,ss2)
  if (is.na(ts[1]))
    return(NA)
  else {
    b1<-ts[2]
    b2<-ts[3]
    mn1main<-s1$seg.mean[b1]
    mn1other<-s1$seg.mean[-b1][1]
    mn2main<-s2$seg.mean[b2]
    mn2other<-s2$seg.mean[-b2][1]
    res2<-rep(NA,Nsim)
### indep
    resi<-0
    res<-rep(NA,Nsim)
    ct<-0
    while (resi<Nsim & ct<500) {
      ct<-ct+1
      if (ts[4]==1) # there is overlap
        {hotspot<-sample(max(1,s1c[b1-1],s2c[b2-1]):min(s1c[b1],s2c[b2]),1)
         breaks1<-c(sort(sample(hotspot,max(0,b1-1))),
                    sort(sample(c((hotspot+1):chrlen),n1-b1)))
         breaks2<-c(sort(sample(hotspot,max(0,b2-1))),
                    sort(sample(c((hotspot+1):chrlen),n2-b2)))
       }
      else
        {hotspot<-sample(max(1,s1c[b1-1]):min(s1c[b1]),1)
         breaks1<-c(sort(sample(hotspot,max(0,b1-1))),
           sort(sample(c((hotspot+1):chrlen),n1-b1)))
         hotspot<-sample(max(1,s2c[b2-1]):min(s2c[b2]),1)
         breaks2<-c(sort(sample(hotspot,max(0,b2-1))),
                    sort(sample(c((hotspot+1):chrlen),n2-b2)))
       }

      mns1<-rep((mn1other+mn2other)/2,n1)
      mns1[b1]<-(mn1main+mn2main)/2

      mns2<-rep((mn1other+mn2other)/2,n2)
      mns2[b2]<-(mn1main+mn2main)/2

      x1<-sample(resid1)+rep(mns1,c(breaks1,chrlen)-c(0,breaks1))
      x2<-sample(resid2)+rep(mns2,c(breaks2,chrlen)-c(0,breaks2))

      sseg1<-segment1(CNA(x1,rep(1,length(x1)),ss1$data$maploc),segmethod=segmethod,segpar=segpar)
      sseg2<-segment1(CNA(x2,rep(1,length(x2)),ss2$data$maploc),segmethod=segmethod,segpar=segpar)

      sseg1$output[,7]<-"Normal"
      sseg2$output[,7]<-"Normal"

      sseg1$output[sseg1$output$seg.mean>(gainthres[1]+gainthres[2])/2,7]<-"Gain"
      sseg1$output[sseg1$output$seg.mean<(lossthres[1]+lossthres[2])/2,7]<-"Loss"

      sseg2$output[sseg2$output$seg.mean>(gainthres[1]+gainthres[2])/2,7]<-"Gain"
      sseg2$output[sseg2$output$seg.mean<(lossthres[1]+lossthres[2])/2,7]<-"Loss"

      names(sseg1$output)[7]<-"state"
      names(sseg2$output)[7]<-"state"

      tsp<-func(sseg1,sseg2)

      if (!is.na(tsp[1])) 
        {resi<-resi+1
         res[resi]<-tsp[1]
       }
    }

    if (sum(!is.na(res))<10) return(NA)

### clonal

    res2<-rep(NA,Nsim)
    resi<-0
    ct<-0
    while (resi<Nsim & ct<500)
      {ct<-ct+1

       if (ts[4]==1) hotspot<-sample(max(1,s1c[b1-1],s2c[b2-1]):min(s1c[b1],s2c[b2]),1)
       else          hotspot<-sample(min(s1c[b1],s2c[b2]):max(s1c[b1-1],s2c[b2-1]),1)

       if (n1==n2 & b1==b2) {
         breaks1<-breaks2<-c(sort(sample(hotspot,max(0,b1-1))),
                             sort(sample(c((hotspot+1):s1c[n1]),n1-b1)))
         mns1<-mns2<-(s1$seg.mean+s2$seg.mean)/2
       }
       else if (n1<n2) 
         {breaks1<-breaks2<-c(sort(sample(hotspot,max(0,b1-1))),
                              sort(sample(c((hotspot+1):s1c[n1]),n1-b1)))
          if (b1==1) mns1<-mns2<-c((mn1main+mn2main)/2,(mn1other+mn2other)/2)
          else mns1<-mns2<-c((mn1other+mn2other)/2,(mn1main+mn2main)/2)
        }
       else if (n2<n1)
         {breaks1<-breaks2<-c(sort(sample(hotspot,max(0,b2-1))),
                              sort(sample(c((hotspot+1):s2c[n2]),n2-b2)))
          if (b2==1) mns1<-mns2<-c((mn1main+mn2main)/2,(mn1other+mn2other)/2)
          else mns1<-mns2<-c((mn1other+mn2other)/2,(mn1main+mn2main)/2)
         }

       lens1<-c(breaks1,chrlen)-c(0,breaks1)
       lens2<-c(breaks2,chrlen)-c(0,breaks2)

       x1<-sample(resid1)+rep(mns1,lens1)
       x2<-sample(resid2)+rep(mns2,lens2)

       sseg1<-segment1(CNA(x1,rep(1,length(x1)),ss1$data$maploc),segmethod=segmethod,segpar=segpar)
       sseg2<-segment1(CNA(x2,rep(1,length(x2)),ss2$data$maploc),segmethod=segmethod,segpar=segpar)

       sseg1$output[,7]<-"Normal"
       sseg2$output[,7]<-"Normal"

       sseg1$output[sseg1$output$seg.mean>(gainthres[1]+gainthres[2])/2,7]<-"Gain"
       sseg1$output[sseg1$output$seg.mean<(lossthres[1]+lossthres[2])/2,7]<-"Loss"

       sseg2$output[sseg2$output$seg.mean>(gainthres[1]+gainthres[2])/2,7]<-"Gain"
       sseg2$output[sseg2$output$seg.mean<(lossthres[1]+lossthres[2])/2,7]<-"Loss"

       names(sseg1$output)[7]<-"state"
       names(sseg2$output)[7]<-"state"

       tsp<-func(sseg1,sseg2)

       if (!is.na(tsp[1])) 
         {resi<-resi+1
          res2[resi]<-tsp[1]
        }
     }
    if (sum(!is.na(res2))<10) return(NA)

    res<-res[!is.na(res)]
    res2<-res2[!is.na(res2)]
    a<-density(res,from=0,to=length(xx1))
    p1<-a$y[sort.list(abs(a$x-ts[1]))[1]]
    a<-density(res2,from=0,to=length(xx1))
    p2<-a$y[sort.list(abs(a$x-ts[1]))[1]]
    pvalue<-mean(res<=ts[1])

    return(list(ts[1],pvalue,p1,p2))
  }
}


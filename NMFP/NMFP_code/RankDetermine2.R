RankDetermine2<-function(V,Cset,alpha0,ntest=20,nave=10,ExonPool,LenRead,Rstart,Rend,input_file){
  
  nRow<-dim(V)[1]
  nCol<-dim(V)[2]
  nExon<-length(ExonPool[1,])
  ELen<-ExonPool[2,]-ExonPool[1,]+1
  
  Wk<-NULL
  WkStar<-matrix(0,ntest,(Rend-Rstart+1))
  Gapk<-NULL
  ek<-NULL
  sdk<-NULL
  
  ReadRange<-apply(V,1,min)
  ReadRange<-cbind(ReadRange,apply(V,1,max))
  
  for(i in Rstart:Rend){
  #Real Error
  Error<-0
  for(j in 1:nave){
    ncNMFresult<-cmp_ncNMF2(V=V,alpha0=alpha0,rank=i,Cset=Cset,input_file=input_file)
    Error<-Error+ncNMFresult$err
  }
  Wk<-c(Wk,Error/nave)
  

  for(k in 1:ntest){
    Vk<-t(apply(t(t(seq(nRow))),1,function(x){return(runif(nCol,min=ReadRange[x,1],max=ReadRange[x,2]))}))
    Error<-0
    for(j in 1:nave){
      ncNMFresult<-cmp_ncNMF2(V=Vk,alpha0=alpha0,rank=i,Cset=Cset,input_file=input_file)
      Error<-Error+ncNMFresult$err
    }
    WkStar[k,(i-Rstart+1)]<-Error/nave
  }
  ek<-c(ek,1/ntest*sum(log(WkStar[,(i-Rstart+1)])))
  Gapk<-c(Gapk,ek[(i-Rstart+1)]-log(Wk[(i-Rstart+1)]))
  sdk<-c(sdk,sqrt(1+1/ntest)*sqrt(1/ntest*sum((log(WkStar[,(i-Rstart+1)])-ek[(i-Rstart+1)])^2)))
  
  if(i>Rstart&&(Gapk[i-Rstart]>=Gapk[i-Rstart+1]-sdk[i-Rstart+1])){
    break
  }
  
  }
  
  rank<-i-1
  
  return(rank)  
}
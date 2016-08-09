Normalization<-function(ReadM,ExonPool,LenRead,TypesBox){
  nSample<-length(ReadM[1,])
  nExon<-length(ExonPool[1,])
  nType<-length(TypesBox)

  #per Ave reads
  AveRead<-mean(apply(ReadM,2,function(x){return(sum(x))}))
  NormM<-apply(ReadM,2,function(x){x/sum(x)})*AveRead
  
  #per 10000 bases
  for(k in 1:nSample){
    
    NormM[,k]<-t(t(apply(t(seq(nType)),2,function(i){return(NormM[i,k]/PairReadLength(ExonPool,LenRead,TypesBox[i]))})))
    
  }

  return(NormM)
}
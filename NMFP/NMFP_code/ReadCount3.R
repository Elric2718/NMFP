ReadCount3<-function(ReadM,ExonBD,TypesBox,LenRead){
  nExon<-length(ExonBD[1,])
  nRead<-length(ReadM[,1])
  ExonLen<-ExonBD[2,]-ExonBD[1,]+1
  SIndex<-which(ExonLen<LenRead)
  SIndex<-setdiff(SIndex,c(1,nExon))
  SIndex<-c(SIndex,0)  #guarantee that SIndex is a vector with length over 1
  
  IndexList1<-vector(mode="list")
  for(i in 1:nExon){
    IndexList1[[i]]<-which(apply(t(ReadM[,1]),2,function(x){
      if(x>=ExonBD[1,i]-1&x<=ExonBD[2,i]+1){return(TRUE)}else{return(FALSE)}})==TRUE)
  }
  IndexList2<-vector(mode="list")
  for(i in 1:nExon){
    IndexList2[[i]]<-which(apply(t(ReadM[,2]),2,function(x){
      if(x>=ExonBD[1,i]-1&x<=ExonBD[2,i]+1){return(TRUE)}else{return(FALSE)}})==TRUE)
  }

  
  TypeNum<-rep(0,length(TypesBox))
  for(i in 1:length(TypesBox)){
    NoSite1<-as.numeric(as.character(strsplit(TypesBox[i],":")[[1]][1]))
    NoSite2<-as.numeric(as.character(strsplit(TypesBox[i],":")[[1]][2]))
    
    ReadList<-intersect(IndexList1[[NoSite1]],IndexList2[[NoSite2]])
    
    
    if(length(ReadList)>0){
      Special0<-intersect(which(NoSite1<SIndex),which(NoSite2>SIndex))
      if(length(Special0)==0){
        TypeNum[i]<-TypeNum[i]+length(ReadList)
      }
    
      
    if(length(Special0)>0){
      for(k in 1:length(Special0)){
        Special<-Special0[k]
        Special<-SIndex[Special]
        SkipList<-ReadList[which(apply(t(ReadList),2,function(x){return(
          ExonBD[2,NoSite1]-ReadM[x,1]+1+ReadM[x,2]-ExonBD[1,NoSite2]+1)})>=LenRead)]
        RetentionList<-setdiff(ReadList,SkipList)
        
        TypeNum[i]<-TypeNum[i]+length(SkipList)
        TypeNum[which(TypesBox==paste(NoSite1,Special,sep=":"))]<-TypeNum[which(TypesBox==paste(NoSite1,Special,sep=":"))]+length(RetentionList)
        TypeNum[which(TypesBox==paste(Special,Special,sep=":"))]<-TypeNum[which(TypesBox==paste(Special,Special,sep=":"))]+length(RetentionList)
        TypeNum[which(TypesBox==paste(Special,NoSite2,sep=":"))]<-TypeNum[which(TypesBox==paste(Special,NoSite2,sep=":"))]+length(RetentionList)
      }
    }
    }
  }

  
  
  
  return(TypeNum)
  
}
VisualizeG5<-function(GSignal,TypesBox,ExonPool,LenRead,IfPrint=FALSE,IfTune=FALSE){
  nExon<-length(ExonPool[1,])
  nIsoform<-length(GSignal[1,])
  nType<-length(TypesBox)
  
  LenExon<-ExonPool[2,]-ExonPool[1,]+1
  Sindex<-which(LenExon<LenRead)
  Sindex<-c(Sindex,0)#guarantee it contains at least one element
  
  #indicator vector
  #weights indicator vector
  Ip1<-rep(0,nExon)
  In1<-rep(0,nExon)
  #non-weights indicator vector
  Ip2<-rep(0,nExon)
  In2<-rep(0,nExon)
  
  for(i in 1:nExon){
    if(i>=4&&i<=nExon-3){
      Ip1[i]<-8
      In1[i]<--8
      
      Ip2[i]<-7
      In2[i]<--4
    }
    
    if(i==3|i==nExon-2){
      Ip1[i]<-7
      In1[i]<--8
      
      Ip2[i]<-6
      In2[i]<--4
    }
    
    if(i==2|i==nExon-1){
      Ip1[i]<-6
      In1[i]<--6
      
      Ip2[i]<-5
      In2[i]<--3
    }
    
    if(i==1|i==nExon){
      if(is.element(i,Sindex)==FALSE){
        Ip1[i]<-5
        In1[i]<--2.5
      
        Ip2[i]<-4
        In2[i]<--1
      }else{
        Ip1[i]<-3
        In1[i]<--5/3
        
        Ip2[i]<-3
        In2[i]<--3
      }
    }
  }
  
  
  
  
  #weights matrix
  Vp1<-matrix(0,nExon,nIsoform)
  Vn1<-matrix(0,nExon,nIsoform)
  
  #non-weights matrix
  Vp2<-matrix(0,nExon,nIsoform)
  Vn2<-matrix(0,nExon,nIsoform)
  for(i in 1:nType){
    for(j in 1:nIsoform){
      Vp1[,j]<-Vp1[,j]+VisualCheck4(ExonPool,LenRead,TypesBox[i],GSignal[i,j])$Positive1
      Vn1[,j]<-Vn1[,j]+VisualCheck4(ExonPool,LenRead,TypesBox[i],GSignal[i,j])$Negative1
      
      Vp2[,j]<-Vp2[,j]+VisualCheck4(ExonPool,LenRead,TypesBox[i],GSignal[i,j])$Positive2
      Vn2[,j]<-Vn2[,j]+VisualCheck4(ExonPool,LenRead,TypesBox[i],GSignal[i,j])$Negative2
    }
  }
  

  
  #trimmed Gmatrix
  TrimG<-matrix(0,nExon,nIsoform)
  for(i in 1:nIsoform){
    for(j in 1:nExon){
      if(Vn2[j,i]<0&&Vp2[j,i]==0){TrimG[j,i]<--1}
      if(Vn2[j,i]==0&&Vp2[j,i]>0){TrimG[j,i]<-1}
      if(Vn2[j,i]<0&&Vp2[j,i]>0){TrimG[j,i]<-0}
      
          
     }
  }
      
  
  Uncertainty<-NULL
  #separate uncertain isoforms
  CertainV<-NULL
  for(i in 1:nIsoform){
    Uindex<-which(TrimG[,i]==0)
    Uncertainty<-c(Uncertainty,length(Uindex))
    if(IfTune==FALSE){
      if(length(Uindex)==0){
        CertainV<-cbind(CertainV,TrimG[,i])}else{
        #v<-TrimG[,i]
        #v[Uindex]<-(-1)
        #CertainV<-cbind(CertainV,t(t(v)))
         # print(length(Uindex))
          v<-matrix(TrimG[,i],nrow=nExon,ncol=2^(length(Uindex)))
          #Findex<-max(which(as.numeric(intToBits(j-1))==1))
          v[Uindex,]<-matrix(as.numeric(intToBits(seq(0,2^(length(Uindex))-1))),nrow=32)[1:length(Uindex),]
          v[which(v==0)]<-(-1)
          CertainV<-cbind(CertainV,t(t(v)))  
      }
    }
  }
  
  
  
  if(IfPrint==TRUE){
    print(Vp1)
    print(Vn1)
    print(Vp2)
    print(Vn2)
    print(Ip1)
    print(In1)
    print(Ip2)
    print(In2)
  }
  
  
  result<-list(CertainV,mean(Uncertainty))
  names(result)<-c("CertainV","Uncertainty")
  return(result)
}
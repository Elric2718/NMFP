PreSelection<-function(Candidates=NULL,SumRead,ExonPool,LenRead,ifPart=FALSE,mode=4){

  nExon<-length(ExonPool[1,])
  LenExon<-ExonPool[2,]-ExonPool[1,]+1
  Sindex<-union(intersect(which(LenExon<=LenRead+10),which(LenExon>=LenRead-10)),which(LenExon<=20))
  Ref<-SumRead
  Ref[which(Ref>0)]<-1
  nType<-length(SumRead)
  CANDIDATES<-Candidates
  
  if(ifPart==FALSE){
    if(length(Candidates)==0){
    Candidates<-matrix(as.numeric(intToBits(seq(1,2^(nExon)-1))),nrow=32)[1:nExon,] 
    #Candidates<-as.vector(rbind(Candidates,":"))  
    #Candidates<-strsplit(paste(Candidates,collapse=""),":")[[1]]
    }else{
      Candidates<-setdiff(Candidates,paste(rep(0,nExon),collapse=""))
      Candidates<-paste(Candidates,collapse="")
      Candidates<-as.numeric(strsplit(Candidates,"")[[1]])
      Candidates<-matrix(Candidates,nrow=nExon)
    }
    nCand<-length(Candidates[1,])
  

    IndexG<-which(apply(Candidates,2,function(x){
     #stage 1
      Match1<-prod(Ref[setdiff(which(x==1),c(Sindex,1,nExon))])
      if(LenExon[1]>=LenRead+5&&x[1]==1){
        Match1<-Match1*Ref[1]
      }
      if(LenExon[nExon]>=LenRead+5&&x[nExon]==1){
        Match1<-Match1*Ref[nExon]
      }
    
      #stage 2
       Match2<-0
      if(Match1==1){
        Match2<-prod(Ref[nExon+which(
        apply(t(seq(nExon-1)),2,function(i){if(x[i]*x[i+1]==1){return(TRUE)}else{return(FALSE)}})==TRUE)])
      }
    
      #stage 3
      Match3<-0
      if(Match2==1){
        Match3<-prod(Ref[2*nExon-1+which(
          apply(t(seq(nExon-2)),2,function(i){if(x[i]*x[i+2]==1&&x[i+1]==0){return(TRUE)}else{return(FALSE)}})==TRUE)])
      }
    
      #stage 4
      if(nExon>=4){
      Match4<-0
        if(Match3==1){
        Match4<-prod(Ref[3*nExon-3+which(
          apply(t(seq(nExon-3)),2,function(i){if(x[i]*x[i+3]==1&&x[i+1]==0&&x[i+2]==0){return(TRUE)}else{return(FALSE)}})==TRUE)])
      }
      }else{Match4<-1}
    
      
      
      if(switch(mode,Match1,Match2,Match3,Match4)==1){
        return(TRUE)
      }else{return(FALSE)}
    
      })==TRUE)
  

      if(length(IndexG)>0){
      SlctCandidates<-Candidates[,IndexG]
        if(length(IndexG)==1){
          SlctCandidates<-t(t(SlctCandidates))
       }
       SlctCandidates<-as.vector(rbind(SlctCandidates,":"))  
       SlctCandidates<-strsplit(paste(SlctCandidates,collapse=""),":")[[1]]
      }else{
        SlctCandidates<-setdiff(CANDIDATES,paste(rep(0,nExon),collapse=""))
      }
    
      }else{
        SlctCandidates<-setdiff(CANDIDATES,paste(rep(0,nExon),collapse=""))
      }
      
  return(SlctCandidates)
}
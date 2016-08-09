###Get the boundary of exons from annotation
ExonBoundary<-function(genelist,Anno,Gene="G8944",LenRead,ifvis=F,Emode=2,Anno_mode=1){


Trans<-genelist[which(as.character(genelist[,1])==Gene),]
nTrans<-length(as.character(Trans[,1]))
chr<-as.character(Trans[1,2])
start<-min(as.numeric(as.character(Trans[,3])))
end<-max(as.numeric(as.character(Trans[,4])))
strand<-as.character(Trans[1,5])
len<-end-start+1

if(Anno_mode==1){
  anno<-Anno[,c(1,3,4,5,7,16,10)]
}else{
  anno<-Anno[,c(1,3,4,5,7,13,10)]
}




TransPool<-union(NULL,as.character(Trans[,6]))

IndexExon<-intersect(which(as.character(anno[,2])=="exon"),which(is.element(as.character(anno[,6]),TransPool)))

ExonStart0<-as.numeric(as.character(anno[IndexExon,3]))-start+1
ExonEnd0<-as.numeric(as.character(anno[IndexExon,4]))-start+1
ExonPool0<-rbind(ExonStart0,ExonEnd0)
ExonPool0<-ExonPool0[,which(duplicated(as.data.frame(t(ExonPool0)))==F)]
if(length(ExonPool0)==2){
  ExonPool0<-t(t(ExonPool0))
}else{
ExonPool0<-ExonPool0[,order(ExonPool0[1,])]
}


Node<-union(ExonStart0,ExonEnd0)
Node<-Node[order(Node)]
ExonPool<-NULL
for(i in 1:(length(Node)-1)){
  if(length(intersect(which(Node[i]>=ExonPool0[1,]),which(Node[i+1]<=ExonPool0[2,])))>0){
    ExonPool<-cbind(ExonPool,c(Node[i],Node[i+1]))
  }
}

ExonPool<-as.matrix(ExonPool,nrow=2)
#for(i in 1:(length(ExonPool[1,])-1)){
#  if(ExonPool[2,i]==ExonPool[1,(i+1)]){
#    ExonPool[1,(i+1)]<-ExonPool[1,(i+1)]+1
#  }
#}

ExonPool1<-ExonPool

nExon<-length(ExonPool[1,])
LenExon<-ExonPool[2,]-ExonPool[1,]+1
if(Emode>0&&nExon>1){

FirstTime=TRUE
Sindex=NULL
while(length(Sindex)>0|FirstTime){
Sindex<-intersect(which(LenExon<LenRead-10),which(apply(t(seq(nExon)),2,function(x){
                          if(x==1&&ExonPool[2,1]>=ExonPool[1,2]-1){return(TRUE)}
                          if(x==nExon&&ExonPool[1,nExon]-1<=ExonPool[2,(nExon-1)]){return(TRUE)}
                          if(x>1&&x<nExon&&(ExonPool[2,x]>=ExonPool[1,(x+1)]-1|ExonPool[1,x]-1<=ExonPool[2,(x-1)])){return(TRUE)}
                          return(FALSE)  
                          })==TRUE))

Sindex<-Sindex[order(Sindex)]
Sindex0<-Sindex

if(Emode==2){
  if(length(Sindex)>1){
    ContIndex<-Sindex[which(Sindex[2:length(Sindex)]-Sindex[1:(length(Sindex)-1)]==1)]
    if(length(ContIndex)>0){
    Sindex<-apply(t(ContIndex),2,FUN=function(x){
      if(LenExon[x]>LenExon[x+1]){return(x+1)}else{return(x)}
    })}}else{
      Sidex<-NULL
    }
  
    Sindex<-c(Sindex,Sindex0[which(LenExon[Sindex0]<10)])
    Sindex<-unique(Sindex)
}

if(length(Sindex)==0){
  break
}
    if(Sindex[1]==1){
        ExonPool[2,1]=ExonPool[2,2]
        ExonPool[1,2]=ExonPool[1,1]
    }

    if(Sindex[1]==nExon){
      ExonPool[2,(nExon-1)]=ExonPool[2,nExon]
      ExonPool[1,nExon]=ExonPool[1,(nExon-1)]  
    }

    if(Sindex[1]<nExon&&Sindex[1]>1){
      if((ExonPool[2,Sindex[1]]>=ExonPool[1,(Sindex[1]+1)]-1)&&(ExonPool[1,Sindex[1]]-1>ExonPool[2,(Sindex[1]-1)])){
        ExonPool[2,Sindex[1]]=ExonPool[2,(Sindex[1]+1)]
        ExonPool[1,(Sindex[1]+1)]=ExonPool[1,Sindex[1]]  
      }
      
      if((ExonPool[2,Sindex[1]]<ExonPool[1,(Sindex[1]+1)]-1)&&(ExonPool[1,Sindex[1]]-1<=ExonPool[2,(Sindex[1]-1)])){
        ExonPool[2,(Sindex[1]-1)]=ExonPool[2,Sindex[1]]
        ExonPool[1,Sindex[1]]=ExonPool[1,(Sindex[1]-1)]  
      }
      
      if((ExonPool[2,Sindex[1]]>=ExonPool[1,(Sindex[1]+1)]-1)&&(ExonPool[1,Sindex[1]]-1<=ExonPool[2,(Sindex[1]-1)])){
        if(LenExon[(Sindex[1]-1)]<=LenExon[(Sindex[1]+1)]){
          ExonPool[2,(Sindex[1]-1)]=ExonPool[2,Sindex[1]]
          ExonPool[1,Sindex[1]]=ExonPool[1,(Sindex[1]-1)] 
        }else{
          ExonPool[2,Sindex[1]]=ExonPool[2,(Sindex[1]+1)]
          ExonPool[1,(Sindex[1]+1)]=ExonPool[1,Sindex[1]]  
        }
      }
      
    }


  ExonPool<-ExonPool[,which(duplicated(as.data.frame(t(ExonPool)))==F)]
  ExonPool<-as.matrix(ExonPool,nrow=2)
  nExon<-length(ExonPool[1,])
  LenExon<-ExonPool[2,]-ExonPool[1,]+1
  Sindex<-setdiff(Sindex,Sindex[1])
  FirstTime=FALSE

  #print(ExonPool)
  }

}


  




#check if there are subexons
#SubIndex<-which(ExonPool[1,(2:length(ExonPool[1,]))]-ExonPool[2,(1:(length(ExonPool[1,])-1))]<0)
#while(length(SubIndex)>0){
#  node<-NULL
#  if(length(SubIndex)>1){
#    node<-which(SubIndex[2:length(SubIndex)]-SubIndex[1:(length(SubIndex)-1)]>1)
#  }
  
#  node<-union(node,length(SubIndex))

#  nodeStart<-1
#  for(i in 1:length(node)){
#    nodeEnd<-node[i]
#    Inter<-NULL
#    for(j in nodeStart:nodeEnd){
#      Inter<-union(c(ExonPool[1,SubIndex[j]],ExonPool[2,SubIndex[j]],ExonPool[1,SubIndex[j]+1]),ExonPool[2,SubIndex[j]+1]) 
#    }
#    Inter<-Inter[order(Inter)]
#    for(j in 1:(length(Inter)-1)){
#      ExonPool<-cbind(ExonPool,c(Inter[j],Inter[j+1]))
#    }
#    nodeStart<-nodeEnd+1
#  }
  
#  ExonPool<-ExonPool[,-union(SubIndex,SubIndex+1)]  
  
#  ExonPool<-ExonPool[,which(duplicated(as.data.frame(t(ExonPool)))==F)]
#  if(length(ExonPool)==2){
#    ExonPool<-t(t(ExonPool))
#  }else{
#    ExonPool<-ExonPool[,order(ExonPool[1,])]
#  }

#  SubIndex<-which(ExonPool[1,(2:length(ExonPool[1,]))]-ExonPool[2,(1:(length(ExonPool[1,])-1))]<0)
#}

##visualization
if(ifvis){
 exonpool<-union(NULL,c(unlist(apply(ExonPool,2,function(x){return(seq(x[1],x[2]))}))))
 plot(x=exonpool,y=rep(0,length(exonpool)),xlim=c(0,ExonPool[2,nExon]+100),ylim=c(-1,nTrans+1),,xlab="bases",ylab="transcript",pch=".")
 abline(v=union(ExonPool[1,],ExonPool[2,]),col="red")
 for(i in 1:nTrans){
   par(new=TRUE)
   indexexon<-intersect(which(as.character(anno[,2])=="exon"),which(as.character(anno[,6])==TransPool[i]))
   exon<-rbind(as.numeric(as.character(anno[indexexon,3]))-start+1,as.numeric(as.character(anno[indexexon,4]))-start+1)
   trans<-union(NULL,c(unlist(apply(exon,2,function(x){return(seq(x[1],x[2]))}))))
   plot(x=trans,y=rep(i,length(trans)),xlim=c(0,ExonPool[2,nExon]+100),ylim=c(-1,nTrans+1),xlab="bases",ylab="transcript",pch=".")
 }
 par(new=FALSE)
}



#get annotated transcripts
annoTranscript<-matrix(0,nExon,nTrans)
annoTrans<-NULL
for(i in 1:nTrans){
  indexexon<-intersect(which(as.character(anno[,2])=="exon"),which(as.character(anno[,6])==TransPool[i]))
  exonstart<-as.numeric(as.character(anno[indexexon,3]))-start+1
  exonend<-as.numeric(as.character(anno[indexexon,4]))-start+1
  for(j in 1:nExon){
    if(Emode==0){
      if(length(intersect(which(exonstart<=ExonPool[1,j]),which(exonend>=ExonPool[2,j])))>0){
        annoTranscript[j,i]<-1
      }
    }else{
      if(length(intersect(which(exonstart<=ExonPool[1,j]),which(exonend>=ExonPool[1,j]+1)))>0|length(intersect(which(exonstart<=ExonPool[2,j]-1),which(exonend>=ExonPool[2,j])))>0){
        annoTranscript[j,i]<-1
      }
    }
  }
  annoTrans<-c(annoTrans,paste(annoTranscript[,i],collapse=""))
}
annoTrans<-unique(annoTrans)




ExonList<-list(chr,strand,start,end,ExonPool,annoTrans)
names(ExonList)<-c("chr","strand","GeneStart","GeneEnd","ExonPool","annoTrans")
return(ExonList)
}
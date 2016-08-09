VoteforIsoform<-function(rank=0,nrun=10,alpha0=1,NormMatrix,TypesBox,ExonPool,LenRead=76,Gcutoff=0.4,mode=1,input_file){
  
  nExon<-length(ExonPool[1,])
  ###apply NMF
  IndexNoneZero<-which(apply(NormMatrix,1,mean)>0*mean(apply(NormMatrix,1,mean)))
  #ContradictTable
  Ctable<-ContradictTable(TypesBox,nExon)
  #aheatmap(Ctable,scale="none",Rowv=NA,Colv=NA,revC=F,col=c("blue","red"))
  Cset<-Ctable[IndexNoneZero,IndexNoneZero]
  #aheatmap(Cset,scale="none",Rowv=NA,Colv=NA,revC=F,col=c("blue","red"))

  #cmp_RankDetermine<-cmpfun(RankDetermine)
  #determine ranks
  if(rank==0){
  rank<-RankDetermine2(V=NormMatrix[IndexNoneZero,],Cset=Cset,
                               alpha0=1,ExonPool=ExonPool,LenRead=LenRead,Rstart=2,Rend=5,input_file=input_file)

  }
  
  #tuning parameters
  if(mode==0){
  alpha0<-c(1,0.1,10)
  rank<-c(rank,rank+1,rank-1)
  Gcutoff<-c(0.4,0.1,0.01)
  
  TuneCertainty<-NULL
  for(i in 1:length(alpha0)){
    for(j in 1:length(rank)){
      for(k in 1:length(Gcutoff)){
        TuneCertainty<-c(TuneCertainty,max(apply(t(seq(10)),2,function(x){ncNMFresult<-cmp_ncNMF2(V=NormMatrix[IndexNoneZero,],alpha0=alpha0[i],rank=rank[j],Cset=Cset,input_file=input_file)
        Basis<-ncNMFresult$basis
        Coef<-ncNMFresult$coefficient
        #stopCluster(cl)
        #rm(cl)
        #aheatmap(Basis,scale="none",Rowv=NA,Colv=NA,revC=F,col=rev(heat.colors(256)))
        
        ##Affinity Matrix
        AffinityMatrix<-diag(apply(Basis,1,max))
        
        
        ##G Matrix
        G<-solve(AffinityMatrix)%*%Basis
        
        #trim G
        G<-apply(G,2,function(x){y<-rep(0,length(x));y[which(x>=Gcutoff[k])]<-1;return(y)})
        #aheatmap(G,scale="none",Rowv=NA,Colv=NA,revC=F,col=rev(heat.colors(256)))
        GMatrix<-matrix(0,dim(NormMatrix)[1],rank)
        GMatrix<-apply(G,2,function(x){y<-rep(0,dim(GMatrix)[1]);y[IndexNoneZero[which(x==1)]]<-1;return(y)})
        
        
        #visualize GMatrix
        GVisual<-VisualizeG5(GMatrix,TypesBox,ExonPool,LenRead,IfPrint=FALSE,IfTune=TRUE)
        return(GVisual$Uncertainty)
        })))
       #if(TuneCertainty[length(TuneCertainty)]<=5){break} 
      }
      #if(TuneCertainty[length(TuneCertainty)]<=5){break} 
    }
    #if(TuneCertainty[length(TuneCertainty)]<=5){break} 
  }
  
  
  if(length(TuneCertainty)<length(alpha0)*length(rank)*length(Gcutoff)){
    alpha0<-alpha0[i]
    rank<-rank[j]
    Gcutoff<-Gcutoff[k]
  }else{
    Tindex<-which.min(TuneCertainty)
    TuneCertainty<-TuneCertainty[Tindex]
    i<-ceiling(TuneCertainty/(length(rank)*length(Gcutoff)))
    j<-ceiling((TuneCertainty%%(length(rank)*length(Gcutoff)))/length(Gcutoff))
    k<-(TuneCertainty%%(length(rank)*length(Gcutoff)))%%length(Gcutoff)+1
    alpha0<-alpha0[i]
    rank<-rank[j]
    Gcutoff<-Gcutoff[k]
  }
  }else{
    alpha0=1
    Gcutoff<-c(0.4,0.1,0.01)
    
    TuneCertainty<-NULL
    for(k in 1:length(Gcutoff)){
      TuneCertainty<-c(TuneCertainty,max(apply(t(seq(10)),2,function(x){
      ncNMFresult<-cmp_ncNMF2(V=NormMatrix[IndexNoneZero,],alpha0=alpha0,rank=rank,Cset=Cset,input_file=input_file)

      Basis<-ncNMFresult$basis
      Coef<-ncNMFresult$coefficient
      #stopCluster(cl)
      #rm(cl)
      #aheatmap(Basis,scale="none",Rowv=NA,Colv=NA,revC=F,col=rev(heat.colors(256)))
      
      ##Affinity Matrix
      AffinityMatrix<-diag(apply(Basis,1,max))
      
      
      ##G Matrix
      G<-solve(AffinityMatrix)%*%Basis
      
      #trim G
      G<-apply(G,2,function(x){y<-rep(0,length(x));y[which(x>=Gcutoff[k])]<-1;return(y)})
      #aheatmap(G,scale="none",Rowv=NA,Colv=NA,revC=F,col=rev(heat.colors(256)))
      GMatrix<-matrix(0,dim(NormMatrix)[1],rank)
      GMatrix<-apply(G,2,function(x){y<-rep(0,dim(GMatrix)[1]);y[IndexNoneZero[which(x==1)]]<-1;return(y)})
      
      
      #visualize GMatrix
      GVisual<-VisualizeG5(GMatrix,TypesBox,ExonPool,LenRead,IfPrint=FALSE,IfTune=TRUE)
      return(GVisual$Uncertainty)
    }))) 
    }
    
    Tindex<-which.min(TuneCertainty)
    k<-Tindex
    Gcutoff<-Gcutoff[k]
      
    }
  #print(rank)
  #print(alpha0)
  #print(Gcutoff)
  
  #multi runs
   Bisoform<-NULL
    #Uncertainty<-0
  for(n in 1:nrun){
    #print(paste("hello",n,sep=""))
    ncNMFresult<-cmp_ncNMF2(V=NormMatrix[IndexNoneZero,],alpha0=alpha0,rank=rank,Cset=Cset,input_file=input_file)
    Basis<-ncNMFresult$basis
    Coef<-ncNMFresult$coefficient
    #stopCluster(cl)
    #rm(cl)
    #aheatmap(Basis,scale="none",Rowv=NA,Colv=NA,revC=F,col=rev(heat.colors(256)))

    ##Affinity Matrix
    AffinityMatrix<-diag(apply(Basis,1,max))

      
    ##G Matrix
    G<-solve(AffinityMatrix)%*%Basis
    
    #trim G
    G<-apply(G,2,function(x){y<-rep(0,length(x));y[which(x>=Gcutoff)]<-1;return(y)})
    #aheatmap(G,scale="none",Rowv=NA,Colv=NA,revC=F,col=rev(heat.colors(256)))
    GMatrix<-matrix(0,dim(NormMatrix)[1],rank)
    GMatrix<-apply(G,2,function(x){y<-rep(0,dim(GMatrix)[1]);y[IndexNoneZero[which(x==1)]]<-1;return(y)})


    #visualize GMatrix
    GVisual<-VisualizeG5(GMatrix,TypesBox,ExonPool,LenRead,IfPrint=FALSE)
    #Uncertainty<-max(Uncertainty,GVisual$Uncertainty)
    Gvisual<-GVisual$CertainV
    Gvisual[which(Gvisual==-1)]<-0
    Gvisual[which(Gvisual>0)]<-1
    Gvisual<-rbind(Gvisual,":")
    Gvisual<-as.vector(Gvisual)
    #Gvisual<-apply(GVisual,2,function(x){y<-rep(0,length(x));y[which(x>0)]<-1;y[which(x==-1)]<-0;return(y)})
    #represented as binary and then transform into decimal
    Bisoform<-c(Bisoform,strsplit(paste(Gvisual,collapse=""),":")[[1]])
    #Bisoform<-c(Bisoform,apply(t(seq(length(Gvisual[1,]))),2,function(x){return(paste(Gvisual[,x],collapse=""))}))
  }

#choose the top isoforms
IsoformSet<-sort(table(Bisoform),decreasing=TRUE)
#Frequency<-IsoformSet[rank]
#Isoform<-names(which(IsoformSet>=Frequency))
  
Isoform<-list(IsoformSet,GMatrix,alpha0,rank,Gcutoff)
names(Isoform)<-c("IsoformSet","GMatrix","alpha0","rank","Gcutoff")
return(Isoform)
}
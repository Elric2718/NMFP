ncNMF2<-function(V,rank,alpha0=10,Cset,iter=300,eps=0.0005,input_file){
  ##load C based function
  #load Update
  dyn.load(paste(input_file,"/Update.so",sep=""))
  dyn.load(paste(input_file,"/CpenaltyObj.so",sep=""))
  
  
  
  ##initial the matrix
  #randomly set intial matrix
  nRow<-dim(V)[1]
  nCol<-dim(V)[2]
  meanV<-sum(V)/(nRow*nCol)
  #  W<-abs(matrix(rnorm(nRow*rank),nrow=nRow,ncol=rank))
  W<-abs(matrix(rnorm(nRow*rank,mean=1/nRow,sd=sqrt(1/nRow)),nrow=nRow,ncol=rank))
  W<-apply(W,2,function(x){return(x/sum(x))})
  H<-abs(matrix(rnorm(rank*nCol,mean=meanV*nRow,sd=sqrt(meanV*nRow)),nrow=rank,ncol=nCol))
  #  H<-abs(matrix(rnorm(rank*nCol,mean=1000,sd=10),nrow=rank,ncol=nCol))
  
  ##name W and H
  rownames(W)<-1:nRow
  colnames(W)<-LETTERS[1:rank]
  
  rownames(H)<-LETTERS[1:rank]
  colnames(H)<-1:nCol
  
  #Code dependent alpha
  alpha<-mean(W%*%H)/mean(W%*%t(W))*alpha0
  
  #initial error
  err0<-.Call("Cpenalty",W,H,V,Cset,alpha)
  err<-err0/2
  count<-0
  # print(c(err,alpha))
  
  ##export values to nodes
  #clusterExport(cl,c("V","nRow","nCol","Cset","alpha"))
  
  #print("Factorization Start")
  ##update Matrix until convergence


  while(abs(err-err0)/err0>=eps&count<iter){
    alpha<-mean(W%*%H)/mean(W%*%t(W))*alpha0
    
    #err<err0&
    count<-count+1
    
    WH<-.Call("Update",W,H,V,Cset,alpha)
    W<-WH[1:nRow,]
    H<-t(WH[(nRow+1):(nRow+nCol),])
    
    W<-t(t(W))
    rownames(W)<-1:nRow
    colnames(W)<-LETTERS[1:rank]
    
    rownames(H)<-LETTERS[1:rank]
    colnames(H)<-1:nCol
    
    if(count%%10==0){
      err0<-err
      err<-.Call("Cpenalty",W,H,V,Cset,alpha)
      #print(c(count,err,alpha))
    }
    
  }
 
  #H<-diag(apply(W,2,sum))%*%H
  #W<-W%*%diag(1/apply(W,2,sum))
  
  
  #print("Factorization ends")
  
  #unload Update  
  dyn.unload(paste(input_file,"/Update.so",sep=""))  
  dyn.unload(paste(input_file,"/CpenaltyObj.so",sep="")) 
  
  ncNMFresult<-list(W,H,err)
  names(ncNMFresult)<-c("basis","coefficient","error")
  return(ncNMFresult)
  
}
ContradictTable<-function(TypesBox,nExon){
  nTypes<-length(TypesBox)
  Ctable<-matrix(0,nTypes,nTypes)
  for(i in 1:nTypes){
    Vtest<-rep(0,nExon)
    Vtest<-VisualCheck(Vtest,TypesBox[i],1)
    IndexCdict<-apply(t(seq(nTypes)),2,function(x){y<-VisualCheck(Vtest,TypesBox[x],1);if(min(y)<(-1)){return(1)}else{return(0)}})
    Ctable[which(IndexCdict==1),i]<-1
  }
  return(Ctable)
}
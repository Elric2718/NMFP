FindTypes2<-function(nExon){
  TypesBox<-NULL
  for(i in 1:nExon){
    TypesBox<-c(TypesBox,paste(c(i,i),collapse=":"))
  }
  
  for(i in 1:(nExon-1)){
    TypesBox<-c(TypesBox,paste(c(i,i+1),collapse=":"))
  }
  
  if(nExon>=3){
  for(i in 1:(nExon-2)){
    TypesBox<-c(TypesBox,paste(c(i,i+2),collapse=":"))
  }
  }
  
  if(nExon>=4){
  for(i in 1:(nExon-3)){
    TypesBox<-c(TypesBox,paste(c(i,i+3),collapse=":"))
  }
  }
  return(TypesBox)
}
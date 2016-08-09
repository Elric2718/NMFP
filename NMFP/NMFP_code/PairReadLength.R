PairReadLength<-function(ExonBD,LenRead,Type){
  LenExon<-ExonBD[2,]-ExonBD[1,]+1
  
     
  if(length(strsplit(Type,":")[[1]])==2){
    NoSite1<-as.numeric(as.character(strsplit(Type,":")[[1]][1]))
    NoSite2<-as.numeric(as.character(strsplit(Type,":")[[1]][2]))
    
    if(NoSite1==NoSite2){
      if(LenRead-1<=LenExon[NoSite1]){
        len<-LenExon[NoSite1]-LenRead+1}else{
          len<-LenRead-LenExon[NoSite1]}
    }else{
      if(LenRead-1<=min(LenExon[NoSite1],LenExon[NoSite2])){
      len<-LenRead-1
      }else{
      len<-min(LenExon[NoSite1],LenExon[NoSite2])
      }
      }
    
    len<-max(len,1)/100
    
  }
  
  
  if(length(strsplit(Type,":")[[1]])==4){
    NoSite1<-as.numeric(as.character(strsplit(Type,":")[[1]][1]))
    NoSite2<-as.numeric(as.character(strsplit(Type,":")[[1]][2]))
    NoSite3<-as.numeric(as.character(strsplit(Type,":")[[1]][3]))
    NoSite4<-as.numeric(as.character(strsplit(Type,":")[[1]][4]))

    
    
    if(NoSite1==NoSite2){
      if(LenRead-1<=LenExon[NoSite1]){
        len1<-(LenExon[NoSite1]-LenRead+1)}else{
          len1<-(LenRead-LenExon[NoSite1])}
    }
    
    if(NoSite1<NoSite2){
      if(LenRead-1<=min(LenExon[NoSite1],LenExon[NoSite2])){
          len1<-LenRead-1
      }else{
        len1<-min(LenExon[NoSite1],LenExon[NoSite2])
      }      
    }
    
    
    if(NoSite3==NoSite4){
      if(LenRead-1<=LenExon[NoSite3]){
        len2<-(LenExon[NoSite3]-LenRead+1)}else{
          len2<-(LenRead-LenExon[NoSite3])}
    }
    
    if(NoSite3<NoSite4){
      if(LenRead-1<=min(LenExon[NoSite3],LenExon[NoSite4])){
        len2<-LenRead-1
      }else{
        len2<-min(LenExon[NoSite3],LenExon[NoSite4])
      }      
    }
    
    if(NoSite1==NoSite3&&NoSite2==NoSite4){
      len<-max(len1,1)*(max(len2,1)+1)/2/(100*100)
    }else{
      len<-max(len1,1)*max(len2,1)/(100*100)
    }
    
  }
  
  if(length(strsplit(Type,":")[[1]])==6){
    NoSite1<-as.numeric(as.character(strsplit(Type,":")[[1]][1]))
    NoSite2<-as.numeric(as.character(strsplit(Type,":")[[1]][2]))
    NoSite3<-as.numeric(as.character(strsplit(Type,":")[[1]][3]))
    NoSite4<-as.numeric(as.character(strsplit(Type,":")[[1]][4]))
    NoSite5<-as.numeric(as.character(strsplit(Type,":")[[1]][5]))
    NoSite6<-as.numeric(as.character(strsplit(Type,":")[[1]][6]))
    
    
  if(NoSite1==NoSite3){
      if(LenRead-1<=LenExon[NoSite1]){
        len1<-(LenExon[NoSite1]-LenRead+1)*(LenExon[NoSite1]-LenRead+1+1)/2}else{
          len1<-(LenRead-LenExon[NoSite1])*(LenRead-LenExon[NoSite1]+1)/2}
    }
  
  if(NoSite1<NoSite3){
    if(LenRead-1<=min(LenExon[NoSite1],LenExon[NoSite3])){
      if(NoSite2==NoSite1){
        len1<-LenRead-1
      }else{
        len1<-LenRead-LenExon[NoSite2]-1
      }     
    }else{
      len1<-min(LenExon[NoSite1],LenExon[NoSite3])
    }
    
  }
  
  
  if(NoSite4==NoSite6){
    if(LenRead-1<=LenExon[NoSite4]){
      len2<-(LenExon[NoSite4]-LenRead+1)*(LenExon[NoSite4]-LenRead+1+1)/2}else{
        len2<-(LenRead-LenExon[NoSite4])*(LenRead-LenExon[NoSite4]+1)/2}
  }
  
  if(NoSite4<NoSite6){
    if(LenRead-1<=min(LenExon[NoSite4],LenExon[NoSite6])){
      if(NoSite5==NoSite4){
        len2<-LenRead-1
      }else{
        len2<-LenRead-LenExon[NoSite5]-1
      }     
    }else{
      len2<-min(LenExon[NoSite4],LenExon[NoSite6])
    }
    
  }
   
  if(NoSite1==NoSite4&&NoSite3==NoSite6){
    len<-max(len1,1)*(max(len2,1)+1)/2/(100*100)
  }else{
    len<-max(len1,1)*max(len2,1)/(100*100)
  }
  
 
  
  }
  
  return(len)
}
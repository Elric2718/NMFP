#####recording votes with consideration about balance of exisiting and non-exisiting 
#####of one exon, the positoin of the exon and short exon. Only consider long isoform
#####with the number of exons no less than 7 and the largest length of bins is 4.

VisualCheck4<-function(ExonPool,LenRead,Type,Signal){
  nExon<-length(ExonPool[1,])
  LenExon<-ExonPool[2,]-ExonPool[1,]+1
  Sindex<-which(LenExon<LenRead)
  Sindex<-c(Sindex,0)#guarantee it contains at least one element
  #weights vector
  VseqP1<-rep(0,nExon)
  VseqN1<-rep(0,nExon)
  
  #non-weights vector
  VseqP2<-rep(0,nExon)
  VseqN2<-rep(0,nExon)
  
  NoSite1<-as.numeric(as.character(strsplit(Type,":")[[1]][1]))
  NoSite2<-as.numeric(as.character(strsplit(Type,":")[[1]][2]))
  
  Ifnorm<-(is.element(NoSite1,Sindex)==FALSE&&is.element(NoSite2,Sindex)==FALSE)|(
    is.element(NoSite1,Sindex)==TRUE&&NoSite1>1&&is.element(NoSite2,Sindex)==FALSE)|(
      is.element(NoSite2,Sindex)==TRUE&&NoSite2<nExon&&is.element(NoSite1,Sindex)==FALSE)|(
      is.element(NoSite1,Sindex)==TRUE&&is.element(NoSite2,Sindex)==TRUE&&NoSite1>1&&NoSite2<nExon)
  
  if(Signal==1){
    
    
  
   # if(Ifnorm==TRUE){
      VseqP1[NoSite1]<-VseqP1[NoSite1]+1
      VseqP1[NoSite2]<-VseqP1[NoSite2]+1

      VseqP2[NoSite1]<-1
      VseqP2[NoSite2]<-1

      if(NoSite1==NoSite2-2){
        VseqN1[NoSite1+1]<-VseqN1[NoSite1+1]-2
        VseqN2[NoSite1+1]<--1
      }
      
      if(NoSite1==NoSite2-3){
        VseqN1[NoSite1+1]<-VseqN1[NoSite1+1]-2
        VseqN2[NoSite1+1]<--1
        VseqN1[NoSite1+2]<-VseqN1[NoSite1+2]-2
        VseqN2[NoSite1+2]<--1
      }
      
      if(NoSite1==NoSite2&&is.element(NoSite1,Sindex)){
        VseqP1[NoSite1-1]<-VseqP1[NoSite1-1]+1
        VseqP1[NoSite1+1]<-VseqP1[NoSite1+1]+1
        
        VseqP2[NoSite1-1]<-1
        VseqP2[NoSite1+1]<-1
      }

    #}
    
    #if((is.element(NoSite1,Sindex)&&NoSite1==1&&NoSite2>1)|(is.element(NoSite2,Sindex)&&NoSite2==nExon&&NoSite1<nExon)){
    #    VseqP1[NoSite1]<-VseqP1[NoSite1]+1
    #    VseqP1[NoSite2]<-VseqP1[NoSite2]+1
        
    #    VseqP2[NoSite1]<-1
    #    VseqP2[NoSite2]<-1
        
    #    if(NoSite1==NoSite2-2){
    #      VseqN1[NoSite1+1]<-VseqN1[NoSite1+1]-2
    #      VseqN2[NoSite1+1]<--1
    #    }
        
    #    if(NoSite1==NoSite2-3){
    #      VseqN1[NoSite1+1]<-VseqN1[NoSite1+1]-2
    #      VseqN2[NoSite1+1]<--1
    #      VseqN1[NoSite1+2]<-VseqN1[NoSite1+2]-2
    #      VseqN2[NoSite1+2]<--1
    #    }
    #}
    
  }else{
    
    if(Ifnorm==TRUE){
      if(NoSite1==NoSite2){
      VseqN1[NoSite1]<-VseqN1[NoSite1]-2.5
      VseqN2[NoSite1]<--1
      }
      
    }else{
      if(is.element(NoSite1,Sindex)==TRUE&&NoSite1==1&&((is.element(NoSite2,Sindex)==FALSE)|
        (is.element(NoSite2,Sindex)==TRUE&&NoSite2>1&&NoSite2<nExon))){
        if(NoSite1+1==NoSite2){
        VseqN1[NoSite1]<-VseqN1[NoSite1]-1
        
        VseqN2[NoSite1]<--1
        }
        
        if(NoSite1+1<NoSite2){
          VseqN1[NoSite1]<-VseqN1[NoSite1]-1/3
          
          VseqN2[NoSite1]<--1
        }
      }
      
      if(is.element(NoSite2,Sindex)==TRUE&&NoSite2==nExon&&((is.element(NoSite1,Sindex)==FALSE)|
         (is.element(NoSite1,Sindex)==TRUE&&NoSite1>1&&NoSite1<nExon))){
        if(NoSite1+1==NoSite2){
          VseqN1[NoSite2]<-VseqN1[NoSite2]-1
          
          VseqN2[NoSite2]<--1
        }
        
        if(NoSite1+1<NoSite2){
          VseqN1[NoSite2]<-VseqN1[NoSite2]-1/3
          
          VseqN2[NoSite2]<--1
        }
      }
      
      if((is.element(NoSite1,Sindex)==TRUE)&&(is.element(NoSite2,Sindex)==TRUE)&&NoSite1==1&&NoSite2==nExon){
        if(NoSite1+1==NoSite2){
          VseqN1[NoSite1]<-VseqN1[NoSite1]-1
          VseqN1[NoSite2]<-VseqN1[NoSite2]-1
          
          VseqN2[NoSite1]<--1
          VseqN2[NoSite2]<--1
        }
        
        if(NoSite1+1<NoSite2){
          VseqN1[NoSite1]<-VseqN1[NoSite1]-1/3
          VseqN1[NoSite2]<-VseqN1[NoSite2]-1/3
          
          VseqN2[NoSite1]<--1
          VseqN2[NoSite2]<--1
        }
      }
      
      
    }
    
  }



Vseq<-list(VseqP1,VseqN1,VseqP2,VseqN2)
names(Vseq)<-c("Positive1","Negative1","Positive2","Negative2")
return(Vseq)
}
AssignMatrix4<-function(chr,strand,GeneStart,GeneEnd,ExonPool,TypesBox,LenRead,input_file){
  nExon<-length(ExonPool[1,])
  nsample<-length(input_file)
  
  ReadMatrix<-matrix(0,length(TypesBox),nsample)
  ###get reads with rbamtools
  for(k in 1:nsample){
    bam<-input_file[k]
    idx<-paste(input_file[k],".bai",sep="")
    reader<-bamReader(bam)
    loadIndex(reader,idx)
    
    ###chr flexibility
    chr<-strsplit(chr,"chr")[[1]]
    chr<-chr[length(chr)]
    chr<-c(chr,paste("chr",chr,sep=""))
    
    ##modify information of chr and strand
    if(length(which(getRefData(reader)$SN==chr[1]))>0){
      chr_index<-which(getRefData(reader)$SN==chr[1])
      chr<-chr[1]
    }else{
      chr_index<-which(getRefData(reader)$SN==chr[2])
      chr<-chr[2]
    }
    chr<-c(chr,getRefData(reader)$ID[chr_index])
    
    ##extract data and process(such process can only be applied to simulation)
    coords<-c(as.numeric(chr[2]),GeneStart,GeneEnd)
    range<-bamRange(reader,coords)
    rdf<-as.data.frame(range)
  
    if(nrow(rdf)==0){next}
    ##extract useful information
    #flag information
    FLAG<-t(matrix(intToBits(rdf$flag),nrow=32)[1:11,])
    #exclude exons which do not coordiate whith strand
    if(strand[1]=="-"){
      SignIndex<-union(
        intersect(which(FLAG[,7]==1),which(as.logical(rdf[,9])==TRUE)),
        intersect(which(FLAG[,8]==1),which(as.logical(rdf[,9])==FALSE)))
    }else{SignIndex<-union(
      intersect(which(FLAG[,7]==1),which(as.logical(rdf[,9])==FALSE)),
      intersect(which(FLAG[,8]==1),which(as.logical(rdf[,9])==TRUE)))
    }
    
    SIGN_CHR<-intersect(SignIndex,which(rdf[,1]==as.numeric(chr[2])))
    if(length(SIGN_CHR)==0){next}
    ReadData<-rdf[SIGN_CHR,c(2,3,4,5,8)]
    FLAG<-FLAG[SIGN_CHR,]
    
    ##Cigar information
    CIGAR_LEN<-apply(t(seq(nrow(ReadData))),2,function(i){
      if(ReadData[i,2]>1){
        CIGAR<-unlist(strsplit(ReadData[i,3],""))
        CIGAR_CLASS<-which(is.element(CIGAR,LETTERS))
        CIGAR_CLASS<-rbind(rbind(CIGAR_CLASS,c(1,CIGAR_CLASS[1:(length(CIGAR_CLASS)-1)]+1)),CIGAR_CLASS-1)
        
        
        CIGAR<-rbind(CIGAR[CIGAR_CLASS[1,]],
                     apply(t(seq(ncol(CIGAR_CLASS))),2,function(x){
                       return(as.numeric(paste(CIGAR[CIGAR_CLASS[2,x]:CIGAR_CLASS[3,x]],collapse="")))
                     }))
        LEN_SEQ<-sum(apply(t(seq(ncol(CIGAR))),2,function(x){
          if(is.element(CIGAR[1,x],c("M","N","D"))){return(as.numeric(CIGAR[2,x]))}else{return(0)}
        }))}else{
          LEN_SEQ<-as.numeric(unlist(strsplit(ReadData[i,3],"M")))
        }
      return(LEN_SEQ)})
    
    ReadData<-cbind(ReadData[,1],ReadData[,1]+CIGAR_LEN-1)

    
    
    ##assign reads
    
    #Type like 1:1:1:1:1:1
    #ReadMatrix[,k]<-ReadCount2(ReadData,ExonPool+GeneStart-1,TypesBox,LenRead)
    #Type like 1:1
    ReadMatrix[,k]<-ReadCount3(ReadData,ExonPool+GeneStart-1,TypesBox,LenRead)
    
    
    ##debug the ReadMatrix
    #IndexReadCount<-NULL
    #for(i in 1:(39*(nExon-4)+41){
    #  
    #    IndexReadCount<-union(IndexReadCount,which(apply(t(seq(nPair)),2,function(x){return(ReadCount(ReadData[x,],ExonPool,TypesBox[i]))})==T))
    
    #IndexReadCount<-IndexReadCount[order(IndexReadCount)]
    #print(paste0("reading sample ",k,"."))
  }
  
  return(ReadMatrix)
}


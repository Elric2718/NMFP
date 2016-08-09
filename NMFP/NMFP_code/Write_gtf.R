Write_gtf<-function(Candidates,ExonPool,chr,strand,Gene_Name,GeneStart){

  
Result_to_gtf<-NULL
for(i in 1:length(Candidates)){
  Boundary<-as.matrix(ExonPool[,which(unlist(strsplit(Candidates[i],split=""))==1)],nrow=2)
  if(ncol(Boundary)>1){
    j<-1
    j2<-1
    while(j<=(ncol(Boundary)-1)){
      if(Boundary[2,j]>=Boundary[1,j+1]-1){
        j2<-j+1
        while(j2<=(ncol(Boundary)-1)&&Boundary[2,j2]>=Boundary[1,j2+1]-1){
          j2<-j2+1
        }
        for(t in j:j2){
          Boundary[2,t]<-Boundary[2,j2]
          Boundary[1,t]<-Boundary[1,j]
        }
      }
      j<-j2+1
      j2<-j
    }
    Boundary<-as.matrix(Boundary,nrow=2)
    Boundary<-Boundary[,which(duplicated(t(Boundary))==FALSE)]
    Boundary<-as.matrix(Boundary,nrow=2)
  }
  Boundary<-Boundary+GeneStart-1
  
  Tx_gtf<-data.frame(
    chrom=chr,
    Source="NMFP",
    feature=c("transcript",rep("exon",ncol(Boundary))),
    start=c(min(Boundary),Boundary[1,]),
    end=c(max(Boundary),Boundary[2,]),
    score=".",
    strand=strand[1],
    frame=".",
    attributes=paste('gene_id \"',Gene_Name,'\"; ','transcript_id \"',Gene_Name,"_",i,'\";',sep="")
  )
  
  Result_to_gtf<-rbind(Result_to_gtf,Tx_gtf)
  
}

return(Result_to_gtf)
}



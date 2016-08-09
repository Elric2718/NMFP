### need complete annotation like that from Ensemble
library(optparse)
makegenelist<-function(file,CHR,mode=1){
  # A function that returns a matrix, each row as a transcript of a gene
  # Args: 
  #  file: the annotation file
  #  CHR: the chromosome on which the genelist is made. if CHR =  NULL, then make the 
  #       the genelist on the whole annotationl file
  #  mode: the mode is used to distinguish different kinds of annotation files.
  anno<-read.table(file,fill=TRUE)
  if(is.null(CHR)==FALSE){
    anno<-anno[which(as.character(anno[,1])==CHR),]
  }
  Gene<-unique(as.character(anno[,10]))
  if(mode==1||mode==3){
    Tx<-unique(as.character(anno[,13]))
    GENE<-t(apply(t(t(seq(length(Tx)))),1,function(i){
      chr<-as.character(anno[which(as.character(anno[,13])==Tx[i])[1],1])
      start_site<-min(as.numeric(as.character(anno[which(as.character(anno[,13])==Tx[i]),4])))
      end_site<-max(as.numeric(as.character(anno[which(as.character(anno[,13])==Tx[i]),5])))
      strand<-as.character(anno[which(as.character(anno[,13])==Tx[i])[1],7])
      gene_id<-as.character(anno[which(as.character(anno[,13])==Tx[i])[1],10])
      print(paste(i,length(Tx),sep="/"))
      return(c(chr,start_site,end_site,strand,Tx[i],gene_id))
    }))
  }else{
    GENE<-anno[which(as.character(anno[,3])=="transcript"),]
    GENE<-GENE[,c(1,4,5,7,16,10)]
  }
 
## Add indexes 
  if(mode<=2){
    GENE<-cbind(apply(t(seq(nrow(GENE))),2,
                      FUN=function(i){return(paste("G",which(Gene==as.character(GENE[i,6])),sep=""))}),GENE)
  }else{
    gene_index1<-1
    gene_index2<-1
    #gene_index_temp<-1
    count<-1
    Gene_Index<-NULL
    while(gene_index1<=nrow(GENE)){
      gene_start<-as.numeric(as.character(GENE[gene_index1,2]))
      gene_end<-as.numeric(as.character(GENE[gene_index2,3]))
      gene_strand<-as.character(GENE[gene_index1,4])
      gene_index2<-intersect(intersect(which(as.numeric(as.character(GENE[,2]))>=gene_start),
                                       which(as.numeric(as.character(GENE[,2]))<=gene_end)),
                             which(as.character(GENE[,4])==gene_strand))
      gene_end<-max(as.numeric(as.character(GENE[gene_index2,3])))
      gene_index2<-gene_index2[length(gene_index2)]
      gene_index_temp<-gene_index2
      ifRepeat<-TRUE
      while(ifRepeat){
        gene_index2<-intersect(intersect(which(as.numeric(as.character(GENE[,2]))>=gene_start),
                                         which(as.numeric(as.character(GENE[,2]))<=gene_end)),
                               which(as.character(GENE[,4])==gene_strand))
        gene_end<-max(as.numeric(as.character(GENE[gene_index2,3])))
        gene_index2<-gene_index2[length(gene_index2)]
        if(gene_index_temp==gene_index2){
          ifRepeat<-FALSE
        }else{
          gene_index_temp<-gene_index2
        }
      }
      Gene_Index<-c(Gene_Index,paste("G",rep(count,(gene_index2-gene_index1+1)),sep=""))
      count<-count+1
      gene_index1<-gene_index2+1
      gene_index2<-gene_index1
      print(gene_index1)
    }
    GENE<-cbind(Gene_Index,GENE)
  }
  
  colnames(GENE)<-c("gene_index","chr","start","end","strand","transcript_id","gene_id")
  return(GENE)
}

## arguments
option_list <- list(
  make_option(c("--anno_path"), type="character", default=NULL,
              help="Set the path to the annotation file.",
              metavar="anno_path"),
  make_option(c("--chr"), type="character", default=NULL,
              help="The chromosome to make gene list on.",
              metavar="chr"),
  make_option(c("--output_path"), type="character", default=NULL,
              help="Set the path to the output of NMFP.",
              metavar="output_path")
)
opt <- parse_args(OptionParser(option_list=option_list))


anno_path <- opt$anno_path
chr <- opt$chr
output_path <- opt$output_path



if(!file.exists(output_path)){
GENE <- makegenelist(anno_path, CHR = chr, mode = 1)
write.table(GENE, file = output_path, quote=F,row.names=F,col.names=F,sep="\t")
}
# load packages and .R files
library(Biobase)
library(parallel)
library(rbamtools)
library(compiler)
library(optparse)

## arguments
option_list <- list(
  make_option(c("--code_path"), type="character", default=NULL,
              help="Set the path to the codes.",
              metavar="code_path"),
  make_option(c("--data_path"), type="character", default=NULL,
              help="Set the path to the data of .bam.",
              metavar="data_path"),
  make_option(c("--anno_path"), type="character", default=NULL,
              help="Set the path to the annotation file.",
              metavar="anno_path"),
  make_option(c("--genelist_path"), type="character", default=NULL,
              help="Set the path to the genelist.",
              metavar="genelist_path"),
  make_option(c("--if_rank_anno"), type="logical", default=TRUE,
              help="If use annotation to estimate the rank of NMFP.",
              metavar="if_rank_anno"),
  make_option(c("--nrun"), type="integer", default=100,
              help="Set the number of repitions for NMFP",
              metavar="nrun"),
  make_option(c("--ncores"), type="integer", default=1,
              help="Number of cores to use [default %default]",
              metavar="ncores"),
  make_option(c("--LenRead"), type="integer", default=76,
              help="Length of the reads.",
              metavar="LenRead"),
  make_option(c("--output_path"), type="character", default=NULL,
              help="Set the path to the output of NMFP.",
              metavar="output_path")
)
opt <- parse_args(OptionParser(option_list=option_list))


code_path <- opt$code_path
data_path <- opt$data_path
anno_path <- opt$anno_path
genelist_path <- opt$genelist_path
if_rank_anno <- opt$if_rank_anno
nrun <- opt$nrun
ncores <- opt$ncores
LenRead <- opt$LenRead
output_path <- opt$output_path


source(paste(code_path,"ExonBoundary.R",sep=""))
source(paste(code_path,"FindTypes2.R",sep=""))
source(paste(code_path,"PreSelection.R",sep=""))
source(paste(code_path,"VoteforIsoform.R",sep=""))
source(paste(code_path,"ncNMF2.R",sep=""))
source(paste(code_path,"VisualizeG5.R",sep=""))
source(paste(code_path,"VisualCheck.R",sep=""))
source(paste(code_path,"VisualCheck4.R",sep=""))
source(paste(code_path,"RankDetermine2.R",sep=""))
source(paste(code_path,"Write_gtf.R",sep=""))
source(paste(code_path,"ContradictTable.R",sep=""))
source(paste(code_path,"Normalization.R",sep=""))
source(paste(code_path,"PairReadLength.R",sep=""))
source(paste(code_path,"AssignMatrix4.R",sep=""))
source(paste(code_path,"ReadCount3.R",sep=""))

cmp_ncNMF2<-cmpfun(ncNMF2)

#setwd("~/Desktop/NMFP/Example2/")



### list files
file_name <- list.files(path=data_path)

### grep .sorted.bam
input_file <- file_name[grep(pattern= "sorted.bam$", file_name)]
input_file <- paste0(data_path, input_file)

###prepare files
anno<-read.table(anno_path,fill=TRUE)
genelist<-read.table(file=genelist_path)

## sets for gene indexes
geneindex <- as.character(unique(genelist[,1]))



NMFResult<-mclapply(X= length(geneindex),mc.cores=ncores,FUN=function(k){
  
  ###NMFP
  Gene_Index <- geneindex[k]
  Gene_Name <- as.character(genelist[which(as.character(genelist[,1])==Gene_Index)[1],7])
  
  ExonList<-ExonBoundary(genelist=genelist, Anno=anno,Gene=Gene_Index,LenRead=LenRead,ifvis=FALSE,Emode=2,Anno_mode=2)
  chr<-ExonList$chr
  strand<-ExonList$strand
  strand<-c(strand,switch(strand,"-"=TRUE,"+"=FALSE))
  GeneStart<-ExonList$GeneStart
  GeneEnd<-ExonList$GeneEnd
  ExonPool<-ExonList$ExonPool
  nExon<-length(ExonPool[1,])
  TypesBox_two<-FindTypes2(nExon)
  LenExon<-ExonPool[2,]-ExonPool[1,]+1
  AnnoTrans<-ExonList$annoTrans
  #AnnoTrans
  
  ReadBin_two<-AssignMatrix4(chr=chr,strand=strand,GeneStart=GeneStart,GeneEnd=GeneEnd,ExonPool=ExonPool,TypesBox=TypesBox_two,LenRead=LenRead,input_file=input_file)
  rownames(ReadBin_two)<-TypesBox_two
  
  ReadBin_two<-ReadBin_two[,which(apply(ReadBin_two,2,sum)>0)]
  
  ###normalization
  NormMatrix<-Normalization(ReadBin_two,ExonPool+GeneStart-1,76,TypesBox_two)
  
  rank0<-length(ExonList$AnnoTrans)
  if(if_rank_anno){
    rank <- rank0
  }else{
    rank <- 0
  }
  
  BIsoform<-VoteforIsoform(rank=rank,nrun=nrun,alpha0=1,
                           NormMatrix=NormMatrix,TypesBox=TypesBox_two,ExonPool=ExonPool,
                           LenRead=LenRead,Gcutoff=0.4,mode=1,
                           input_file=code_path)
  IsoformSet<-BIsoform$IsoformSet
  rank<-BIsoform$rank
  
  
  ### post selection
  if(length(names(IsoformSet))>10){
    PreSelection1_1<-PreSelection(Candidates=names(IsoformSet),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,mode=1)
    PreSelection1_2<-PreSelection(Candidates=names(IsoformSet),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,mode=2)
    PreSelection1_3<-PreSelection(Candidates=names(IsoformSet),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,mode=3)
    PreSelection1_4<-PreSelection(Candidates=names(IsoformSet),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,mode=4)
  }else{
    PreSelection1_1<-PreSelection(Candidates=names(IsoformSet),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,ifPart=TRUE)
    PreSelection1_2<-PreSelection(Candidates=names(IsoformSet),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,ifPart=TRUE)
    PreSelection1_3<-PreSelection(Candidates=names(IsoformSet),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,ifPart=TRUE)
    PreSelection1_4<-PreSelection(Candidates=names(IsoformSet),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,ifPart=TRUE)
  }
  if(length(which(IsoformSet>=2))>10){
    PreSelection2_1<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=2)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,mode=1)
    PreSelection2_2<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=2)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,mode=2)
    PreSelection2_3<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=2)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,mode=3)
    PreSelection2_4<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=2)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,mode=4)
  }else{
    PreSelection2_1<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=2)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,ifPart=TRUE)
    PreSelection2_2<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=2)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,ifPart=TRUE)
    PreSelection2_3<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=2)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,ifPart=TRUE)
    PreSelection2_4<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=2)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,ifPart=TRUE)
  }
  if(length(which(IsoformSet>=5))>10){
    PreSelection3_1<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=5)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,mode=1)
    PreSelection3_2<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=5)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,mode=2)
    PreSelection3_3<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=5)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,mode=3)
    PreSelection3_4<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=5)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,mode=4)
  }else{
    PreSelection3_1<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=5)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,ifPart=TRUE)
    PreSelection3_2<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=5)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,ifPart=TRUE)
    PreSelection3_3<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=5)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,ifPart=TRUE)
    PreSelection3_4<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=5)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,ifPart=TRUE)
  }
  if(length(which(IsoformSet>=10))>10){
    PreSelection4_1<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=10)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,mode=1)
    PreSelection4_2<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=10)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,mode=2)
    PreSelection4_3<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=10)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,mode=3)
    PreSelection4_4<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=10)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,mode=4)
  }else{
    PreSelection4_1<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=10)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,ifPart=TRUE)
    PreSelection4_2<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=10)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,ifPart=TRUE)
    PreSelection4_3<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=10)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,ifPart=TRUE)
    PreSelection4_4<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=10)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,ifPart=TRUE)
  }
  if(length(which(IsoformSet>=25))>10){
    PreSelection5_1<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=25)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,mode=1)
    PreSelection5_2<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=25)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,mode=2)
    PreSelection5_3<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=25)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,mode=3)
    PreSelection5_4<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=25)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,mode=4)
  }else{
    PreSelection5_1<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=25)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,ifPart=TRUE)
    PreSelection5_2<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=25)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,ifPart=TRUE)
    PreSelection5_3<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=25)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,ifPart=TRUE)
    PreSelection5_4<-PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=25)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead,ifPart=TRUE)
  }
  
  
  #length(PreSelection1)
  preslct1_1<-NULL
  preslct1_2<-NULL
  preslct1_3<-NULL
  preslct1_4<-NULL
  for(i in 1:length(AnnoTrans)){
    preslct1_1<-c(preslct1_1,PreSelection1_1[which(is.element(PreSelection1_1,AnnoTrans[i]))])
    preslct1_2<-c(preslct1_2,PreSelection1_2[which(is.element(PreSelection1_2,AnnoTrans[i]))])
    preslct1_3<-c(preslct1_3,PreSelection1_3[which(is.element(PreSelection1_3,AnnoTrans[i]))])
    preslct1_4<-c(preslct1_4,PreSelection1_4[which(is.element(PreSelection1_4,AnnoTrans[i]))])
  }
  #print(preslct1)
  
  
  
  #length(PreSelection2)
  preslct2_1<-NULL
  preslct2_2<-NULL
  preslct2_3<-NULL
  preslct2_4<-NULL
  for(i in 1:length(AnnoTrans)){
    preslct2_1<-c(preslct2_1,PreSelection2_1[which(is.element(PreSelection2_1,AnnoTrans[i]))])
    preslct2_2<-c(preslct2_2,PreSelection2_2[which(is.element(PreSelection2_2,AnnoTrans[i]))])
    preslct2_3<-c(preslct2_3,PreSelection2_3[which(is.element(PreSelection2_3,AnnoTrans[i]))])
    preslct2_4<-c(preslct2_4,PreSelection2_4[which(is.element(PreSelection2_4,AnnoTrans[i]))])
  }
  #print(preslct2)
  
  
  
  #length(PreSelection3)
  preslct3_1<-NULL
  preslct3_2<-NULL
  preslct3_3<-NULL
  preslct3_4<-NULL
  for(i in 1:length(AnnoTrans)){
    preslct3_1<-c(preslct3_1,PreSelection3_1[which(is.element(PreSelection3_1,AnnoTrans[i]))])
    preslct3_2<-c(preslct3_2,PreSelection3_2[which(is.element(PreSelection3_2,AnnoTrans[i]))])
    preslct3_3<-c(preslct3_3,PreSelection3_3[which(is.element(PreSelection3_3,AnnoTrans[i]))])
    preslct3_4<-c(preslct3_4,PreSelection3_4[which(is.element(PreSelection3_4,AnnoTrans[i]))])
  }
  #print(preslct3)
  
  #length(PreSelection4)
  preslct4_1<-NULL
  preslct4_2<-NULL
  preslct4_3<-NULL
  preslct4_4<-NULL
  for(i in 1:length(AnnoTrans)){
    preslct4_1<-c(preslct4_1,PreSelection4_1[which(is.element(PreSelection4_1,AnnoTrans[i]))])
    preslct4_2<-c(preslct4_2,PreSelection4_2[which(is.element(PreSelection4_2,AnnoTrans[i]))])
    preslct4_3<-c(preslct4_3,PreSelection4_3[which(is.element(PreSelection4_3,AnnoTrans[i]))])
    preslct4_4<-c(preslct4_4,PreSelection4_4[which(is.element(PreSelection4_4,AnnoTrans[i]))])
  }
  #print(preslct4)
  
  preslct5_1<-NULL
  preslct5_2<-NULL
  preslct5_3<-NULL
  preslct5_4<-NULL
  for(i in 1:length(AnnoTrans)){
    preslct5_1<-c(preslct5_1,PreSelection5_1[which(is.element(PreSelection5_1,AnnoTrans[i]))])
    preslct5_2<-c(preslct5_2,PreSelection5_2[which(is.element(PreSelection5_2,AnnoTrans[i]))])
    preslct5_3<-c(preslct5_3,PreSelection5_3[which(is.element(PreSelection5_3,AnnoTrans[i]))])
    preslct5_4<-c(preslct5_4,PreSelection5_4[which(is.element(PreSelection5_4,AnnoTrans[i]))])
  }
  
  
  ##postselection for nExon<=10
  PreSelection_all<-list(PreSelection1_1,PreSelection1_2,PreSelection1_3,PreSelection1_4,
                         PreSelection2_1,PreSelection2_2,PreSelection2_3,PreSelection2_4,
                         PreSelection3_1,PreSelection3_2,PreSelection3_3,PreSelection3_4,
                         PreSelection4_1,PreSelection4_2,PreSelection4_3,PreSelection4_4,
                         PreSelection5_1,PreSelection5_2,PreSelection5_3,PreSelection5_4)
  preslct_all<-list(preslct1_1,preslct1_2,preslct1_3,preslct1_4,
                    preslct2_1,preslct2_2,preslct2_3,preslct2_4,
                    preslct3_1,preslct3_2,preslct3_3,preslct3_4,
                    preslct4_1,preslct4_2,preslct4_3,preslct4_4,
                    preslct5_1,preslct5_2,preslct5_3,preslct5_4)
  
  PreLen<-listLen(PreSelection_all)
  PreOrder<-order(PreLen)
  
  
  
  IndexCand<-which(PreLen[PreOrder]<=length(AnnoTrans)*10)
  if(length(IndexCand)==0){IndexCand<-1}else{IndexCand<-IndexCand[length(IndexCand)]}
  Candidates<-PreSelection_all[[PreOrder[IndexCand]]]
  candidates<-preslct_all[[PreOrder[IndexCand]]]
  
  
  
  
  print(paste("Process: Gene ",k,"/",length(geneindex),
              ". The number of subexons is ",nExon,". The number of annotated isoforms is ",length(AnnoTrans),". ",
              "The number of candidates is ",length(Candidates),".",sep=""))
  
  
  ## NMFP.gtf
  Result_to_gtf<-Write_gtf(Candidates=Candidates,ExonPool=ExonPool,chr=chr,strand=strand,Gene_Name=Gene_Name,GeneStart=GeneStart)
  write.table(Result_to_gtf,file=output_path,quote=F,row.names=F,col.names=F,sep="\t",append=TRUE)
  
})


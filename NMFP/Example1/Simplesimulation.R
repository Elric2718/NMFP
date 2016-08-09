###############################################################################
######################### NMFP for lowly expressed gene #######################
###############################################################################
###### load package
library(dplyr)
library(compiler)
library(data.table)
library(reshape2)
library(ggplot2)
library(grid)


###load codes
########### please first set the path to the 'NMFP_code/NMFP' dir ############
code_file <- "~/Desktop/NMFP_code/NMFP"
source(paste(code_file,"/ExonBoundary.R",sep=""))
source(paste(code_file,"/FindTypes2.R",sep=""))
source(paste(code_file,"/PreSelection.R",sep=""))
source(paste(code_file,"/VoteforIsoform.R",sep=""))
source(paste(code_file,"/ncNMF2.R",sep=""))
source(paste(code_file,"/VisualizeG5.R",sep=""))
source(paste(code_file,"/VisualCheck.R",sep=""))
source(paste(code_file,"/VisualCheck4.R",sep=""))
source(paste(code_file,"/RankDetermine2.R",sep=""))
source(paste(code_file,"/Write_gtf.R",sep=""))
source(paste(code_file,"/ContradictTable.R",sep=""))
source(paste(code_file,"/Normalization.R",sep=""))
source(paste(code_file,"/PairReadLength.R",sep=""))
source(paste(code_file,"/AssignMatrix4.R",sep=""))
source(paste(code_file,"/ReadCount3.R",sep=""))

cmp_ncNMF2<-cmpfun(ncNMF2)

################## construct the gene ##################
nExon <- 8
nIsof <- 3
LenRead <- 76

exon_length <- seq(100,500, by = 50)
GenerateExon <- function(nExon, exon_length){
##### Genereate ExonPool (generating exons)
  nchoice <- length(exon_length)
  ExonPool <- rbind(rep(1, nExon), sample(exon_length, nExon, replace = TRUE))
  for(i in 1 : (nExon-1)){
    ExonPool[,(i+1)] <- ExonPool[,(i+1)] + ExonPool[2,i]
  }
  return(ExonPool)
}

ExonPool <- GenerateExon(nExon, exon_length)
ELen <- ExonPool[2,] - ExonPool[1,] + 1

## types of bins, like (1,1), (1,2)
TypesBox <- FindTypes2(nExon)

## get the effective length of each bin
effective_len <- sapply(TypesBox, function(type) 
                                PairReadLength(ExonPool, LenRead, Type = type))

## generating the gene with three isoforms. The length of gene is of 8 exons.
ISOF_exon <- matrix(0, nrow = nExon , ncol = nIsof)
ISOF_exon[, 1] <- c(1, 1, 0, rep(1, 5))
ISOF_exon[, 2] <- c(0, 1, 1, rep(1, 5))
ISOF_exon[, 3] <- c(1, 0, 0, rep(1, 5))

## The profile for true isoforms
AnnoTrans <- apply(ISOF_exon, 2, function(isof){
  paste(isof, collapse="")
})

## adding junction bins
ISOF <- apply(ISOF_exon, 2, function(isof){
  nExon <- length(isof)
  bin_level1 <- isof
  bin_level2 <- isof[-nExon] * isof[-1]
  bin_level3 <- isof[-c(nExon-1, nExon)] * isof[-c(1, 2)]
  bin_level4 <- isof[-seq(nExon-2, nExon)] * isof[-seq(1, 3)]
  return(c(bin_level1, bin_level2, bin_level3, bin_level4))
}) 

rownames(ISOF) <- TypesBox

## The probablity of reads falling into one bin given the isoform
length_ratio <- apply(ISOF, 2, function(isof){
  isof * effective_len/sum(isof * effective_len)
})

GenerateReps <- function(ISOF, nrep = 10, length_ratio, depth){
  ##### Generate nrep replicates of genes given ISOF
  stopifnot(nrep == length(depth))
  nbin <- nrow(ISOF)
  nIsof <- ncol(ISOF)
  
  sapply(seq(nrep), function(i){
    ### simulatre gene expression
    isof_ratio <- runif(nIsof)
    isof_ratio <- isof_ratio/sum(isof_ratio)
    express_level <- ISOF %*% diag(isof_ratio) * depth[i]

    ### add noise
    bin_noise <- matrix(rnorm(n = nbin * nIsof , mean = 0, sd = 1), 
                        nrow = nbin) %*% 
                 diag(isof_ratio)
    noise_express <- express_level + bin_noise
    noise_express[noise_express < 0] <- 0 
    ### simulate sequencing
    reads_level <- rowSums(noise_express * length_ratio)
    return(reads_level)
  })
}



## generate the expression level
DepthSelect <- function(ntotal = 10, high_level = 5){
  depth <- c(3, rep(high_level, (ntotal-1)))
  return(depth)
}


### simulation part
RECALL <- NULL
nCand <- NULL
for(k in c(0.5,seq(10))){
recall <- NULL
nCandidates<- NULL

for(j in 1:50){
  
  depth <- DepthSelect(ntotal = 10, high_level = k)

  print(depth)
  REPs <- GenerateReps(ISOF, 10, length_ratio, depth)
  
  NormMatrix<-Normalization(REPs, ExonPool, LenRead, TypesBox)
  
  
  BIsoform<-try(VoteforIsoform(rank=3,nrun=100,alpha0=0.1,
                               NormMatrix=NormMatrix,TypesBox=TypesBox,ExonPool=ExonPool,
                               LenRead=LenRead,Gcutoff=0.4,mode=1,
                               input_file=code_file))
  if(inherits(BIsoform, "try-error"))
  {
    #error handling code, maybe just skip this iteration using
    next
  }
  IsoformSet<-BIsoform$IsoformSet
  #rank<-BIsoform$rank
  
  Candidates <- names(IsoformSet[1:30])
  #Candidates <- PreSelection(Candidates=names(IsoformSet[which(IsoformSet>=70)]),SumRead=apply(NormMatrix,1,mean),ExonPool=ExonPool,LenRead=LenRead)
  candidates <- NULL
  for(i in 1:length(AnnoTrans)){
    candidates<-c(candidates,Candidates[which(is.element(Candidates,AnnoTrans[i]))])
  }
  
  nCandidates <- c(nCandidates, max(which(is.element(names(IsoformSet),AnnoTrans))))
  recall <- c(recall,length(candidates)/nIsof)
  print(j)
}
nCand[[as.character(k)]] <- nCandidates
RECALL[[as.character(k)]] <- recall

}



candmat <-
apply(t(c(0.5,seq(10))),2,function(i){
  vec <- rep(NA,50)
  vec[1:length(nCand[[as.character(i)]])] <- nCand[[as.character(i)]]
  return(vec)
  
}) %>%
  data.table() %>%
  setnames(paste(c(0.5,seq(10)))) %>%
  melt() %>%
  setnames(c("exp_level", "n_candidates"))


ggplot(candmat, aes(x= exp_level, y = n_candidates)) +
  geom_boxplot() + 
  xlab("Expression Level") +
  ylab("Number of Isoform Candidates") + 
  theme_bw() +
  theme(
    axis.text = element_text(size = rel(0.9), colour = "grey50"),
    axis.title = element_text(size = rel(1.2), colour = "black"),
    axis.ticks = element_line(size = rel(1), colour = "black"),
    axis.ticks.length = unit(0.15, "cm"), 
    axis.ticks.margin = unit(0.1, "cm"),
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text = element_text(colour="black", size=12, face="bold"),  
    plot.title = element_text(size = rel(1.2), colour = "grey50"))



recallmat <-
  apply(t(c(0.5,seq(10))),2,function(i){
    vec <- rep(NA,50)
    vec[1:length(RECALL[[as.character(i)]])] <- RECALL[[as.character(i)]]
    return(vec)
    
  }) %>%
  data.table() %>%
  setnames(paste(c(0.5,seq(9)))) %>%
  melt() %>%
  setnames(c("exp_level", "recall"))


ggplot(recallmat, aes(x= exp_level, y = recall)) +
  geom_boxplot() + 
  xlab("Expression Level") +
  ylab("Recall") + 
  theme_bw() +
  theme(
    axis.text = element_text(size = rel(0.9), colour = "grey50"),
    axis.title = element_text(size = rel(1.5), colour = "black"),
    axis.ticks = element_line(size = rel(1), colour = "black"),
    axis.ticks.length = unit(0.15, "cm"), 
    axis.ticks.margin = unit(0.1, "cm"),
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text = element_text(colour="black", size=12, face="bold"),  
    plot.title = element_text(size = rel(1.5), colour = "grey50"))


##########################################################################################
##################################  Annotation Processing ################################
##########################################################################################

# ------------------ need complete annotation like that from Ensemble ------------------ #
library(optparse)
library(dplyr)

# ----------------------------------- useful functions --------------------------------- #
AttrExtract <- function(attr, id){
  # Extract a vector given the id from the attribute field
  #  Args:
  #   attr: attribute field in annotation file
  #   id: the id in the attribute field
  extract_term <- sapply(attr, function(x){
    y = grep(id, unlist(strsplit(x, ";")), value = TRUE) %>%
      unlist %>%
      strsplit(split = paste0(" ", id, " ")) %>%
      unlist 
    y[2]
  })
  names(extract_term) <- NULL
  extract_term
}


TrimGTF <- function(file){
  # Returns a trimed.gtf which only contains the information needed for NMFP
  # Args:
  #  file: annotation file
  anno <- read.table(file,fill=TRUE, sep = "\t")
  
  attribute <- as.character(anno[, 9])
  
  data.frame(seqname = as.character(anno[, 1]), 
             feature = as.character(anno[, 3]),
             start = as.numeric(as.character(anno[, 4])),
             end = as.numeric(as.character(anno[, 5])),
             strand = as.character(anno[, 7]),
             transcript_id = AttrExtract(attribute, "transcript_id"),
             gene_id = AttrExtract(attribute, "gene_id"))
  
}


MakeGeneList <- function(file){
  # Returns a gene list, each row is a transcript
  # Args: 
  #  file: the trimmed annotation file from TrimGTF
  
  anno <- read.table(file, fill=TRUE, sep = "\t")
  
  gene <- unique(as.character(anno[, 7]))
  
  if(!"transcript" %in% unique(as.character(anno[, 2]))){
    tx_id <- unique(as.character(anno[, 6])) 
    
    genelist <- sapply(seq(length(tx_id)), function(i){
      tx_index <- which(as.character(anno[ ,6]) == tx_id[i])
      chr <- as.character(anno[tx_index[1], 1])
      start <- min(as.numeric(as.character(anno[tx_index, 3])))
      end <- max(as.numeric(as.character(anno[tx_index, 4])))
      strand <- as.character(anno[tx_index[1], 5])
      gene_id <- as.character(anno[tx_index[1], 7])
      
      print(paste0("Making gene list: ",i,"/",length(tx_id),"."))
      return(c(chr, start, end, strand, tx_id[i], gene_id))
    }) %>% t()
  }else{
    genelist <- anno[as.character(anno[, 2]) == "transcript", -2]
  }
  
  ## Add indexes 
  ## Note: For some bad annotations which give the same gene_id as the transcript_id,
  ##       like USCS annotation, its gene_id is of no use. Thus it's needed to first 
  ##       figure out which transcripts are from the same genes. 
  genelist <- cbind(sapply(seq(nrow(genelist)), FUN = function(i){ paste0("G", which(gene == as.character(genelist[i, 6])))}),
                    genelist)
  
  colnames(genelist)<-c("gene_index", "chr", "start", "end", "strand", "transcript_id", "gene_id")
  return(genelist)
}

#------------------------------------- Main --------------------------------------------#
## arguments
option_list <- list(
  make_option(c("--anno_path"), type="character", default=NULL,
              help="Set the path to the annotation file.",
              metavar="anno_path")
)
opt <- parse_args(OptionParser(option_list=option_list))


anno_path <- opt$anno_path
trim_path <- paste0(unlist(strsplit(anno_path, ".gtf")), "_trimmed.txt")
genelist_path <- paste0(unlist(strsplit(anno_path, ".gtf")), "_genelist.txt")

## make the trimmed annoation
if(!file.exists(trim_path)){
  trimmed_anno <- TrimGTF(anno_path)
  write.table(trimmed_anno, file = trim_path, quote=F, row.names=F, col.names=F, sep="\t")
}

## make the gene list
if(!file.exists(genelist_path)){
  genelist <- MakeGeneList(trim_path)
  write.table(genelist, file = genelist_path, quote=F,row.names=F,col.names=F,sep="\t")
}
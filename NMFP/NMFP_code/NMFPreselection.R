# load packages and .R files
library(parallel)
library(rbamtools)
library(compiler)
library(optparse)
library(dplyr)

## arguments
option_list <- list(
  make_option(c("--code_path"), type="character", default = NULL,
              help = "Set the path to the codes.",
              metavar = "code_path"),
  make_option(c("--data_path"), type = "character", default = NULL,
              help = "Set the path to the data of .bam.",
              metavar = "data_path"),
  make_option(c("--anno_path"), type = "character", default = NULL,
              help="Set the path to the annotation file.",
              metavar = "anno_path"),
  make_option(c("--genelist_path"), type = "character", default = NULL,
              help = "Set the path to the genelist.",
              metavar = "genelist_path"),
  make_option(c("--if_rank_anno"), type = "logical", default = TRUE,
              help = "If use annotation to estimate the rank of NMFP.",
              metavar = "if_rank_anno"),
  make_option(c("--nrun"), type="integer", default = 100,
              help = "Set the number of repitions for NMFP",
              metavar = "nrun"),
  make_option(c("--ncores"), type = "integer", default = 1,
              help = "Number of cores to use [default %default]",
              metavar = "ncores"),
  make_option(c("--len_read"), type = "integer", default = 76,
              help = "Length of the reads.",
              metavar = "len_read"),
  make_option(c("--output_path"), type = "character", default = NULL,
              help = "Set the path to the output of NMFP.",
              metavar = "output_path")
)
opt <- parse_args(OptionParser(option_list = option_list))


code_path <- opt$code_path
data_path <- opt$data_path
anno_path <- opt$anno_path
genelist_path <- opt$genelist_path
if_rank_anno <- opt$if_rank_anno
nrun <- opt$nrun
ncores <- opt$ncores
len_read <- opt$len_read
output_path <- opt$output_path

### load the source codes
source(paste(code_path, "helper.R", sep = ""))
cmp_ncNMF2 <- cmpfun(ncNMF2)

### list files
file_name <- list.files(path = data_path)

### grep .sorted.bam
input_file <- file_name[grep(pattern = "sorted.bam$", file_name)]
input_file <- paste0(data_path, input_file)

### prepare files
anno <- read.table(anno_path, fill = TRUE)
genelist <- read.table(file = genelist_path)

### sets for gene indexes
genes_index <- as.character(unique(genelist[, 1]))

nmf_result <- mclapply(X = length(genes_index), mc.cores = ncores, FUN = function(k){
  ## NMFP
  gindex <- genes_index[k]
  gname <- as.character(genelist[which(as.character(genelist[, 1]) == gindex)[1], 7])
  
  exonlist <- ExonBoundary(genelist = genelist, anno = anno, gene_index = gindex, len_read = len_read)
  chr <- exonlist$chr
  strand <- exonlist$strand
  strand <- c(strand, switch(strand, "-" = TRUE, "+" = FALSE))
  gene_start <- exonlist$gene_start
  gene_end <- exonlist$gene_end
  exon_pool <- exonlist$exon_pool
  nexon <- ncol(exon_pool)
  typesbox <- FindTypes(nexon)
  len_exon <- exon_pool[2, ] - exon_pool[1, ] + 1
  anno_trans <- exonlist$anno_trans
  
  readbin <- AssignMatrix(chr = chr, strand = strand, gene_start = gene_start, gene_end = gene_end, exon_pool = exon_pool, typesbox = typesbox, len_read = len_read, input_file = input_file)
  rownames(readbin) <- typesbox
  
  readbin <- readbin[, colSums(readbin) > 0]
  
  ## normalization
  nmatrix <- Normalization(readbin, exon_pool + gene_start - 1, len_read, typesbox)
  
  rank0 <- length(exonlist$anno_trans)
  if(if_rank_anno){
    rank <- rank0
  }else{
    rank <- 0
  }
  
  bisoform <- VoteforIsoform(nmatrix = nmatrix, rank = rank, nrun = nrun, typesbox = typesbox, exon_pool = exon_pool, len_read = len_read, code_path = code_path)
  isoform_set <- bisoform$isoform_set
  rank <- bisoform$rank
  sum_read <- rowMeans(nmatrix)
  
  ## post selection
  # class 1
  if(length(names(isoform_set)) > 10){
    preselection1 <- PreFilter(names(isoform_set), sum_read, exon_pool, len_read, FALSE, anno_trans)
  }else{
    preselection1 <- PreFilter(names(isoform_set), sum_read, exon_pool, len_read, TRUE, anno_trans)
  }

  # class 2
  index_class2 <- which(isoform_set >= 2)
  if(length(index_class2) > 10){
    preselection2 <- PreFilter(names(isoform_set[index_class2]), sum_read, exon_pool, len_read, FALSE, anno_trans)
  }else{
    preselection2 <- PreFilter(names(isoform_set[index_class2]), sum_read, exon_pool, len_read, TRUE, anno_trans)
  }

  # class 3
  index_class3 <- which(isoform_set >= 5)
  if(length(index_class3) > 10){
    preselection3 <- PreFilter(names(isoform_set[index_class3]), sum_read, exon_pool, len_read, FALSE, anno_trans)
  }else{
    preselection3 <- PreFilter(names(isoform_set[index_class3]), sum_read, exon_pool, len_read, TRUE, anno_trans)
  }

  # class 4
  index_class4 <- which(isoform_set >= 10)
  if(length(index_class4) > 10){
    preselection4 <- PreFilter(names(isoform_set[index_class4]), sum_read, exon_pool, len_read, FALSE, anno_trans)
  }else{
    preselection4 <- PreFilter(names(isoform_set[index_class4]), sum_read, exon_pool, len_read, TRUE, anno_trans)
  }

  # class 5
  index_class5 <- which(isoform_set >= 25)
  if(length(index_class5) > 10){
    preselection5 <- PreFilter(names(isoform_set[index_class5]), sum_read, exon_pool, len_read, FALSE, anno_trans)
  }else{
    preselection5 <- PreFilter(names(isoform_set[index_class5]), sum_read, exon_pool, len_read, TRUE, anno_trans)
  }
  
  preselection_all <- list(preselection1$preselection1, preselection1$preselection2, preselection1$preselection3, preselection1$preselection4,
                           preselection2$preselection1, preselection2$preselection2, preselection2$preselection3, preselection2$preselection4,
                           preselection3$preselection1, preselection3$preselection2, preselection3$preselection3, preselection3$preselection4,
                           preselection4$preselection1, preselection4$preselection2, preselection4$preselection3, preselection4$preselection4,
                           preselection5$preselection1, preselection5$preselection2, preselection5$preselection3, preselection5$preselection4)
  preslct_all <- list(preselection1$preslct1, preselection1$preslct2, preselection1$preslct3, preselection1$preslct4,
                      preselection2$preslct1, preselection2$preslct2, preselection2$preslct3, preselection2$preslct4,
                      preselection3$preslct1, preselection3$preslct2, preselection3$preslct3, preselection3$preslct4,
                      preselection4$preslct1, preselection4$preslct2, preselection4$preslct3, preselection4$preslct4,
                      preselection5$preslct1, preselection5$preslct2, preselection5$preslct3, preselection5$preslct4)
  
  len_pre <- lapply(preselection_all, function(elem){length(elem)}) %>% do.call(what = c)
  order_pre <- order(len_pre) 
  
  ## choose a filtering level
  index_cand <- which(len_pre[order_pre] <= length(anno_trans) * 10)
  if(length(index_cand) == 0){index_cand <- 1}else{index_cand <- index_cand[length(index_cand)]}
  isoform_candidates <- preselection_all[[order_pre[index_cand]]]
  matched_candidates <- preslct_all[[order_pre[index_cand]]]
  
  print(paste("Process: Gene ", k, "/", length(genes_index),
              ". The number of subexons is ", nexon, ". The number of annotated isoforms is ", length(anno_trans), ". ",
              "The number of candidates is ", length(isoform_candidates), ", among which, ", length(matched_candidates), 
              " match the annotation.", sep = ""))
  
  ## NMFP.gtf
  result_to_gtf <- WriteGTF(candidates = isoform_candidates, exon_pool = exon_pool, chr = chr, strand = strand, gene_name = gname, gene_start = gene_start)
  write.table(result_to_gtf, file = output_path, quote = F, row.names = F, col.names = F, sep = "\t", append = TRUE)
})


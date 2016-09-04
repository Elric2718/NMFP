#########################################################################################
##################################  useful functions ####################################
#########################################################################################

ExonBoundary<-function(genelist, anno, gene_index, len_read, ifvis = F){
  # Returns a list containing the boundary information of the given gene
  #  Args:
  #   genelist: the gene list made from the original annotation
  #   anno: the trimmed annotation
  #   gene_index: the index of the given gene
  #   len_read: length of the read
  #   ifvis: a logical value. Display the profile of the gene if True. Default is False
  #  Returns:
  #   chr: the chromosome of the given gene
  #   strand: the strand of the given gene
  #   gene_start: the start site of the gene
  #   gene_end: the end site of the gene
  #   exon_pool: a pool of exons in this gene
  #   anno_trans: the annotated transcripts
  
 
  ## basic information
  trans <- genelist[as.character(genelist[, 1]) == gene_index, ]
  ntrans <- nrow(trans)
  chr <- as.character(trans[1, 2])
  start <- min(as.numeric(as.character(trans[, 3])))
  end <- max(as.numeric(as.character(trans[, 4])))
  strand <- as.character(trans[1, 5])
  len <- end - start + 1
  
  ## Build the exon pool
  trans_pool <- unique(as.character(trans[, 6]))
  
  index_exon <- intersect(which(as.character(anno[, 2]) == "exon"),
                         which(is.element(as.character(anno[, 6]), trans_pool)))
  
  exon_start0 <- as.numeric(as.character(anno[index_exon, 3])) - start + 1
  exon_end0 <- as.numeric(as.character(anno[index_exon, 4])) - start + 1
  exon_pool0 <- data.frame(cbind(exon_start0,exon_end0)) %>%
                unique() %>%
                arrange(exon_start0) %>%
                as.matrix() %>%
                t()
                
  
  node <- union(exon_start0, exon_end0) %>% sort

  exon_pool <- sapply(seq(length(node) - 1), function(i){
    if(length(intersect(which(node[i] >= exon_pool0[1, ]), which(node[i+1] <= exon_pool0[2, ]))) > 0){
      c(node[i],node[i+1])
    }
  }) %>% do.call(what = cbind) %>% as.matrix(nrow = 2)
    
  
  ## merge small subexons of which the lengths are below len_read - 10  
  nexon <- ncol(exon_pool)
  len_exon <- exon_pool[2, ]-exon_pool[1, ] + 1
  
  if(nexon > 1){
    
    first_time = TRUE
    sindex = NULL
  
    while(length(sindex) > 0 |first_time){
      # Find all subexons which can be merged into the previous or the next subexon
      sindex <- intersect(which(len_exon < len_read - 10), 
                          which(sapply(seq(nexon), function(x){
        if(x == 1 && exon_pool[2, 1] >= exon_pool[1, 2] - 1){return(TRUE)}
        if(x == nexon && exon_pool[1, nexon] - 1 <= exon_pool[2, (nexon - 1)]){return(TRUE)}
        if(x > 1 && x < nexon && (exon_pool[2, x] >= exon_pool[1, (x + 1)] - 1 | exon_pool[1, x] - 1 <= exon_pool[2, (x - 1)])){return(TRUE)}
        return(FALSE)  
      }) == TRUE)) %>% sort()
      sindex0 <- sindex
      
      # Only find continuous short subexons or extremely short subexons. 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!Note: More considerations are needed here!!!!!!!!!!!!!!!!!!!!!!!#
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
      if(length(sindex) > 1){
        cont_index <- sindex[which(sindex[-1] - sindex[-length(sindex)] == 1)]
        if(length(cont_index) > 0){
          sindex <- sapply(cont_index, function(x){
            if(len_exon[x] > len_exon[x + 1]){return(x + 1)}else{return(x)}
          })}}else{
            sindex <- NULL
          }
      sindex <- c(sindex, sindex0[len_exon[sindex0] < 10]) %>% unique()
      
      # Merge short subexons
      if(length(sindex)==0){
        break
      }
      if(sindex[1] == 1){
        exon_pool[2, 1] = exon_pool[2,2]
        exon_pool[1, 2] = exon_pool[1,1]
      }
      
      if(sindex[1] == nexon){
        exon_pool[2, (nexon-1)]=exon_pool[2, nexon]
        exon_pool[1, nexon]=exon_pool[1, (nexon-1)]  
      }
      
      if(sindex[1] < nexon && sindex[1] > 1){
        if((exon_pool[2, sindex[1]] >= exon_pool[1, (sindex[1] + 1)] - 1) && (exon_pool[1, sindex[1]] - 1 > exon_pool[2, (sindex[1]-1)])){
          exon_pool[2, sindex[1]] = exon_pool[2, (sindex[1] + 1)]
          exon_pool[1, (sindex[1] + 1)] = exon_pool[1, sindex[1]]  
        }
        
        if((exon_pool[2, sindex[1]] < exon_pool[1, (sindex[1] + 1)] -1) && (exon_pool[1, sindex[1]] - 1 <= exon_pool[2, (sindex[1] - 1)])){
          exon_pool[2, (sindex[1]-1)] = exon_pool[2, sindex[1]]
          exon_pool[1, sindex[1]] = exon_pool[1, (sindex[1] - 1)]  
        }
        
        if((exon_pool[2, sindex[1]] >= exon_pool[1, (sindex[1] + 1)] - 1) && (exon_pool[1, sindex[1]] - 1 <= exon_pool[2, (sindex[1] - 1)])){
          if(len_exon[(sindex[1] - 1)] <= len_exon[(sindex[1] + 1)]){
            exon_pool[2, (sindex[1] - 1)] = exon_pool[2, sindex[1]]
            exon_pool[1, sindex[1]] = exon_pool[1, (sindex[1] - 1)] 
          }else{
            exon_pool[2, sindex[1]] = exon_pool[2, (sindex[1] + 1)]
            exon_pool[1, (sindex[1] + 1)] = exon_pool[1, sindex[1]]  
          }
        }
        
      }
      
                   
      exon_pool <- exon_pool[, duplicated(as.data.frame(t(exon_pool))) == F] %>%
                   as.matrix(nrow = 2)
      nexon <- ncol(exon_pool)
      len_exon <- exon_pool[2, ] - exon_pool[1, ] + 1
      first_time=FALSE
      
    }
    
  }
  
  
  ## visualization of the gene profile
  if(ifvis){
    exonpool <- unique(c(unlist(apply(exon_pool, 2, function(x){seq(x[1], x[2])}))))
    plot(x= exonpool, y = rep(0, length(exonpool)), xlim=c(0, exon_pool[2, nexon] + 100), ylim=c(-1, ntrans + 1), xlab="bases", ylab="transcript", pch=".")
    abline(v = union(exon_pool[1, ], exon_pool[2, ]), col = "red")
    for(i in 1:ntrans){
      par(new = TRUE)
      indexexon <- intersect(which(as.character(anno[, 2]) == "exon"), which(as.character(anno[, 6]) == trans_pool[i]))
      exon <- rbind(as.numeric(as.character(anno[indexexon, 3])) - start + 1, as.numeric(as.character(anno[indexexon, 4])) - start + 1)
      trans <- unique(c(unlist(apply(exon, 2, function(x){return(seq(x[1],x[2]))}))))
      plot(x = trans, y = rep(i, length(trans)), xlim = c(0, exon_pool[2, nexon] + 100), ylim = c(-1, ntrans + 1), xlab = "bases", ylab = "transcript", pch = ".")
    }
    par(new = FALSE)
  }
  
  
  ## get annotated transcripts
  anno_transcript <- matrix(0, nexon, ntrans)
  anno_trans <- NULL
  for(i in 1 : ntrans){
    indexexon <- intersect(which(as.character(anno[, 2]) == "exon"),which(as.character(anno[, 6]) == trans_pool[i]))
    exonstart <- as.numeric(as.character(anno[indexexon, 3])) - start + 1
    exonend <- as.numeric(as.character(anno[indexexon, 4])) - start + 1
    for(j in 1 : nexon){
        if(length(intersect(which(exonstart <= exon_pool[1, j]), which(exonend >= exon_pool[1, j] + 1))) > 0 | length(intersect(which(exonstart <= exon_pool[2, j] - 1), which(exonend >= exon_pool[2, j]))) > 0){
          anno_transcript[j, i] <- 1
        }
    }
    anno_trans <- c(anno_trans, paste(anno_transcript[, i], collapse = ""))
  }
  anno_trans <- unique(anno_trans)
  
  
  exon_list <- list(chr, strand, start, end, exon_pool, anno_trans)
  names(exon_list) <- c("chr", "strand", "gene_start", "gene_end", "exon_pool", "anno_trans")
  return(exon_list)
}


FindTypes<-function(nexon){
  # Returns a vector containing all the wanted bins
  #  Args:
  #   nexon: the number of the exons
  
  typesbox <- c(sapply(seq(nexon), function(i) paste(c(i, i), collapse = ":")),
                sapply(seq(nexon - 1), function(i) paste(c(i, i + 1), collapse = ":")))
  if(nexon >= 3){
   typesbox <- c(typesbox, sapply(seq(nexon - 2), function(i) paste(c(i, i + 2), collapse = ":")))
  }
  if(nexon >= 4){
    typesbox <- c(typesbox, sapply(seq(nexon - 3), function(i) paste(c(i, i + 3), collapse = ":")))
  }
}


AssignMatrix <- function(chr, strand, gene_start, gene_end, exon_pool, typesbox, len_read, input_file){
  # Returns a read count matrix with rows as the bins and columns as the samples
  #  Args:
  #   chr: the chromosome of the given gene
  #   strand: the strand of the given gene
  #   gene_start: the start site of the gene
  #   gene_end: the end site of the gene
  #   exon_pool: a pool of exons in this gene
  #   typesbox: a vector containing all the bins
  #   len_read: the length of the read
  #   input_file: the files referring to .bam file
  nexon <- ncol(exon_pool)
  nsample <- length(input_file)
 
  read_matrix <- matrix(0, length(typesbox), nsample)
  
  ## makes chr into the format like chr1
  chr <- unlist(strsplit(chr,"chr"))
  chr <- chr[length(chr)]
  chr <- c(chr, paste("chr", chr, sep = ""))
  chr0 <- chr
  
  ## get reads with rbamtools
  for(k in 1:nsample){
    bam <- input_file[k]
    idx <- paste(input_file[k], ".bai", sep = "")
    reader <- bamReader(bam)
    loadIndex(reader, idx)
    
    ## modify information of chr and strand
    if(sum(getRefData(reader)$SN == chr0[1]) > 0){
      chr_index <- which(getRefData(reader)$SN == chr0[1])
      chr <- chr0[1]
    }else{
      chr_index <- which(getRefData(reader)$SN == chr0[2])
      chr <- chr0[2]
    }
    chr <- c(chr, getRefData(reader)$ID[chr_index])
    
    ## extract data and process (such process can only be applied to simulation)
    coords <- c(as.numeric(chr[2]), gene_start, gene_end)
    range <- bamRange(reader, coords)
    rdf <- as.data.frame(range)
    
    if(nrow(rdf) == 0){next}
    ## extract useful information
    
    # flag information
    flag <- t(matrix(intToBits(rdf$flag), nrow = 32)[1:11, ])
    
    #exclude exons which do not coordiate whith strand
    if(strand[1] == "-"){
      sign_index <- union(
        intersect(which(flag[, 7] == 1), which(as.logical(rdf$revstrand) == TRUE)),
        intersect(which(flag[, 8] == 1), which(as.logical(rdf$revstrand) == FALSE)))
    }else{
      sign_index <- union(
        intersect(which(flag[, 7] == 1), which(as.logical(rdf$revstrand) == FALSE)),
        intersect(which(flag[, 8] == 1), which(as.logical(rdf$revstrand) == TRUE)))
    }
    
    sign_chr <- intersect(sign_index, which(rdf$refid == as.numeric(chr[2])))
    if(length(sign_chr) == 0){next}
    read_data <- rdf[sign_chr, c(2,3,4,5,8)]
    flag <- flag[sign_chr, ]
    
    ## Use Cigar information to calculate the length of seq 
    cigar_len <- sapply(seq(nrow(read_data)), function(i){
      if(read_data$nCigar[i] > 1){
        # split the letters in Cigar
        cigar <- unlist(strsplit(read_data$cigar[i], ""))
        cigar_class <- which(is.element(cigar, LETTERS))
        cigar_class <- rbind(rbind(cigar_class, c(1, cigar_class[1 : (length(cigar_class) - 1)] + 1)), cigar_class - 1)
        
        # construct a matrix indicating different parts of the Cigar
        cigar <- rbind(cigar[cigar_class[1, ]],
                     sapply(seq(ncol(cigar_class)), function(x){
                       return(as.numeric(paste(cigar[cigar_class[2, x] : cigar_class[3, x]], collapse = "")))
                     }))
        # take the sum to get the length between the end and start of the read
        len_seq <- sum(sapply(seq(ncol(cigar)), function(x){
          if(is.element(cigar[1, x], c("M", "N", "D"))){return(as.numeric(cigar[2, x]))}else{return(0)}
        }))}else{
          len_seq <- as.numeric(unlist(strsplit(read_data[i, 3], "M")))
        }
      return(len_seq)})
    
    read_data <- data.frame(start = read_data$position, end = read_data$position + cigar_len - 1)
    
    ##assign reads
    read_matrix[, k] <- ReadCount(read_data, exon_pool + gene_start - 1, typesbox, len_read)
  }
  
  return(read_matrix)
}


ReadCount <- function(read_data, exonbd, typesbox, len_read){
  # Returns a vector containing the read counts in each bin
  #  Args:
  #   read_data: a two-column matrix with each row as a read, and the first column as
  #              the start site, the second column as the end site.
  #   exonbd: a two-row matrix with each column as a subexon, and the first row as the start
  #           site, the second row as the end site.
  #   typesbox: a vector containg all the bins.
  #   len_read: the length of the read.
  nexon <- ncol(exonbd)
  nread <- nrow(read_data)
  len_exon <- exonbd[2, ] - exonbd[1, ] + 1
  sindex <- which(len_exon < len_read)
  sindex <- setdiff(sindex, c(1, nexon))
  sindex <- c(sindex, 0)  #guarantee that sindex is a vector with length over 1
  
  # a list: each element is a vector of the read indexes with the first sites in the
  #         corresponding exon
  indexlist1 <- lapply(seq(nexon), function(i){
    intersect(which(read_data$start >= exonbd[1, i] - 1),
              which(read_data$start <= exonbd[2, i] + 1))
  })
  # a list: each element is a vector of the read indexes with the end sites in the
  #         corresponding exon
  indexlist2 <- lapply(seq(nexon), function(i){
    intersect(which(read_data$end >= exonbd[1, i] - 1),
              which(read_data$end <= exonbd[2, i] + 1))
  })
  
  
  num_type <- rep(0, length(typesbox))
  for(i in 1:length(typesbox)){
    site1 <- as.numeric(as.character(strsplit(typesbox[i], ":")[[1]][1]))
    site2 <- as.numeric(as.character(strsplit(typesbox[i], ":")[[1]][2]))
    
    read_list <- intersect(indexlist1[[site1]],indexlist2[[site2]])
    
    if(length(read_list) > 0){
      special0 <- intersect(which(site1 < sindex), which(site2 > sindex))
      if(length(special0) == 0){
        num_type[i] <- num_type[i] + length(read_list)
      }else if(length(special0) > 0){
        for(k in 1:length(special0)){
          special <- sindex[special0[k]]
          skip_list <- read_list[which(sapply(read_list, function(x){
            return(exonbd[2, site1] - read_data[x, 1] + 1 + read_data[x,2] - exonbd[1, site2] + 1)}) >= len_read)]
          retain_list <- setdiff(read_list, skip_list)
          
          num_type[i] <- num_type[i] + length(skip_list)
          
          num_type[which(typesbox == paste(site1, special, sep = ":"))] <- 
            num_type[which(typesbox == paste(site1, special, sep = ":"))] + length(retain_list)
          num_type[which(typesbox == paste(special, special, sep = ":"))] <- 
            num_type[which(typesbox == paste(special, special, sep = ":"))] + length(retain_list)
          num_type[which(typesbox == paste(special, site2, sep = ":"))] <- 
            num_type[which(typesbox == paste(special, site2, sep = ":"))] + length(retain_list)
        }
      }
    }
  }
  return(num_type)
}


Normalization <- function(read_matrix, exon_pool, len_read, typesbox){
  # Returns a normalized read-count matrix
  #  Args:
  #   read_matrix: the original read-count matrix
  #   exon_pool: a pool of exons in this gene
  #   typesbox: a vector containing all the bins
  #   len_read: the length of the read
  nsample <- ncol(read_matrix)
  nexon <- ncol(exon_pool)
  ntype <- length(typesbox)
  
  ## per Ave reads
  ave_read <- mean(colSums(read_matrix))
  nmatrix <- apply(read_matrix, 2, function(x){x / sum(x)}) * ave_read
  
  ## bin length
  len_bin <- sapply(seq(ntype), function(i){PairReadLength(exon_pool,len_read,typesbox[i])})
  
  ## per 10000 bases
  nmatrix <- nmatrix / len_bin
  
  return(nmatrix)
}


PairReadLength <- function(exonbd, len_read, type){
  # Returns the effective length of the bin
  #  Args:
  #   exonbd: a two-row matrix with each column as a subexon, and the first row as the start
  #           site, the second row as the end site.
  #   len_read: the length of the read.
  #   type: the type of the inputed bin.
  len_exon <- exonbd[2, ] - exonbd[1, ] + 1
  
  if(length(strsplit(type, ":")[[1]]) == 2){
    site1 <- as.numeric(as.character(strsplit(type, ":")[[1]][1]))
    site2 <- as.numeric(as.character(strsplit(type, ":")[[1]][2]))
    
    if(site1 == site2){
      if(len_read - 1 <= len_exon[site1]){
        len <- len_exon[site1] - len_read + 1}else{
          len <- len_read - len_exon[site1] - 1}
    }else{
      if(len_read - 1<= min(len_exon[site1], len_exon[site2])){
        len <- len_read - 1
      }else{
        len <- min(len_exon[site1], len_exon[site2])
      }
    }
    
    len <- max(len, 1) / 100
  }
  return(len)
}


VoteforIsoform <- function(nmatrix, rank = 0, nrun = 100, typesbox, exon_pool, alpha0 = c(1, 10, 0.1), len_read = 76, gcutoff = c(0.4, 0.1, 0.01), code_path){
  # Returns the result of NMFP
  #  Args:
  #   nmatrix: the normalized read-count matrix
  #   typesbox: a vector containing all the bins
  #   exon_pool: a pool of exons in this gene
  #   rank: the NMF rank. If it's a vector or 0, then the optimal one will be tuned out
  #   nrun: the number of runs of NMFP
  #   alpha0 : the coefficient of L1 penalty. If it's a vector, then the optimal one will be tuned out
  #   len_read: the length of the read
  #   gcutoff: the cutoff to binarize the NMF basis matrix. If it's a vector, then the optimal one will be tuned out
  #   code_path: the path containing the NMF code implemented by C++.
  nexon <- ncol(exon_pool)
  index_nzero <- which(rowMeans(nmatrix) > 0 * rowMeans(nmatrix))
  
  # Contradict Table
  ctable <- ContradictTable(typesbox)
  cset <- ctable[index_nzero, index_nzero]
  
  # Tune parameters
  par <- TunePar(nmatrix = nmatrix, rank = rank, alpha0 = alpha0, gcutoff = gcutoff, cset = cset, index_nzero = index_nzero, exon_pool = exon_pool, len_read = len_read, typesbox = typesbox, code_path = code_path)
  
  # multi runs
  bisoform <- replicate(nrun, expr = {
    ncNMFresult <- cmp_ncNMF2(V = nmatrix[index_nzero, ], alpha0 = par$alpha0, rank = par$rank, cset=cset, code_path = code_path)
    basis <- ncNMFresult$basis
    coef <- ncNMFresult$coefficient
    
    # Affinity Matrix
    affinity_matrix <- diag(apply(basis, 1, max))
    
    # G Matrix 
    G <- solve(affinity_matrix) %*% basis
    
    # trimmed G
    G <- apply(G, 2, function(x){y <- rep(0, length(x)); y[which(x >= gcutoff)] <- 1; return(y)})
    Gmatrix <- apply(G, 2, function(x){y <- rep(0, dim(nmatrix)[1]); y[index_nzero[which(x == 1)]] <- 1; return(y)})
    
    # visualize Gmatrix
    Gvisualization <- VisualizeGmatrix(Gmatrix, If_print = FALSE, typesbox = typesbox, exon_pool = exon_pool, len_read = len_read)
    
    Gvisual <- Gvisualization$certain_vec
    Gvisual[which(Gvisual == -1)] <- 0
    Gvisual[which(Gvisual > 0)] <- 1
    Gvisual <- rbind(Gvisual, ":")
    Gvisual <- as.vector(Gvisual)
    
    #represented as binary and then transform into decimal
    strsplit(paste(Gvisual, collapse = ""), ":")[[1]]
  }) %>% do.call(what = c)
  
  #choose the top isoforms
  isoform_set <- sort(table(bisoform), decreasing = TRUE)
  isoform <- list(isoform_set, alpha0, rank, gcutoff)
  names(isoform) <- c("isoform_set", "alpha0", "rank", "gcutoff")
  return(isoform)
}


ContradictTable <- function(typesbox){
  # Returns a table indicating all the contradicting bins. Bin i and Bin j
  # contradicts with each other if (i, j) = 1.
  #  Args:
  #   typesbox: a vector containing all the bins
  ntype <- length(typesbox)
  
  ctable <- lapply(seq(ntype), function(i){
    sapply(seq(ntype), function(j){
      FindConflict(typesbox[i], typesbox[j])
    })
  }) %>% do.call(what = cbind)
  
  return(ctable)
}


FindConflict <- function(type1, type2){
  # Returns if Bin 1 conflicts with Bin2
  #  Args:
  #   type1, 2: Bin type of Bin 1, and Bin 2
  
  site1 <- as.numeric(as.character(strsplit(type1, ":")[[1]]))
  site2 <- as.numeric(as.character(strsplit(type2, ":")[[1]]))

  diff11 <- site1[1] - site2[1]
  diff21 <- site1[2] - site2[1]
  diff12 <- site1[1] - site2[2]
  diff22 <- site1[2] - site2[2]
  
  if((diff11 < 0 && diff21 > 0) | (diff12 < 0 && diff22 > 0) | (diff11 > 0 && diff12 < 0) | (diff21 > 0 && diff22 < 0)){
    return(1)
  }else{
    return(0)
  }
}


ncNMF2 <- function(V, rank, alpha0 = 10, cset, iter = 300, eps = 0.0005, code_path){
  # Returns the NMF result
  #  Args:
  #   V: the matrix NMF to be applied on
  #   rank: the NMF rank
  #   alpha0: the coefficient of L1 penalty
  #   cset: a table containing all the pairs of conflicting bins
  #   iter: the maximal number of iterations
  #   eps: the convergence accuracy
  #   code_path: the path to NMF codes implemented in C++
  
  ## load C based function
  #load Update
  dyn.load(paste(code_path,"/Update.so",sep=""))
  dyn.load(paste(code_path,"/CpenaltyObj.so",sep=""))
  
  ##initial the matrix
  #randomly set intial matrix
  nrow <- dim(V)[1]
  ncol <- dim(V)[2]
  meanV <- sum(V)/(nrow*ncol)
  W <- abs(matrix(rnorm(nrow * rank, mean = 1 / nrow, sd = sqrt(1 / nrow)), nrow = nrow, ncol = rank))
  W <- apply(W, 2, function(x){return(x / sum(x))})
  H <- abs(matrix(rnorm(rank * ncol, mean = meanV * nrow, sd = sqrt(meanV * nrow)), nrow = rank, ncol = ncol))
  
  ## name W and H
  rownames(W) <- 1 : nrow
  colnames(W) <- LETTERS[1 : rank]
  
  rownames(H) <- LETTERS[1 : rank]
  colnames(H) <- 1 : ncol
  
  # Code dependent alpha
  alpha <- mean(W %*% H) / mean(W %*% t(W)) * alpha0
  
  # initial error
  err0 <- .Call("Cpenalty", W, H, V, cset, alpha)
  err <- err0 / 2
  count <- 0
  
  while(abs(err - err0) / err0 >= eps & count < iter){
    alpha <- mean(W %*% H) / mean(W %*% t(W)) * alpha0

    count <- count + 1
    
    WH <- .Call("Update", W, H, V, cset, alpha)
    W <- WH[1 : nrow,]
    H <- t(WH[(nrow + 1) : (nrow + ncol), ])
    
    W <- t(t(W))
    rownames(W) <- 1 : nrow
    colnames(W) <- LETTERS[1 : rank]
    
    rownames(H) <- LETTERS[1 : rank]
    colnames(H) <- 1 : ncol
    
    if(count %% 10 == 0){
      err0 <- err
      err <- .Call("Cpenalty", W, H, V, cset, alpha)
    }
  }
  
  #unload Update  
  dyn.unload(paste(code_path, "/Update.so", sep = ""))  
  dyn.unload(paste(code_path, "/CpenaltyObj.so", sep = "")) 
  
  ncNMFresult <- list(basis = W, coefficient = H, error = err)
  return(ncNMFresult)
}


TunePar <- function(nmatrix, rank, alpha0, gcutoff, cset, index_nzero, ...){
  # Returns the optimal parameter set including the rank, the alpha0, the gcutoff
  #  Args:
  #   V: the matrix NMF to be applied on
  #   rank: the NMF rank to be tuned if inputing a vector or 0
  #   alpha0: the coefficeint of L1 penalty in NMF to be tuned if inputing a vetor
  #   gcutoff: the cutoff to binarize the basis matrix to be tuned if inputing a vctor
  #   cset: a table containing all the pairs of conflicting bins
  #   index_nzero: the indexes for non-zero bins
  #   ...: arguments passed from VoteforIsoform
  optional_par <- list(...)
  exon_pool <- optional_par$exon_pool
  len_read <- optional_par$len_read
  code_path <- optional_par$code_path
  typesbox <- optional_par$typesbox
  
  V <- nmatrix[index_nzero, ]
  ## determine ranks
  if(length(rank) == 1 && rank == 0){
    rank <- RankDetermine(V = V, cset = cset,
                          alpha0 = 1, rank_start = 2, rank_end = 5,  exon_pool = exon_pool, len_read = len_read, code_path = code_path)
  }else{
    rank <- rank[rank > 0]
  }
  
  parameter <- expand.grid(rank = rank, alpha0 = alpha0, gcutoff = gcutoff) 
  
  ## tuning parameters of alpha0 and gcutoff
  uncertainty <- apply(parameter, 1, function(par){
    replicate(10, expr = {
      ncNMFresult <- cmp_ncNMF2(V = V, rank = par[1], alpha0 = par[2], cset = cset, code_path = code_path)
      basis <- ncNMFresult$basis
      coef <- ncNMFresult$coefficient
      
      # Affinity Matrix
      affinity_matrix <- diag(apply(basis, 1, max))
      
      # G Matrix
      G <- solve(affinity_matrix) %*% basis
      
      # trimmed G
      G <- apply(G, 2, function(x){y <- rep(0, length(x)); y[which(x >= par[3])] <- 1; return(y)})
      Gmatrix <- apply(G, 2, function(x){y <- rep(0, dim(nmatrix)[1]); y[index_nzero[which(x == 1)]] <- 1; return(y)})
      
      # visualize Gmatrix
      Gvisualization <- VisualizeGmatrix(Gmatrix, FALSE, typesbox = typesbox, exon_pool = exon_pool, len_read = len_read)
      return(Gvisualization$uncertainty)
    }) %>% mean
  })
  
  index_tune <- which.min(uncertainty)
  
  return(list(rank = parameter$rank[index_tune], alpha0 = parameter$alpha0[index_tune], 
              gcutoff = parameter$gcutoff[index_tune]))
}


RankDetermine <- function(V, cset, alpha0, ntest = 20, nave = 10, rank_start, rank_end, ...){
  # Returns the rank determined by gap stat
  #  Args: 
  #   V: the matrix NMF to be applied on
  #   cset: a table containing all the pairs of conflicting bins
  #   ntest: the number of tests for gap stat
  #   nave: the number of times to calculate the NMF error for average
  #   rank_start: the starting value for the rank to be determined
  #   rank_end: the ending value for the rank to be determined
  #   ...: arguments passed from TunePar
  optional_par <- list(...)
  exon_pool <- optional_par$exon_pool
  len_read <- optional_par$len_read
  code_path <- optional_par$code_path
  
  nrow <- dim(V)[1]
  ncol <- dim(V)[2]
  nexon <- length(exon_pool[1, ])
  len_exon <- exon_pool[2, ] - exon_pool[1, ] + 1
  
  Wk <- NULL
  Wk_star <- matrix(0, ntest, (rank_end - rank_start + 1))
  gapk <- NULL
  ek <- NULL
  sdk <- NULL
  
  read_range <- cbind(apply(V, 1, min), apply(V, 1, max))
  
  for(i in rank_start : rank_end){
    # Real Error
    error <- replicate(nave, expr = {
      ncNMFresult <- cmp_ncNMF2(V = V, alpha0 = alpha0, rank = i, cset = cset, code_path = code_path)
      ncNMFresult$err
    }) %>% mean
    Wk <- c(Wk, error)
    
    index_col <- i - rank_start + 1

    Wk_star[, index_col] <- replicate(ntest, expr = {
      Vk <- lapply(seq(nrow), function(j){return(runif(ncol, min = read_range[j, 1], max = read_range[j, 2]))}) %>%
        do.call(what = rbind)
      
      replicate(nave, expr = {
        ncNMFresult <- cmp_ncNMF2(V = Vk, alpha0 = alpha0, rank = i, cset = cset, code_path = code_path)
        ncNMFresult$err
      }) %>% mean
    })
      
    ek <- c(ek, 1 / ntest * sum(log(Wk_star[, (index_col)])))
    gapk <- c(gapk, ek[(index_col)] - log(Wk[(index_col)]))
    sdk <- c(sdk, sqrt(1 + 1 / ntest) * sqrt(1 / ntest * sum((log(Wk_star[, (index_col)]) - ek[(index_col)])^2)))
    
    if(i > rank_start && (gapk[index_col - 1] >= gapk[index_col] - sdk[index_col])){
      break
    }
  }
  rank <- i
  
  return(rank)  
}


VisualizeGmatrix <- function(Gsignal, if_print = FALSE, ...){
  # Returns a list after splitting the conflicitng bins in isoforms
  #  Args:
  #   Gsignal: a matrix indicating if Bin i exists in Isoform j
  #   if_print: a logical value. Print the bin weights in each isoform if True. Default is False
  #   ...: arguments passed from VoteforIsoform
  #  Returns:
  #   certain_vec: all the isoform vectors after splitting the isoforms with conflicting bins
  #   uncertainty: the average of conflicting subexons
  optional_par <- list(...)
  exon_pool <- optional_par$exon_pool
  len_read <- optional_par$len_read
  typesbox <- optional_par$typesbox
  
  nexon <- ncol(exon_pool)
  nisoform <- ncol(Gsignal)
  ntype <- length(typesbox)
  
  len_exon <- exon_pool[2, ] - exon_pool[1, ] + 1
  sindex <- which(len_exon < len_read)
  sindex <- c(sindex,0) # guarantee it contains at least one element
  
  visual_mat <- lapply(seq(nisoform), function(j){
    lapply(seq(ntype), function(i){
      visual_result <- VisualCheck(exon_pool, len_read, typesbox[i], Gsignal[i, j])
      c(visual_result$Positive1, visual_result$Negative1, visual_result$Positive2, visual_result$Negative2)
    }) %>% do.call(what = cbind) %>% rowSums()
  }) %>% do.call(what = cbind)
  
  # weights matrix
  Vp1 <- visual_mat[1 : nexon, ]
  Vn1 <- visual_mat[(nexon + 1) : (2 * nexon), ]
  # non-weights matrix
  Vp2 <- visual_mat[(2 * nexon + 1) : (3 * nexon), ]
  Vn2 <- visual_mat[(3 * nexon + 1) : (4 * nexon), ]
  
  # trimmed Gmatrix
  Gtrimmed <- lapply(seq(nisoform), function(i){
    lapply(seq(nexon), function(j){
      if(Vn2[j, i] < 0 && Vp2[j, i] == 0){return(-1)}
      if(Vn2[j, i] == 0 && Vp2[j, i] > 0){return(1)}
      if(Vn2[j, i] < 0 && Vp2[j, i] > 0){return(0)}
    }) %>% do.call(what = rbind)
  }) %>% do.call(what = cbind)
  
  uncertainty <- apply(Gtrimmed, 2, function(vec){
    length(vec == 0)
  })
  
  # separate uncertain isoforms
  certain_vec <- lapply(seq(nisoform), function(i){
    uindex <- which(Gtrimmed[, i] == 0)
    if(length(uindex) == 0){
      Gtrimmed[, i]
    }else{
      v <- matrix(Gtrimmed[, i], nrow = nexon, ncol = 2^(length(uindex)))
      v[uindex, ] <- matrix(as.numeric(intToBits(seq(0, 2^(length(uindex)) - 1))), nrow = 32)[1 : length(uindex), ]
      v[v == 0] <- -1
      return(v)
    }
  }) %>% do.call(what = cbind)
  
  if(if_print == TRUE){
    ## indicator vector
    # weights indicator vector for existence isoform
    Ip1 <- rep(0, nexon)
    # weights indicator vector for non-existence isoform
    In1 <- rep(0, nexon)
    # non-weights indicator vector for existence isoform
    Ip2<-rep(0,nexon)
    # non-weights indicator vector for non-existence isoform
    In2<-rep(0,nexon)
    
    for(i in 1:nexon){
      if(i >= 4 && i <= nexon - 3){
        Ip1[i] <- 8
        In1[i] <- -8
        
        Ip2[i] <- 7
        In2[i] <- -4
      }
      if(i == 3 | i == nexon - 2){
        Ip1[i] <- 7
        In1[i] <- -8
        
        Ip2[i] <- 6
        In2[i] <- -4
      }
      if(i == 2 | i == nexon - 1){
        Ip1[i] <- 6
        In1[i] <- -6
        
        Ip2[i] <- 5
        In2[i] <- -3
      }
      if(i == 1 | i == nexon){
        if(is.element(i, sindex) == FALSE){
          Ip1[i] <- 5
          In1[i] <- -2.5
          
          Ip2[i] <- 4
          In2[i] <- -1
        }else{
          Ip1[i] <- 3
          In1[i] <- -5/3
          
          Ip2[i] <- 3
          In2[i] <- -3
        }
      }
    }
    
    print(Vp1)
    print(Vn1)
    print(Vp2)
    print(Vn2)
    print(Ip1)
    print(In1)
    print(Ip2)
    print(In2)
  }
  
  result <- list(certain_vec, mean(uncertainty))
  names(result) <- c("certain_vec", "uncertainty")
  return(result)
}


VisualCheck <- function(exon_pool, len_read, type, signal){
  # Record votes with consideration about balance of exisiting and non-exisiting 
  # of one exon, the positoin of the exon and short exon. Only consider long isoform
  # with the number of exons no less than 7 and the largest length of bins is 4.
  #  Args:
  #   exon_pool: a pool of exons in this gene
  #   len_read: the length of the read
  #   type: the type of the bin
  #   signal: a binary value. It indicates the existence of the bin if it's 1, otherwise
  #           implies its non-existence.
  nexon <- ncol(exon_pool)
  len_exon <- exon_pool[2, ] - exon_pool[1, ] + 1
  sindex <- which(len_exon < len_read)
  sindex <- c(sindex, 0)#guarantee it contains at least one element
  # weights vector
  VseqP1 <- rep(0, nexon)
  VseqN1 <- rep(0, nexon)
  
  # non-weights vector
  VseqP2 <- rep(0, nexon)
  VseqN2 <- rep(0, nexon)
  
  site1 <- as.numeric(as.character(strsplit(type, ":")[[1]][1]))
  site2 <- as.numeric(as.character(strsplit(type, ":")[[1]][2]))
  
  if_norm <- (is.element(site1, sindex) == FALSE && is.element(site2, sindex) == FALSE) | 
    (is.element(site1, sindex) == TRUE && site1 > 1 && is.element(site2, sindex) == FALSE)|
      (is.element(site2, sindex) == TRUE && site2 < nexon && is.element(site1, sindex) == FALSE)|
        (is.element(site1, sindex) == TRUE && is.element(site2, sindex) == TRUE && site1 > 1 && site2 < nexon)
  
  if(signal == 1){
    VseqP1[site1] <- VseqP1[site1] + 1
    VseqP1[site2] <- VseqP1[site2] + 1
    
    VseqP2[site1] <- 1
    VseqP2[site2] <- 1
    
    if(site1 == site2 - 2){
      VseqN1[site1 + 1] <- VseqN1[site1 + 1] - 2
      VseqN2[site1 + 1] <- -1
    }
    
    if(site1 == site2 - 3){
      VseqN1[site1 + 1] <- VseqN1[site1+1] - 2
      VseqN2[site1 + 1] <- -1
      VseqN1[site1 + 2] <- VseqN1[site1 + 2] - 2
      VseqN2[site1 + 2] <- -1
    }
    
    if(site1 == site2 && is.element(site1, sindex)){
      VseqP1[site1 - 1] <- VseqP1[site1 - 1] + 1
      VseqP1[site1 + 1] <- VseqP1[site1 + 1] + 1
      
      VseqP2[site1 - 1] <- 1
      VseqP2[site1 + 1] <- 1
    }
  }else{
    if(if_norm == TRUE){
      if(site1 == site2){
        VseqN1[site1] <- VseqN1[site1] - 2.5
        VseqN2[site1] <- -1
      }
    }else{
      if(is.element(site1, sindex) == TRUE && site1 == 1 && ((is.element(site2, sindex) == FALSE)|
                                                        (is.element(site2, sindex) == TRUE && site2 > 1 && site2 < nexon))){
        if(site1 + 1 == site2){
          VseqN1[site1] <- VseqN1[site1] - 1
          
          VseqN2[site1] <- -1
        }
        
        if(site1 + 1 < site2){
          VseqN1[site1] <- VseqN1[site1] -1 / 3
          
          VseqN2[site1] <- -1
        }
      }
      
      if(is.element(site2, sindex) == TRUE && site2 == nexon && ((is.element(site1, sindex) == FALSE)|
                                                            (is.element(site1, sindex) == TRUE && site1 > 1 && site1 < nexon))){
        if(site1 + 1 == site2){
          VseqN1[site2] <- VseqN1[site2] - 1
          
          VseqN2[site2] <- -1
        }
        
        if(site1 + 1 < site2){
          VseqN1[site2] <- VseqN1[site2] - 1 / 3
          
          VseqN2[site2] <- -1
        }
      }
      
      if((is.element(site1, sindex) == TRUE) && (is.element(site2, sindex) == TRUE) && site1 == 1 && site2 == nexon){
        if(site1 + 1 == site2){
          VseqN1[site1] <- VseqN1[site1] - 1
          VseqN1[site2] <- VseqN1[site2] - 1
          
          VseqN2[site1] <- -1
          VseqN2[site2] <- -1
        }
        
        if(site1 + 1 < site2){
          VseqN1[site1] <- VseqN1[site1] - 1 / 3
          VseqN1[site2] <- VseqN1[site2] - 1 / 3
          
          VseqN2[site1] <- -1
          VseqN2[site2] <- -1
        }
      }
    }
  }
  Vseq <- list(VseqP1, VseqN1, VseqP2, VseqN2)
  names(Vseq) <- c("Positive1","Negative1","Positive2","Negative2")
  return(Vseq)
}


PreSelection <- function(candidates = NULL, sum_read, exon_pool, len_read, if_part = FALSE, mode = 4){
  # Returns the pre-filtered isoform candidates which are restained if there are reads 
  # supporting them.
  #  Args:
  #   candidates: a vector of isoforms to be filtered like "11110000"
  #   sum_read: the sum of reads in all the samples for each bin
  #   exon_pool: a pool of exons in this gene
  #   len_read: the length of the read
  #   if_part: a logical value. Filter the isoforms in terms of the supporting reads if False.
  #           The default is False.
  #   mode: the level of filtering.
  
  nexon <- ncol(exon_pool)
  len_exon <- exon_pool[2, ] - exon_pool[1, ] + 1
  sindex <- union(intersect(which(len_exon <= len_read + 10), which(len_exon >= len_read - 10)), 
                  which(len_exon <= 20))
  reference <- sum_read
  reference[reference > 0] <- 1
  ntype <- length(sum_read)
  candidates0 <- candidates
  
  if(if_part == FALSE){
    if(length(candidates) == 0){
      candidates <- matrix(as.numeric(intToBits(seq(1, 2^(nexon) - 1))), nrow = 32)[1 : nexon, ] 
    }else{
      candidates <- setdiff(candidates, paste(rep(0, nexon), collapse=""))
      candidates <- paste(candidates, collapse = "")
      candidates <- as.numeric(strsplit(candidates, "")[[1]])
      candidates <- matrix(candidates, nrow = nexon)
    }
    
    ncand <- ncol(candidates)
    
    gindex <- which(apply(candidates, 2, function(x){
      #stage 1
      match1 <- prod(reference[setdiff(which(x == 1), c(sindex, 1, nexon))])
      if(len_exon[1] >= len_read + 5 && x[1] == 1){
        match1 <- match1 * reference[1]
      }
      if(len_exon[nexon] >= len_read + 5 && x[nexon] == 1){
        match1 <- match1 * reference[nexon]
      }
      
      #stage 2
      match2 <- 0
      if(match1 == 1){
        match2 <- prod(reference[nexon + which(
          sapply(seq(nexon-1), function(i){if(x[i] * x[i + 1] == 1){return(TRUE)}else{return(FALSE)}}))])
      }
      
      #stage 3
      match3 <- 0
      if(match2 == 1){
        match3 <- prod(reference[2 * nexon - 1 + which(
          sapply(seq(nexon - 2), function(i){if(x[i] * x[i+2] == 1 && x[i + 1] == 0){return(TRUE)}else{return(FALSE)}}))])
      }
      
      #stage 4
      if(nexon >= 4){
        match4 <- 0
        if(match3 == 1){
          match4 <- prod(reference[3 * nexon - 3 + which(
            sapply(seq(nexon-3), function(i){if(x[i] * x[i + 3] == 1 && x[i + 1] == 0 && x[i + 2] == 0){return(TRUE)}else{return(FALSE)}}))])
        }
      }else{match4 <- 1}
      
      
      
      if(switch(mode, match1, match2, match3, match4) == 1){
        return(TRUE)
      }else{return(FALSE)}
      
    }))
    
    
    if(length(gindex) > 0){
      select_cands <- candidates[, gindex]
      if(length(gindex) == 1){
        select_cands <- t(t(select_cands))
      }
      select_cands <- as.vector(rbind(select_cands, ":"))  
      select_cands <- strsplit(paste(select_cands, collapse=""), ":")[[1]]
    }else{
      select_cands <- setdiff(candidates0, paste(rep(0, nexon), collapse = ""))
    }
    
  }else{
    select_cands <- setdiff(candidates0, paste(rep(0, nexon), collapse = ""))
  }
  return(select_cands)
}


PreFilter <- function(candidates, sum_read, exon_pool, len_read, if_part = FALSE, anno_trans){
  # Returns a list containing four cases of filtered isoforms
  #  Args:
  #   candidates: a vector of isoform candidates from NMFP
  #   sum_read: the sum of reads in all the samples for each bin
  #   exon_pool: a pool of exons in this gene
  #   len_read: the length of the read
  #   if_part: a logical value. Filter the isoforms in terms of the supporting reads if False.
  #           The default is False.
  #   anno_trans: the annotated transcripts
  
  
  preselection1 <- PreSelection(candidates = candidates, sum_read=sum_read, exon_pool = exon_pool, len_read = len_read, if_part = if_part, mode = 1)
  preselection2 <- PreSelection(candidates = candidates, sum_read=sum_read, exon_pool = exon_pool, len_read = len_read, if_part = if_part, mode = 2)
  preselection3 <- PreSelection(candidates = candidates, sum_read=sum_read, exon_pool = exon_pool, len_read = len_read, if_part = if_part, mode = 3)
  preselection4 <- PreSelection(candidates = candidates, sum_read=sum_read, exon_pool = exon_pool, len_read = len_read, if_part = if_part, mode = 4)
  
  preslct1 <- sapply(seq(length(anno_trans)), function(i){preselection1[which(is.element(preselection1, anno_trans[i]))]}) %>% unlist
  preslct2 <- sapply(seq(length(anno_trans)), function(i){preselection2[which(is.element(preselection2, anno_trans[i]))]}) %>% unlist
  preslct3 <- sapply(seq(length(anno_trans)), function(i){preselection3[which(is.element(preselection3, anno_trans[i]))]}) %>% unlist
  preslct4 <- sapply(seq(length(anno_trans)), function(i){preselection4[which(is.element(preselection4, anno_trans[i]))]}) %>% unlist
  
  
  list(preselection1 = preselection1,
       preselection2 = preselection2,
       preselection3 = preselection3,
       preselection4 = preselection4,
       preslct1 = preslct1,
       preslct2 = preslct2,
       preslct3 = preslct3,
       preslct4 = preslct4)
}


WriteGTF <- function(candidates, exon_pool, chr, strand, gene_name, gene_start){
  # Returns a data frame containing the NMFP result in the format of gtf
  #  Args:
  #   candidates: a vector of isoform candidates from NMFP
  #   exon_pool: a pool of exons in this gene
  #   chr: the chromosome of the given gene
  #   strand: the strand of the given gene
  #   gene_name: the name of the given gene
  #   gene_start: the start site of the gene
  result_to_gtf <- NULL
  for(i in 1 : length(candidates)){
    boundary <- as.matrix(exon_pool[, which(unlist(strsplit(candidates[i], split = "")) == 1)], nrow = 2)
    if(ncol(boundary) > 1){
      j <- 1
      j2 <- 1
      while(j <= (ncol(boundary) - 1)){
        if(boundary[2, j] >= boundary[1, j + 1] - 1){
          j2 <- j + 1
          while(j2 <= (ncol(boundary) - 1) && boundary[2, j2] >= boundary[1, j2 + 1] - 1){
            j2 <- j2 + 1
          }
          for(t in j : j2){
            boundary[2, t] <- boundary[2, j2]
            boundary[1, t] <- boundary[1, j]
          }
        }
        j <- j2 + 1
        j2 <- j
      }
      boundary <- as.matrix(boundary, nrow = 2)
      boundary <- boundary[, which(duplicated(t(boundary)) == FALSE)]
      boundary <- as.matrix(boundary, nrow = 2)
    }
    boundary <- boundary + gene_start - 1
    
    tx_gtf <- data.frame(
      chrom = chr,
      Source = "NMFP",
      feature = c("transcript", rep("exon", ncol(boundary))),
      start = c(min(boundary), boundary[1, ]),
      end = c(max(boundary), boundary[2, ]),
      score = ".",
      strand = strand[1],
      frame = ".",
      attributes = paste('gene_id \"', gene_name, '\"; ', 'transcript_id \"', gene_name, "_", i, '\";', sep = "")
    )
    result_to_gtf <- rbind(result_to_gtf, tx_gtf)
  }
  
  return(result_to_gtf)
}



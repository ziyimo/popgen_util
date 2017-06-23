calc_dis <- function(X, Y, diploid = TRUE){
  # This function calculates the pairwise genetic distance (a derivative of Gronau et al. 2011)
  # For haplotypes, this is simply the observed no. of differences between the pair (i.e. hamming distance)
  # The input are two vectors of genotypes (X, Y)
  # The (unphased) genotype should conform to the standard VCF specification (0/0, 0/1 etc.)
  
  dis_vec <- numeric(length(X))
  
  if (diploid){
    for (i in 1:length(X)){
      if (is.na(X[i]) | is.na(Y[i])){
        dis_vec[i] <- NA
        next
      }
      gtX <- unlist(strsplit(X[i],'/'))
      gtY <- unlist(strsplit(Y[i],'/'))
      if (gtX[1]==gtY[1] & gtX[2]==gtY[2]){
        dis_vec[i] <- 0
      }else if (gtX[1]==gtY[1] | gtX[2]==gtY[2]){
        dis_vec[i] <- 1/2
      }else{
        dis_vec[i] <- 1
      }
    }
  }else{
    for (i in 1:length(X)){
      if (is.na(X[i]) | is.na(Y[i])){
        dis_vec[i] <- NA
        next
      }
      gtX <- X[i]
      gtY <- Y[i]
      if (gtX==gtY){
        dis_vec[i] <- 0
      }else{
        dis_vec[i] <- 1
      }
    }
  }
  
  dis <- sum(dis_vec, na.rm = TRUE) # for now, the NA are discarded from distance calculation
  
  return(dis)
}

pairwise_dis <- function(gt_Mtx, diploid = TRUE){
  # Given a matrix of genotypes extracted from VCF file (col: sample, row: SNP),
  # this function outputs all pairwise distance (observed # of differences)
  # The output is a list with the following two objects
  # 1) dataframe with pairwise distance in each row
  # 2) distance matrix
  
  dimension <- ncol(gt_Mtx)
  
  binNK <- choose(dimension, 2) # pre-allocate space for better efficency
  samp1 <- character(binNK)
  samp2 <- character(binNK)
  dis <- numeric(binNK)

  pwDis_mtx <- matrix(nrow = dimension, ncol = dimension)
  colnames(pwDis_mtx) <- colnames(gt_Mtx)
  rownames(pwDis_mtx) <- colnames(gt_Mtx)
  
  counter <- 0
  for (i in 1:(dimension-1)){
    for (j in (i+1):dimension){
      counter <- counter+1
      distance <- calc_dis(gt_Mtx[,i], gt_Mtx[,j], diploid)
      samp1[counter] <- colnames(gt_Mtx)[i]
      samp2[counter] <- colnames(gt_Mtx)[j]
      dis[counter] <- distance
      pwDis_mtx[i,j] <- distance
    }
  }
  
  upperTri <- pwDis_mtx[upper.tri(pwDis_mtx)]
  pwDis_mtx <- t(pwDis_mtx)
  pwDis_mtx[upper.tri(pwDis_mtx)] <- upperTri
  diag(pwDis_mtx) <- 0
  
  pwDis_df <- data.frame(samp1, samp2, dis)
  
  pwDis <- list(pwDis_df, pwDis_mtx)
  return(pwDis)
}

collapseGT <- function(vcf_gt){
  # collapses '0/0' to '0', '1/1' to '1'
  gt_vec <- unlist(strsplit(vcf_gt,'/'))
  if(gt_vec[1] != gt_vec[2]){
    print("CANNOT collapse heterozygous genotype to haplotype!")
    return(NULL)
  }
  return(gt_vec[1])
}


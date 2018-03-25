#!/usr/bin/env Rscript

collapseGT <- function(vcf_gt){
  # collapses '0/0' to '0', '1/1' to '1'
  # randomly collapses '0/1' to '0' or '1'

  if (is.na(vcf_gt)){
    return(NA)
  }
  gt_vec <- unlist(strsplit(vcf_gt,'/'))
  if (gt_vec[1] == gt_vec[2]){
    return(gt_vec[1])
  } else {
    return(sample(gt_vec, 1))
  }
}

vcf2Phylip <- function(mode, inPath, outPath){
  # Use 'inPath' to specify a vcf file
  # This function outputs a phylip file at 'outPath' 
  # The way to produce haplotypes from heterozygous site depends on 'mode'
  # if mode == "-r", the heterozygous sites are randomly sampled
  # if mode == "-x", the heterozygous sites are excluded
  
  library(vcfR)
  library(ape)
  library(phangorn)
  
  vcfFile <- read.vcfR(inPath)
  
  gt <- extract.gt(vcfFile, return.alleles = TRUE)
  tot_cnt <- nrow(gt)
  # apply all filters
  flt <- getFILTER(vcfFile)
  passAll <- which(flt=="PASS")
  gt <- gt[passAll, ]
  pass_cnt <- length(passAll)
  cat(pass_cnt, "sites out of", tot_cnt, "sites pass all filters applied.", "\n")

  # get rid of sites with missing data
  gt <- gt[!apply(gt, 1, function(y) "." %in% y),]
  cat(nrow(gt), "sites remaining after sites with missing genotypes excluded.", "\n")

  if(mode == '-x'){
    # get rid of heterozygous sites
    het_sites <- apply(is_het(gt, na_is_false = FALSE), 1, function(x) TRUE %in% x)
    het_cnt <- length(which(het_sites))
    gt <- gt[!het_sites, ]
    cat(het_cnt, "heterozygous sites are discarded.", "\n")
    cat(nrow(gt), "sites remain", "\n")
    }

  gt_res <- apply(gt, 1:2, collapseGT)
  gt_res <- gt_res[apply(gt_res, 1, function(xx) length(unique(xx))) != 1, ]
  cat(nrow(gt_res), "sites remain after non variant sites filtered", "\n")
  
  cp_char_mtx <- t(gt_res)
  seq_as_phyDat <- phyDat(cp_char_mtx)
  
  cat("writing to a phylip file...", "\n")
  write.phyDat(seq_as_phyDat, outPath, format = "phylip")
}

main <- function(arguments){

  syn1 <- length(arguments) != 4
  syn2 <- arguments[1]!='-r' & arguments[1]!='-x'
  bad_syntax <- syn1 | syn2
  if(bad_syntax){
    cat("Usage: $ ./vcf2phylip.R [het_res] vcf_file out_phylip [random seed]", "\n")
    cat("\n")
    cat("[het_res] has two options", "\n")
    cat("-r randomly sample haplotypes", "\n")
    cat("-x exclude all heterozygous sites", "\n")
    return(-1)
  }
  set.seed(strtoi(arguments[4]))
  vcf2Phylip(arguments[1], arguments[2], arguments[3])
  return(0)
}

cmd_args <- commandArgs(trailingOnly = TRUE)
main(cmd_args)
#!/usr/bin/env Rscript

cmd_args <- commandArgs()

vcf2MultiAlignFile <- function(inPath, refPath, outPath, fileType){
  # Use 'inPath' to specify a vcf file, 'refPath' to specify the reference fasta file
  # This function outputs a NEXUS/phylip file at 'outPath' containing the haplotypes of each sample (heterozygous sites filtered out) in the input vcf file.
  # fileType can be 'NEXUS' or 'Phylip'
  
  library(vcfR)
  library(ape)
  library(phangorn)
  
  ref <- read.FASTA(refPath)
  vcfFile <- read.vcfR(inPath)
  
  gt <- extract.gt(vcfFile)
  tot_cnt <- nrow(gt)
  # apply all filters
  flt <- getFILTER(vcfFile)
  passAll <- which(flt=="PASS")
  vcfFile <- vcfFile[passAll]
  pass_cnt <- length(passAll)
  cat(pass_cnt, "sites out of", tot_cnt, "sites pass all filters applied.")

  # get rid of heterozygous sites  
  het_gt <- is_het(extract.gt(vcfFile), na_is_false = FALSE)
  het_sites <- apply(het_gt, 1, function(x) TRUE %in% x)
  het_cnt <- length(which(het_sites))
  vcfFile <- vcfFile[!het_sites]
  cat(het_cnt, "heterozygous sites out of", nrow(het_gt), "sites are discarded.")
  # get rid of sites with missing data
  vcfFile <- vcfFile[complete.cases(extract.gt(vcfFile))]
  cat(nrow(extract.gt(vcfFile)), "sites remaining after sites with missing genotypes excluded.")
  
  sampleSeq <- vcfR2DNAbin(vcfFile, ref.seq = ref, unphased_as_NA = FALSE, extract.haps = FALSE, consensus = TRUE)
  seq_as_phyDat <- phyDat(sampleSeq)
  
  if(fileType==".nex"){
    cat("writing to a NEXUS file...")
    write.phyDat(seq_as_phyDat, outPath, format = "nexus")
  }else if(fileType==".phy"){
    cat("writing to a phylip file...")
    write.phyDat(seq_as_phyDat, outPath, format = "phylip")
  }
}

main <- function(arguments){
  usr_args <- tail(arguments, 3)
  last_arg <- tail(usr_args, 1)
  f_ext <- substr(last_arg, nchar(last_arg)-3, nchar(last_arg))
  #print(usr_args)
  syn1 <- substr(last_arg, 1,1)=='-'
  syn2 <- f_ext!='.nex' & f_ext!='.phy'
  bad_syntax <- syn1 | syn2
  if(bad_syntax){
    cat("Usage: $ ./makeNEXUS.R vcf_file ref_fasta out.ext")
    cat("Output file extension can be either \'nex\' or \'phy\'")
    return()
  }
  vcf2MultiAlignFile(usr_args[1], usr_args[2], usr_args[3], f_ext)
}

main(cmd_args)
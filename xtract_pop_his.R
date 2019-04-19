#!/usr/bin/env Rscript

library(dplyr)

cmd_args <- commandArgs(trailingOnly = TRUE)

get_Ne_by_gen<-function(psmc_file,it=30,mu=2.5e-8,s=100,maxGBP=2e5){
  X<-scan(file=psmc_file,what="",sep="\n",quiet=TRUE)
  
  START<-grep("^RD",X)
  END<-grep("^//",X)
  
  X<-X[START[it+1]:END[it+1]]
  
  TR<-grep("^TR",X,value=TRUE)
  RS<-grep("^RS",X,value=TRUE)
  
  t_r <- as.double(unlist(strsplit(TR, '\t'))[2:3])
  #write(TR,"temp.psmc.result")
  #theta0<-as.numeric(read.table("temp.psmc.result")[1,2])
  N0<-t_r[1]/4/mu/s
  r = t_r[2]/(4*N0)
  
  write(RS,"temp.psmc.result")
  a<-read.table("temp.psmc.result")
  Generation<-round(as.numeric(2*N0*a[,3]))
  Ne<-round(as.numeric(N0*a[,4]))
  
  mask <- Generation < maxGBP
  Generation<-Generation[mask]
  Ne<-Ne[mask]
  
  file.remove("temp.psmc.result")
  
  goi <- c(maxGBP, rev(Generation[Ne!=lag(Ne)][-1]), 0)
  dNe <- c(tail(Ne, n=1), rev(Ne[Ne!=lead(Ne)])[-1], Ne[1])
  
  goi <- maxGBP - goi
  
  gen_Ne <- cbind(goi, dNe)
  #data.frame(round(Generation),round(Ne))
  return(list(c(mu, r), gen_Ne))
}

main <- function(arguments){

  if(length(arguments) != 3){
    cat("Usage: $ ./xtract_pop_his.R <psmc_file> <out_file> <max_gen>", "\n")
    return(-1)
  }
  #usr_args <- tail(arguments, 3)
  usr_args <- arguments
  
  params <- get_Ne_by_gen(usr_args[1], maxGBP=strtoi(usr_args[3]))
  write.table(params[[1]], file = usr_args[2], row.names = FALSE, col.names = FALSE)
  write.table(params[[2]], file = usr_args[2], append = TRUE, row.names = FALSE, col.names = FALSE)
  
  return(0)
}

main(cmd_args)
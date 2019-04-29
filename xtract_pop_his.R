#!/usr/bin/env Rscript

library(dplyr)

cmd_args <- commandArgs(trailingOnly = TRUE)

get_Ne_by_gen<-function(psmc_file, scheme='argw',it=30,mu=2.5e-8,s=100,maxGBP=2e5,scaling=1){
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
  r = t_r[2]/(4*N0*s)
  
  write(RS,"temp.psmc.result")
  a<-read.table("temp.psmc.result")
  Generation<-round(as.numeric(2*N0*a[,3]))
  Ne<-round(as.numeric(N0*a[,4]))
  
  mask <- Generation < maxGBP
  Generation<-Generation[mask]
  Ne<-Ne[mask]
  
  file.remove("temp.psmc.result")
  
  if(scheme=='slim'){
    goi <- c(maxGBP, rev(Generation[Ne!=lag(Ne)][-1]), 0)
    dNe <- c(tail(Ne, n=1), rev(Ne[Ne!=lead(Ne)])[-1], Ne[1])
    goi <- maxGBP - goi
  } else if (scheme=='argw'){
    goi <- Generation[Ne!=lag(Ne)]
    goi[1] <- Generation[1]
    dNe <- Ne[Ne!=lag(Ne)]
    dNe[1] <- Ne[1]
  }
  
  gen_Ne <- cbind(round(goi/scaling), round(dNe/scaling))
  #data.frame(round(Generation),round(Ne))
  return(list(c(mu*scaling, r*scaling), gen_Ne))
}

main <- function(arguments){

  if(length(arguments) != 5){
    cat("Usage: $ ./xtract_pop_his.R <psmc_file> <out_file> <`argw` or `slim`> <scaling> <max_gen>", "\n")
    return(-1)
  }
  #usr_args <- tail(arguments, 3)
  usr_args <- arguments
  
  params <- get_Ne_by_gen(usr_args[1], usr_args[3], scaling=strtoi(usr_args[4]), maxGBP=strtoi(usr_args[5]))
  write.table(params[[1]], file = usr_args[2], row.names = FALSE, col.names = FALSE)
  write.table(params[[2]], file = usr_args[2], append = TRUE, row.names = FALSE, col.names = FALSE)
  
  return(0)
}

main(cmd_args)
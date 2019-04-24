#!/usr/bin/env Rscript

library(dplyr)

psmc2discoal<-function(psmc_file,it=30,mu=2.5e-8,s=100,maxGBP=2e5){
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
  
  N_curr <- Ne[1]
  
  goi <- Generation[Ne!=lag(Ne)][-1]
  goi <- round(goi/(4*N_curr), digits = 6)
  dNe <- Ne[Ne!=lag(Ne)][-1]
  dNe <- round(dNe/N_curr, digits = 6)
  flag <- rep('-en', length(goi))
  popnID <- rep(0, length(goi))
  gen_Ne <- cbind(flag, goi, popnID, dNe)
  
  PopSizeChange <- paste(apply(gen_Ne, 1, function(x) paste(x, collapse = ' ')), collapse = ' ')
  
  discoal_PATH <- "discoal/discoal"
  sampSize <- 48
  noRep <- 5000
  nSites <- 1e5
  theta <- signif(4*N_curr*mu, 6)
  rho <- signif(4*N_curr*r, 4)
  tau <- 0
  selcoef_min <- 0.05
  selcoef_max <- 0.15
  alpha_min <- 2*N_curr*selcoef_min
  alpha_max <- 2*N_curr*selcoef_max
  AF_min <- 0.2
  AF_max <- 1
  
  discoal_cmd <- paste(discoal_PATH, sampSize, noRep, format(nSites, scientific = FALSE),
        '-t', theta, '-r', rho, '-N', N_curr,
        '-ws', tau, '-x', 0.5, '-i', 4,
        '-Pa', alpha_min, alpha_max, '-Pc', AF_min, AF_max,
        PopSizeChange, '-T\n', sep = ' ')
  
  return(discoal_cmd)
  
}

main <- function(arguments){

  if(length(arguments) != 1){
    cat("Usage: $ ./psmc2discoal.R <psmc_file>", "\n")
    return(-1)
  }
  #usr_args <- tail(arguments, 3)
  usr_args <- arguments
  
  cat(psmc2discoal(usr_args[1]))
  
  return(0)
}

cmd_args <- commandArgs(trailingOnly = TRUE)
main(cmd_args)
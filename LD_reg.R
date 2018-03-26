#!/usr/bin/env Rscript

rm(list = ls())

formatData <- function(LDstat){
  # LDstat - a dataframe read from vcftools --geno-r2 output
  # returns data points of distance-r2 pairs
  
  dist <- abs(LDstat$POS2 - LDstat$POS1)
  R2 <- LDstat$R.2
  
  return(cbind(dist, R2))
}

nl_fit <- function(dist_r2, sample_size){
  ## non-linear model fit adapted from
  ## https://fabiomarroni.wordpress.com/2011/08/09/estimate-decay-of-linkage-disequilibrium-with-distance/
  distance <- dist_r2[,1]
  LD.data <- dist_r2[,2]
  n <- sample_size
  HW.st<-c(C=0.1)
  
  HW.nonlinear<-nls(LD.data~((10+C*distance)/((2+C*distance)*(11+C*distance)))*(1+((3+C*distance)*(12+12*C*distance+(C*distance)^2))/(n*(2+C*distance)*(11+C*distance))),start=HW.st,control=nls.control(maxiter=100))
  
  tt<-summary(HW.nonlinear)
  new.rho<-tt$parameters[1]
  fit_dist <- seq(1, 100000,20)
  fpoints<-((10+new.rho*fit_dist)/((2+new.rho*fit_dist)*(11+new.rho*fit_dist)))*(1+((3+new.rho*fit_dist)*(12+12*new.rho*fit_dist+(new.rho*fit_dist)^2))/(n*(2+new.rho*fit_dist)*(11+new.rho*fit_dist)))
  fit_curve <- cbind(fit_dist, fpoints)
  
  h.decay <- max(fpoints)/2
  half.decay.distance <- fit_dist[which.min(abs(fpoints-h.decay))]
  
  return(list(HW.nonlinear, fit_curve, c(half.decay.distance, h.decay)))
}


main <- function(arguments){
  
  inval1 = length(arguments) != 3
  inval2 = arguments[1] != '-m' & arguments[1] != '-c'
  inval = inval1 | inval2
  if(inval){
    cat("Usage: $ ./LD_reg.R [option] in.geno.ld out.Rds", "\n")
    cat("[option]:", "\n")
    cat("\t -m \t save fitted model and curve", "\n")
    cat("\t -c \t save curve only", "\n")
    return(-1)
  }

  intrachromLD <- read.table(arguments[2], header = TRUE, sep = "\t")
  dist_LD <- formatData(intrachromLD)
  LD_fit <- nl_fit(dist_LD, max(intrachromLD$N_INDV))

  if (arguments[1] == '-c'){
  	LD_fit[[1]] <- NULL
  }

  saveRDS(LD_fit, file = arguments[3])
  
  return(0)
}

cmd_args <- commandArgs(trailingOnly = TRUE)
main(cmd_args)
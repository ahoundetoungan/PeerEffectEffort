#' Identifying peer effects on academic achievements through students’ effort
#' Elysée Aristide Houndetoungan and Cristelle Kouame
#' 
#' This file computes B-splines basis for the endogenous estimation

rm(list = ls())
library(PartialNetwork)
library(splines)

proot <- c("~/GPAeffort",
           "~/Dropbox/Papers - In progress/EffortGPA/Code-EffortGPA")
root  <- sapply(proot, dir.exists)
root  <- proot[root][1]
setwd(root)

# load objects
Rcpp::sourceCpp("codefiles/SourceCpp.cpp")
source("codefiles/SourceR.R")
load(file = "AHD.rda")
load("_output/Net.FE.rda")
load("_output/mu.RE.rda")


# data
data           <- data.frame(muout.fe = homoFE$estimate$mu, muin.fe = homoFE$estimate$nu, muout.re = muout, muin.re = muin)
rm(list = c("Xlogit", "homoFE", "muin", "muout", "mydata", "G"))
gc()


fbs            <- function(x, k){
  knots        <- NULL
  if(k > 0){
    knots      <- quantile(x, probs = seq(0, 1, length.out = k + 2))
    knots      <- knots[-c(1, k + 2)]
  }
  bs(x, degree = 3L, knots = knots)
}

for (k in 0:19){
  cat(k, "/19\n")
  cnames       <- colnames(data)
  smu          <- lapply(1:4, function(x){
    out        <- fbs(data[,x], k); colnames(out) <- paste0(cnames[x], ".k", k, ".b.", 1:(k + 3))
    out
  })
  
  iF           <- lapply(1:(k + 3), function(x) lapply(1:(k + 3), function(y) smu[[1]][,x]*smu[[2]][,y]))
  iR           <- lapply(1:(k + 3), function(x) lapply(1:(k + 3), function(y) smu[[3]][,x]*smu[[4]][,y]))
  iF           <- as.data.frame(iF); 
  iR           <- as.data.frame(iR); 
  colnames(iF) <- paste0("muinout.fe.k", k, ".b.", sapply(1:(k + 3), function(x) sapply(1:(k + 3), function(y) paste0(x, ".", y))))
  colnames(iR) <- paste0("muinout.re.k", k, ".b.", sapply(1:(k + 3), function(x) sapply(1:(k + 3), function(y) paste0(x, ".", y))))
  out          <- do.call(cbind, c(smu, list(iF, iR)))
  saveRDS(out, file = paste0("_output/endogeneity/data.np", k, ".RDS"))
}

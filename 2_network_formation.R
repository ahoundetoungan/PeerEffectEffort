#' Identifying peer effects on academic achievements through students’ effort
#' Elysée Aristide Houndetoungan and Cristelle Kouame
#' 
#' This file estimates the network formation model.

rm(list = ls())
library(CDatanet)
library(AER)
library(dplyr)

proot <- c("~/GPAeffort",
           "~/Dropbox/Papers - In progress/EffortGPA/Code-EffortGPA")
root  <- sapply(proot, dir.exists)
root  <- proot[root][1]
setwd(root)

# load objects
Rcpp::sourceCpp("codefiles/SourceCpp.cpp")
source("codefiles/SourceR.R")
load(file = "AHD.rda")
gc()


# fixed effects
Xlogit         <- as.data.frame(Xlogit)
form.net       <- Formula::as.Formula(paste0(c("~ -1 +", paste0(colnames(Xlogit), collapse = " + ")), collapse = ""))
Gnet           <- lapply(G, function(x) 1*(x > 0))
init           <- list(beta = rep(0, ncol(Xlogit)),
                       mu   = rep(0, nrow(mydata)),
                       nu   = rep(0, nrow(mydata) - nsch))

rm(list = ls()[!(ls() %in% c("Gnet", "form.net", "Xlogit", "init"))])
gc()

homoFE         <- homophily.FE(network =  Gnet, formula = form.net, data = Xlogit, init = init, 
                               opt.ctr = list(maxit = 1e9, eps_f = 1e-20, eps_g = 1e-20))

save(homoFE, file = "_output/Net.FE.rda")


# Random effects
rm(list = ls()[!(ls() %in% c("Gnet", "form.net", "Xlogit"))])
gc()

set.seed(123)
out            <- homophily(network =  Gnet, formula = form.net, fixed.effects = TRUE,
                            iteration = 2e4, data = Xlogit)

save(out, file = "_output/Net.RE.rda")
muout          <- colMeans(tail(out$posterior$mu, 1e4))
muin           <- colMeans(tail(out$posterior$nu, 1e4))
save(muout, muin , file = "_output/mu.RE.rda")



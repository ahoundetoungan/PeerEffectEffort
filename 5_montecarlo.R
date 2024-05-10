#' Identifying Peer Effects with Unobserved Effort and Isolated Students
#' Aristide Houndetoungan, Cristelle Kouame, and Michael Vlassopoulos
#' 
#' This file replicates Monte Carlo Simulation results

# Libraries
rm(list = ls())
library(AER)
library(PartialNetwork)
library(MASS)
library(doParallel)

# Please, add your working directory in proot
proot <- c("~/GPAeffort",
           "~/Dropbox/Papers - In progress/EffortGPA/Code-EffortGPA")
root  <- sapply(proot, dir.exists)
root  <- proot[root][1]
setwd(root)

# load objects
Rcpp::sourceCpp("codefiles/SourceCpp.cpp")
source("codefiles/SourceR.R")

# Parameters and sizes
lambd <- 0.7
beta  <- c(1, 1.5)
gamma <- c(5, -3)
s2eta <- 15
s2eps <- 8
rho   <- 0.4
M     <- 20
nvec  <- rep(50, M)
nsum  <- sum(nvec)
ncum  <- c(0, cumsum(nvec))

# Function that replicate one Monte Carlo simulation
fMC   <- function(...){
  # Simulate exogenous variables 
  x1  <- rnorm(nsum, rep(runif(M, 0, 10), nvec), 4)
  x2  <- rpois(nsum, rep(runif(M, 0, 10), nvec))
  X   <- cbind(x1, x2)
  
  # Simulate Networks
  A   <- lapply(1:M, function(x){
    mat <- matrix(0, nvec[x], nvec[x])
    for (i in 1:nvec[x]) {
      nfr <- sample(0:10, 1, prob = (1/((1:11)^0.6))/sum(1/((1:11)^0.6)))
      if(nfr > 0){
        mat[i, sample((1:nvec[x])[-i], nfr)] <- 1
      }
    }
    mat
  })
  G      <- norm.network(A)
  
  # Simulate the random errors
  Sigma  <- matrix(c(s2eps, sqrt(s2eps*s2eta)*rho, sqrt(s2eps*s2eta)*rho, s2eta), 2)
  epseta <- mvrnorm(nsum, mu = c(0, 0), Sigma = Sigma)
  eps    <- epseta[,1]
  eta    <- epseta[,2]
    
  # Simulate effort
  c0     <- -1.5*rep(sapply(1:M, function(x) quantile(x2[(ncum[x] + 1):ncum[x + 1]], 0.9)), nvec)
  effort <- c(peer.avg(lapply(1:M, function(x) solve(diag(nvec[x]) - lambd*G[[x]])), 
                     c0 + X%*%beta + peer.avg(G, X%*%gamma) + eps))
  
  # Compute GPA
  gpaA    <- effort + eta # DGP A: alpha = 0
  alpha0  <- 10*rep(sapply(1:M, function(x) quantile(x1[(ncum[x] + 1):ncum[x + 1]], 0.9)), nvec)
  gpaC    <- alpha0 + effort + eta # DGP C: alpha vary across schools
  gpaB    <- mean(alpha0) + effort + eta # DGP B: alpha does not vary across schools
  
  # Estimation following the file 1_exogenous_networks.R
  YX            <- cbind(gpaA, gpaB, gpaC, X); colnames(YX) <- c("gpaA", "gpaB", "gpaC", "x1", "x2")
  GXY           <- peer.avg(G, YX); colnames(GXY) <- paste0("G_", colnames(YX))
  GGX           <- peer.avg(G, GXY[,-(1:3)]); colnames(GGX) <- paste0("G", colnames(GXY[,-(1:3)]))
  hasfriends    <- unlist(lapply(G, rowSums))
  Ghasfriends   <- peer.avg(G, hasfriends)
  
  # matrix J
  J1            <- lapply(1:M, function(x) diag(nvec[x]) - matrix(1, nvec[x], nvec[x])/nvec[x])
  J3            <- lapply(1:M, function(x) {
    onehat      <- rowSums(G[[x]])
    onecheck    <- 1 - onehat
    sonehat     <- sum(onehat); sonehat     <- ifelse(sonehat == 0, 1, sonehat)
    sonecheck   <- sum(onecheck); sonecheck <-  ifelse(sonecheck == 0, 1, sonecheck)
    diag(nvec[x]) - onehat %*% t(onehat)/sonehat - onecheck %*% t(onecheck)/sonecheck
  })
  F1            <- fdataFs(J1)
  F3            <- fdataFs(J3)
  
  # Multiply by J
  F1XY    <- peer.avg(F1, cbind(YX, hasfriends)); colnames(F1XY) <- paste0("F1_", c(colnames(YX), "hasfriends"))
  F3XY    <- peer.avg(F3, YX); colnames(F3XY) <- paste0("F3_", colnames(YX))
  F1GXY   <- peer.avg(F1, cbind(GXY, Ghasfriends)); colnames(F1GXY) <- paste0("F1", c(colnames(GXY), "G_hasfriends"))
  F3GXY   <- peer.avg(F3, GXY); colnames(F3GXY) <- paste0("F3", colnames(GXY))
  F1GGX   <- peer.avg(F1, GGX); colnames(F1GGX) <- paste0("F1", colnames(GGX))
  F3GGX   <- peer.avg(F3, GGX); colnames(F3GGX) <- paste0("F3", colnames(GGX))
  mydata0 <- data.frame(YX, GXY, GGX)
  mydata1 <- data.frame(F1XY, F1GXY, F1GGX)
  mydata3 <- data.frame(F3XY, F3GXY, F3GGX)
  
  # GMM
  ## formula and instruments
  va.exo  <- c("x1", "x2")
  formA0  <- as.formula(paste("gpaA ~", paste(c("G_gpaA", va.exo, paste0("G_", va.exo)), collapse = "+")))
  formB0  <- as.formula(paste("gpaB ~", paste(c("G_gpaB", va.exo, paste0("G_", va.exo)), collapse = "+")))
  formC0  <- as.formula(paste("gpaC ~", paste(c("G_gpaC", va.exo, paste0("G_", va.exo)), collapse = "+")))
  instr0  <- as.formula(paste("~", paste(c(va.exo, paste0("G_", va.exo), paste0("GG_", va.exo)), collapse = "+")))
  formA1  <- as.formula(paste("F1_gpaA ~", paste(c(-1, "F1G_gpaA", paste0("F1_", va.exo), paste0("F1G_", va.exo)), collapse = "+")))
  formB1  <- as.formula(paste("F1_gpaB ~", paste(c(-1, "F1G_gpaB", paste0("F1_", va.exo), paste0("F1G_", va.exo)), collapse = "+")))
  formC1  <- as.formula(paste("F1_gpaC ~", paste(c(-1, "F1G_gpaC", paste0("F1_", va.exo), paste0("F1G_", va.exo)), collapse = "+")))
  instr1  <- as.formula(paste("~", paste(c(-1, paste0("F1_", va.exo), paste0("F1G_", va.exo), paste0("F1GG_", va.exo)), collapse = "+")))
  formA2  <- as.formula(paste("F1_gpaA ~", paste(c(-1, "F1G_gpaA", "F1_hasfriends", paste0("F1_", va.exo), paste0("F1G_", va.exo)), collapse = "+")))
  formB2  <- as.formula(paste("F1_gpaB ~", paste(c(-1, "F1G_gpaB", "F1_hasfriends", paste0("F1_", va.exo), paste0("F1G_", va.exo)), collapse = "+")))
  formC2  <- as.formula(paste("F1_gpaC ~", paste(c(-1, "F1G_gpaC", "F1_hasfriends", paste0("F1_", va.exo), paste0("F1G_", va.exo)), collapse = "+")))
  instr2  <- as.formula(paste("~", paste(c(-1, paste0("F1_", va.exo), paste0("F1G_", va.exo), paste0("F1GG_", va.exo)), collapse = "+"), "+ F1_hasfriends + F1G_hasfriends"))
  formA3  <- as.formula(paste("F3_gpaA ~", paste(c(-1, "F3G_gpaA", paste0("F3_", va.exo), paste0("F3G_", va.exo)), collapse = "+")))
  formB3  <- as.formula(paste("F3_gpaB ~", paste(c(-1, "F3G_gpaB", paste0("F3_", va.exo), paste0("F3G_", va.exo)), collapse = "+")))
  formC3  <- as.formula(paste("F3_gpaC ~", paste(c(-1, "F3G_gpaC", paste0("F3_", va.exo), paste0("F3G_", va.exo)), collapse = "+")))
  instr3  <- as.formula(paste("~", paste(c(-1, paste0("F3_", va.exo), paste0("F3G_", va.exo), paste0("F3GG_", va.exo)), collapse = "+")))
  

  out     <- c(ivreg(formula = formA0, instruments = instr0, data = mydata0)$coefficients,
               ivreg(formula = formB0, instruments = instr0, data = mydata0)$coefficients,
               ivreg(formula = formC0, instruments = instr0, data = mydata0)$coefficients,
               ivreg(formula = formA1, instruments = instr1, data = mydata1)$coefficients,
               ivreg(formula = formB1, instruments = instr1, data = mydata1)$coefficients,
               ivreg(formula = formC1, instruments = instr1, data = mydata1)$coefficients)
  
  outA2   <- ivreg(formula = formA2, instruments = instr2, data = mydata1)
  outB2   <- ivreg(formula = formB2, instruments = instr2, data = mydata1)
  outC2   <- ivreg(formula = formC2, instruments = instr2, data = mydata1)
  
  mlA2    <- foptim(outA2$residuals, outA2$coefficients["F1G_gpaA"], G, fixed.effects = TRUE, F1, start = c(1.37, 0.4))
  mlB2    <- foptim(outB2$residuals, outB2$coefficients["F1G_gpaB"], G, fixed.effects = TRUE, F1, start = c(1.37, 0.4))
  mlC2    <- foptim(outA2$residuals, outC2$coefficients["F1G_gpaC"], G, fixed.effects = TRUE, F1, start = c(1.37, 0.4))

  outA3   <- ivreg(formula = formA3, instruments = instr3, data = mydata3)
  outB3   <- ivreg(formula = formB3, instruments = instr3, data = mydata3)
  outC3   <- ivreg(formula = formC3, instruments = instr3, data = mydata3)
  
  mlA3    <- foptim(outA3$residuals, outA3$coefficients["F3G_gpaA"], G, fixed.effects = TRUE, F3, start = c(1.37, 0.4))
  mlB3    <- foptim(outB3$residuals, outB3$coefficients["F3G_gpaB"], G, fixed.effects = TRUE, F3, start = c(1.37, 0.4))
  mlC3    <- foptim(outA3$residuals, outC3$coefficients["F3G_gpaC"], G, fixed.effects = TRUE, F3, start = c(1.37, 0.4))
  c(out, outA2$coefficients, unlist(mlA2[c("sigma2epsilon", "sigma2eta", "rho")]), 
    outB2$coefficients, unlist(mlB2[c("sigma2epsilon", "sigma2eta", "rho")]), 
    outC2$coefficients, unlist(mlC2[c("sigma2epsilon", "sigma2eta", "rho")]), 
    outA3$coefficients, unlist(mlA3[c("sigma2epsilon", "sigma2eta", "rho")]), 
    outB3$coefficients, unlist(mlB3[c("sigma2epsilon", "sigma2eta", "rho")]), 
    outC3$coefficients, unlist(mlC3[c("sigma2epsilon", "sigma2eta", "rho")]))
}

RNGkind("L'Ecuyer-CMRG")
set.seed(1234)

outMC <- t(apply(do.call(cbind, mclapply(1:1000, function(x){cat("Iteration: ", x, "\n", sep = ""); fMC()}, mc.cores = 10)), 1, 
               function(x) c(mean = mean(x), sderr = sd(x), "Pctl 25%" = quantile(x, 0.25), "Pctl 75%" = quantile(x, 0.75))))
outMC

write.csv(outMC, file = "_output/simu.alpha.csv")

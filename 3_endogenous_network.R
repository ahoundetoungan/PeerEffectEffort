#' Identifying Peer Effects with Unobserved Effort and Isolated Students
#' Aristide Houndetoungan, Cristelle Kouame, and Michael Vlassopoulos
#' 
#' This file replicates the peer effect model estimation controlling for network endogeneity:
#' the nonparametric approach.

rm(list = ls())
library(PartialNetwork)
library(AER)
library(dplyr)
library(splines)

proot <- c("~/GPAeffort",
           "~/Dropbox/Papers - In progress/EffortGPA/Code-EffortGPA")
root  <- sapply(proot, dir.exists)
root  <- proot[root][1]
setwd(root)

################### Computation of B-spline
# load objects
Rcpp::sourceCpp("codefiles/SourceCpp.cpp")
source("codefiles/SourceR.R")
load(file = "../../../Data/AHdata/PEEffort/AHDgpa.rda")
load("_output/Net.FE.rda")
load("_output/mu.RE.rda")


# data
data         <- data.frame(muout.fe = homoFE$estimate$mu, muin.fe = homoFE$estimate$nu, muout.re = muout, muin.re = muin)
rm(list = c("Xlogit", "homoFE", "muin", "muout", "mydata", "G"))
gc()


fbs          <- function(x, k){
  knots      <- NULL
  if(k > 0){
    knots    <- quantile(x, probs = seq(0, 1, length.out = k + 2))
    knots    <- knots[-c(1, k + 2)]
  }
  bs(x, degree = 3L, knots = knots)
}

k            <- 10
cnames       <- colnames(data)
smu          <- lapply(1:4, function(x){
  out        <- fbs(data[,x], k); colnames(out) <- paste0(cnames[x], ".k", k, ".b.", 1:(k + 3))
  out
})

# Interaction between in and out
# iF           <- lapply(1:(k + 3), function(x) lapply(1:(k + 3), function(y) smu[[1]][,x]*smu[[2]][,y]))
# iR           <- lapply(1:(k + 3), function(x) lapply(1:(k + 3), function(y) smu[[3]][,x]*smu[[4]][,y]))
# iF           <- as.data.frame(iF); 
# iR           <- as.data.frame(iR); 
# colnames(iF) <- paste0("muinout.fe.k", k, ".b.", sapply(1:(k + 3), function(x) sapply(1:(k + 3), function(y) paste0(x, ".", y))))
# colnames(iR) <- paste0("muinout.re.k", k, ".b.", sapply(1:(k + 3), function(x) sapply(1:(k + 3), function(y) paste0(x, ".", y))))
# out          <- do.call(cbind, c(smu, list(iF, iR)))
out          <- do.call(cbind, smu)
saveRDS(out, file = paste0("_output/data.np", k, ".RDS"))


######################## Model estimation
# load objects
rm(list = ls())
Rcpp::sourceCpp("codefiles/SourceCpp.cpp")
source("codefiles/SourceR.R")
load(file = "../../../Data/AHdata/PEEffort/AHDgpa.rda")
rm("Xlogit")
gc()

# Network
G             <- norm.network(G)

# matrices J and F
J1            <- lapply(1:nsch, function(x) diag(sch.size[x]) - matrix(1, sch.size[x], sch.size[x])/sch.size[x])
J3            <- lapply(1:nsch, function(x) {
  onehat      <- rowSums(G[[x]])
  onecheck    <- 1 - onehat
  sonehat     <- sum(onehat); sonehat     <- ifelse(sonehat == 0, 1, sonehat)
  sonecheck   <- sum(onecheck); sonecheck <-  ifelse(sonecheck == 0, 1, sonecheck)
  diag(sch.size[x]) - onehat %*% t(onehat)/sonehat - onecheck %*% t(onecheck)/sonecheck
})
F1            <- fdataFs(J1)
F3            <- fdataFs(J3)

rm(list = c("J1", "J3"))
gc()

festim        <- function(k = 10, eff = "fe", dvar = "gpa"){
  data        <- readRDS(file = paste0("_output/data.np", k, ".RDS"))
  data        <- cbind(mydata, data)
  
  newvar      <- c(paste0("muout.", eff, ".k", k, ".b.", 1:(k + 3)), paste0("muin.", eff, ".k", k, ".b.", 1:(k + 3)))
  va.exo      <- va.names[-length(va.names)]
  va.names.n  <- c(va.names, newvar)
  va.exo.n    <- c(va.exo, newvar)
  iInst       <- 1:length(va.exo)
  
  ### With fixed effects and without dummy for isolated students
  form.FE1    <- as.formula(paste0("F1_", dvar, " ~", paste(c(-1, paste0("F1G_", dvar), paste0("F1_", va.exo), paste0("F1G_", va.exo), paste0("F1_", newvar)), collapse = "+")))
  instr.FE1   <- as.formula(paste0("~", paste(c(-1, paste0("F1_", va.exo.n), paste0("F1G_", va.exo), paste0("F1GG_", va.exo)), collapse = "+")))
  ### With fixed effects and with dummy for isolated students
  form.FE2    <- as.formula(paste0("F1_", dvar, " ~", paste(c(-1, paste0("F1G_", dvar), "F1_hasfriends", paste0("F1_", va.exo), paste0("F1G_", va.exo), paste0("F1_", newvar)), collapse = "+")))
  instr.FE2   <- as.formula(paste0("~", paste(c(-1, paste0("F1_", va.exo.n), paste0("F1G_", va.exo), paste0("F1GG_", va.exo)), collapse = "+"), "+ F1_hasfriends + F1G_hasfriends"))
  ### Without fixed effects and with dummy for isolated students per school
  form.FE3    <- as.formula(paste0("F3_", dvar, " ~", paste(c(-1, paste0("F3G_", dvar), paste0("F3_", va.exo), paste0("F3G_", va.exo), paste0("F3_", newvar)), collapse = "+")))
  instr.FE3   <- as.formula(paste0("~", paste(c(-1, paste0("F3_", va.exo.n), paste0("F3G_", va.exo), paste0("F3GG_", va.exo)), collapse = "+")))
  
  XY          <- data[, c(va.exo.n, va.names[length(va.names)])]
  GXY         <- peer.avg(G, XY); colnames(GXY) <- paste0("G_", colnames(XY))
  GGX         <- peer.avg(G, GXY[,-ncol(GXY)]); colnames(GGX) <- paste0("G", colnames(GXY[,-ncol(GXY)]))
  hasfriends  <- unlist(lapply(G, rowSums))
  Ghasfriends <- peer.avg(G, hasfriends)
  
  F1XY        <- peer.avg(F1, cbind(XY, hasfriends)); colnames(F1XY) <- paste0("F1_", c(colnames(XY), "hasfriends"))
  F3XY        <- peer.avg(F3, XY); colnames(F3XY) <- paste0("F3_", colnames(XY))
  F1GXY       <- peer.avg(F1, cbind(GXY, Ghasfriends)); colnames(F1GXY) <- paste0("F1", c(colnames(GXY), "G_hasfriends"))
  F3GXY       <- peer.avg(F3, GXY); colnames(F3GXY) <- paste0("F3", colnames(GXY))
  F1GGX       <- peer.avg(F1, GGX); colnames(F1GGX) <- paste0("F1", colnames(GGX))
  F3GGX       <- peer.avg(F3, GGX); colnames(F3GGX) <- paste0("F3", colnames(GGX))
  data        <- cbind(mydata, GXY, GGX, "hasfriends" = hasfriends, "G_hasfriends" = Ghasfriends)
  dataFE1     <- data.frame(F1XY, F1GXY, F1GGX)
  dataFE3     <- data.frame(F3XY, F3GXY, F3GGX)
  
  rm(list = c("GXY", "GGX", "hasfriends", "Ghasfriends", "F1XY", "F3XY", "F1GXY", "F3GXY", "F1GGX", "F3GGX"))
  gc()
  
  
  ## IV estimation without fixed effects
  ### With fixed effects and without dummy for isolated students
  iv.FE1       <- ivreg(formula = form.FE1, instruments = instr.FE1, data = dataFE1)
  (siv.FE1     <- summary(iv.FE1, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F1), function(x) rep(x, nrow(F1[[x]]))))))
  saveRDS(siv.FE1, file = paste0("_output/ennp.siv2", ".", eff, ".RDS"))
  siv.FE1      <- readRDS(file = paste0("_output/ennp.siv2", ".", eff, ".RDS"))
  mean(siv.FE1$residuals^2)
  write.csv(siv.FE1$coefficients, file = paste0("_output/ennp.siv2", ".", eff, ".csv"))
  
  ### With fixed effects and with dummy for isolated students
  iv.FE2       <- ivreg(formula = form.FE2, instruments = instr.FE2, data = dataFE1)
  (siv.FE2     <- summary(iv.FE2, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F1), function(x) rep(x, nrow(F1[[x]]))))))
  mean(siv.FE2$residuals^2)
  saveRDS(iv.FE2, file = paste0("_output/ennp.iv3", ".", eff, ".RDS"))
  iv.FE2       <- readRDS(file = paste0("_output/ennp.iv3", ".", eff, ".RDS"))
  ml.iv.FE2    <- foptim(iv.FE2$residuals, iv.FE2$coefficients[paste0("F1G_", dvar)], 
                         G, fixed.effects = TRUE, F1, start = c(2.491167, 0.6053381))
  va.iv.FE2    <- fvariance.iv(iv.FE2, ml.iv.FE2, dvar)  
  va.iv.FE2$output
  va.iv.FE2$sargan.stat
  va.iv.FE2$sargan.pvalue
  va.iv.FE2$sigma2eta
  va.iv.FE2$sigma2epsilon
  va.iv.FE2$rho
  saveRDS(va.iv.FE2, file = paste0("_output/ennp.va.iv3", ".", eff, ".RDS"))
  va.iv.FE2    <- readRDS(file = paste0("_output/ennp.va.iv3", ".", eff, ".RDS"))
  write.csv(va.iv.FE2$output, file = paste0("_output/ennp.va.iv3", ".", eff, ".csv"))
  
  ### With fixed effects and with dummy for isolated students per school
  iv.FE3       <- ivreg(formula = form.FE3, instruments = instr.FE3, data = dataFE3)
  (siv.FE3     <- summary(iv.FE3, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F3), function(x) rep(x, nrow(F3[[x]]))))))
  saveRDS(iv.FE3, file = paste0("_output/ennp.iv4", ".", eff, ".RDS"))
  iv.FE3       <- readRDS(file = paste0("_output/ennp.iv4", ".", eff, ".RDS"))
  ml.iv.FE3    <- foptim(iv.FE3$residuals, iv.FE3$coefficients[paste0("F3G_", dvar)], 
                         G, fixed.effects = TRUE, F3, start = c(2.491167, 0.6053381))
  va.iv.FE3    <- fvariance.iv(iv.FE3, ml.iv.FE3, dvar)  
  va.iv.FE3$output
  va.iv.FE3$sargan.stat
  va.iv.FE3$sargan.pvalue
  va.iv.FE3$sigma2eta
  va.iv.FE3$sigma2epsilon
  va.iv.FE3$rho
  saveRDS(va.iv.FE3, file = paste0("_output/ennp.va.iv4", ".", eff, ".RDS"))
  va.iv.FE3    <- readRDS(file = paste0("_output/ennp.va.iv4", ".", eff, ".RDS"))
  write.csv(va.iv.FE3$output, file = paste0("_output/ennp.va.iv4", ".", eff, ".csv"))
}

festim(k = 10, eff = "fe", dvar = "gpa")
festim(k = 10, eff = "re", dvar = "gpa")
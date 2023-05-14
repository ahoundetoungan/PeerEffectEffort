#' Identifying peer effects on academic achievements through students’ effort
#' Elysée Aristide Houndetoungan and Cristelle Kouame
#' 
#' This file replicates the peer effect model estimation assuming that the network is exogenous

rm(list = ls())
library(PartialNetwork)
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
rm("Xlogit")
gc()

mydata        <- mydata %>% mutate(cst = 1)
va.names.cst  <- c("cst", va.names)
va.exo        <- va.names[-length(va.names)]
va.exo.cst    <- c("cst", va.exo)

# Econometric part
# data
G             <- norm.network(G)
GXY           <- peer.avg(G, mydata[,va.names]); colnames(GXY) <- paste0("G_", va.names)
GGX           <- peer.avg(G, GXY[,-ncol(GXY)]); colnames(GGX) <- paste0("G", colnames(GXY[,-ncol(GXY)]))
hasfriends    <- unlist(lapply(G, rowSums))
Ghasfriends   <- peer.avg(G, hasfriends)

# matrix J
J1            <- lapply(1:nsch, function(x) diag(sch.size[x]) - matrix(1, sch.size[x], sch.size[x])/sch.size[x])
J2            <- lapply(1:nsch, function(x) {
  onehat      <- rowSums(G[[x]])
  onecheck    <- 1 - onehat
  sonehat     <- sum(onehat); sonehat     <- ifelse(sonehat == 0, 1, sonehat)
  sonecheck   <- sum(onecheck); sonecheck <-  ifelse(sonecheck == 0, 1, sonecheck)
  diag(sch.size[x]) - onehat %*% t(onehat)/sonehat - onecheck %*% t(onecheck)/sonecheck
})
F1            <- fdataFs(J1)
F2            <- fdataFs(J2)

# Multiply by J
F1XY          <- peer.avg(F1, mydata[,va.names]); colnames(F1XY) <- paste0("F1_", va.names)
F2XY          <- peer.avg(F2, mydata[,va.names]); colnames(F2XY) <- paste0("F2_", va.names)
F1GXY         <- peer.avg(F1, GXY); colnames(F1GXY) <- paste0("F1", colnames(GXY))
F2GXY         <- peer.avg(F2, GXY); colnames(F2GXY) <- paste0("F2", colnames(GXY))
F1GGX         <- peer.avg(F1, GGX); colnames(F1GGX) <- paste0("F1", colnames(GGX))
F2GGX         <- peer.avg(F2, GGX); colnames(F2GGX) <- paste0("F2", colnames(GGX))
mydata        <- cbind(mydata, GXY, GGX, "hasfriends" = hasfriends, "G_hasfriends" = Ghasfriends)
mydataFE1     <- data.frame(F1XY, F1GXY, F1GGX)
mydataFE2     <- data.frame(F2XY, F2GXY, F2GGX)

rm(list = c("GXY", "GGX", "hasfriends", "Ghasfriends", "F1XY", "F2XY", "F1GXY", "F2GXY", "F1GGX", "F2GGX", "J1", "J2"))
gc()

# GMM
## formula and instruments
form.noFE1    <- as.formula(paste("gpa ~", paste(c(-1, va.exo.cst, paste0("G_", va.exo)), collapse = "+"), "+ G_gpa"))
instr.noFE1   <- as.formula(paste("~", paste(c(-1, va.exo.cst, paste0("G_", va.exo), paste0("GG_", va.exo)), collapse = "+")))
form.noFE2    <- as.formula(paste("gpa ~", paste(c(-1, va.exo.cst, paste0("G_", va.exo)), collapse = "+"), "+ hasfriends + G_gpa"))
instr.noFE2   <- as.formula(paste("~", paste(c(-1, va.exo.cst, paste0("G_", va.exo), paste0("GG_", va.exo)), collapse = "+"), "+ hasfriends + G_hasfriends"))
form.FE1      <- as.formula(paste("F1_gpa ~", paste(c(-1, paste0("F1_", va.exo), paste0("F1G_", va.exo)), collapse = "+"), "+ F1G_gpa"))
instr.FE1     <- as.formula(paste("~", paste(c(-1, paste0("F1_", va.exo), paste0("F1G_", va.exo), paste0("F1GG_", va.exo)), collapse = "+")))
form.FE2      <- as.formula(paste("F2_gpa ~", paste(c(-1, paste0("F2_", va.exo), paste0("F2G_", va.exo)), collapse = "+"), "+ F2G_gpa"))
instr.FE2     <- as.formula(paste("~", paste(c(-1, paste0("F2_", va.exo), paste0("F2G_", va.exo), paste0("F2GG_", va.exo)), collapse = "+")))

## OLS estimation without fixed effects
ols.noFE1    <- lm(formula = form.noFE1, data = mydata)
(sols.noFE1  <- summary(ols.noFE1, vcov = vcovCL, cluster = unlist(lapply(1:length(G), function(x) rep(x, nrow(G[[x]]))))))
saveRDS(sols.noFE1, file = "_output/sols.noFE1.RDS")
sols.noFE1   <- readRDS(file = "_output/sols.noFE1.RDS")

ols.noFE2    <- lm(formula = form.noFE2, data = mydata)
(sols.noFE2  <- summary(ols.noFE2, vcov = vcovCL, cluster = unlist(lapply(1:length(G), function(x) rep(x, nrow(G[[x]]))))))
saveRDS(ols.noFE2, file = "_output/ols.noFE2.RDS")
ols.noFE2    <- readRDS(file = "_output/ols.noFE2.RDS")
ml.ols.noFE2 <- foptim(ols.noFE2$residuals, ols.noFE2$coefficients["G_gpa"], 
                       G, fixed.effects = FALSE, F2, start = c(1, 0))
va.ols.noFE2 <- fvariance.ols(ols.noFE2, ml.ols.noFE2)  
va.ols.noFE2$output
saveRDS(va.ols.noFE2, file = "_output/va.ols.noFE2.RDS")
va.ols.noFE2 <- readRDS(file = "_output/va.ols.noFE2.RDS")

## OLS estimation with fixed effects
ols.FE1      <- lm(formula = form.FE1, data = mydataFE1)
(sols.FE1    <- summary(ols.FE1, vcov = vcovCL, cluster = unlist(lapply(1:length(F1), function(x) rep(x, nrow(F1[[x]]))))))
saveRDS(sols.FE1, file = "_output/sols.FE1.RDS")
sols.FE1     <- readRDS(file = "_output/sols.FE1.RDS")

ols.FE2      <- lm(formula = form.FE2, data = mydataFE2)
(sols.FE2    <- summary(ols.FE2, vcov = vcovCL, cluster = unlist(lapply(1:length(F2), function(x) rep(x, nrow(F2[[x]]))))))
saveRDS(ols.FE2, file = "_output/ols.FE2.RDS")
ols.FE2      <- readRDS(file = "_output/ols.FE2.RDS")
ml.ols.FE2   <- foptim(ols.FE2$residuals, ols.FE2$coefficients["F2G_gpa"], 
                       G, fixed.effects = TRUE, F2, start = c(1, 0))
va.ols.FE2   <- fvariance.ols(ols.FE2, ml.ols.FE2)  
va.ols.FE2$output
saveRDS(va.ols.FE2, file = "_output/va.ols.FE2.RDS")
va.ols.FE2   <- readRDS(file = "_output/va.ols.FE2.RDS")

## IV estimation without fixed effects
iv.noFE1     <- ivreg(formula = form.noFE1, instruments = instr.noFE1, data = mydata)
(siv.noFE1   <- summary(iv.noFE1, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(G), function(x) rep(x, nrow(G[[x]]))))))
saveRDS(siv.noFE1, file = "_output/siv.noFE1.RDS")
siv.noFE1    <- readRDS(file = "_output/siv.noFE1.RDS")
write.csv(siv.noFE1$coefficients, file = "_output/siv.noFE1.csv")

iv.noFE2     <- ivreg(formula = form.noFE2, instruments = instr.noFE2, data = mydata)
(siv.noFE2   <- summary(iv.noFE2, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(G), function(x) rep(x, nrow(G[[x]]))))))
saveRDS(iv.noFE2, file = "_output/iv.noFE2.RDS")
iv.noFE2     <- readRDS(file = "_output/iv.noFE2.RDS")
ml.iv.noFE2  <- foptim(iv.noFE2$residuals, iv.noFE2$coefficients["G_gpa"], 
                       G, fixed.effects = FALSE, F2, start = c(1, 0))
va.iv.noFE2  <- fvariance.iv(iv.noFE2, ml.iv.noFE2)  
va.iv.noFE2$output
va.iv.noFE2$sargan.stat
va.iv.noFE2$sargan.pvalue
va.iv.noFE2$sigma2eta
va.iv.noFE2$sigma2epsilon
va.iv.noFE2$rho
saveRDS(va.iv.noFE2, file = "_output/va.iv.noFE2.RDS")
va.iv.noFE2  <- readRDS(file = "_output/va.iv.noFE2.RDS")
write.csv(va.iv.noFE2$output, file = "_output/va.iv.noFE2.csv")

## IV estimation with fixed effects
iv.FE1       <- ivreg(formula = form.FE1, instruments = instr.FE1, data = mydataFE1)
(siv.FE1     <- summary(iv.FE1, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F1), function(x) rep(x, nrow(F1[[x]]))))))
saveRDS(siv.FE1, file = "_output/siv.FE1.RDS")
siv.FE1      <- readRDS(file = "_output/siv.FE1.RDS")
write.csv(siv.FE1$coefficients, file = "_output/siv.FE1.csv")

iv.FE2       <- ivreg(formula = form.FE2, instruments = instr.FE2, data = mydataFE2)
(siv.FE2     <- summary(iv.FE2, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F2), function(x) rep(x, nrow(F2[[x]]))))))
saveRDS(iv.FE2, file = "_output/iv.FE2.RDS")
iv.FE2       <- readRDS(file = "_output/iv.FE2.RDS")
ml.iv.FE2    <- foptim(iv.FE2$residuals, iv.FE2$coefficients["F2G_gpa"], 
                       G, fixed.effects = TRUE, F2, start = c(1, 0))
va.iv.FE2    <- fvariance.iv(iv.FE2, ml.iv.FE2)  
va.iv.FE2$output
va.iv.FE2$sargan.stat
va.iv.FE2$sargan.pvalue
va.iv.FE2$sigma2eta
va.iv.FE2$sigma2epsilon
va.iv.FE2$rho
saveRDS(va.iv.FE2, file = "_output/va.iv.FE2.RDS")
va.iv.FE2    <- readRDS(file = "_output/va.iv.FE2.RDS")
write.csv(va.iv.FE2$output, file = "_output/va.iv.FE2.csv")

## IV estimation with fixed effects removing completely isolated students
rm(list = ls())
Rcpp::sourceCpp("codefiles/SourceCpp.cpp")
source("codefiles/SourceR.R")
load(file = "AHD.rda")
rm("Xlogit")
gc()

mydata        <- mydata %>% mutate(cst = 1) 
va.names.cst  <- c("cst", va.names)
va.exo        <- va.names[-length(va.names)]
va.exo.cst    <- c("cst", va.exo)

### We remove completely isolated students
nhafr         <- lapply(G, rowSums) # number of friends
nisfr         <- lapply(G, colSums) # number of times the student is a friend
keep          <- lapply(1:nsch, function(x) (nhafr[[x]] > 0) | (nisfr[[x]] > 0)) #
mydata        <- mydata %>% filter(unlist(keep)) %>% group_by(sschlcde) %>% mutate(nstudent = n()) %>% ungroup() %>% filter(nstudent > 2)
G             <- norm.network(lapply(1:nsch, function(x) G[[x]][keep[[x]], keep[[x]]]))
Gnrow         <- sapply(G, nrow); G <- G[Gnrow > 2] # We remove school with completely isolated students
nsch          <- length(G)
sch.size      <- Gnrow[Gnrow > 2]

### Data
GXY           <- peer.avg(G, mydata[,va.names]); colnames(GXY) <- paste0("G_", va.names)
GGX           <- peer.avg(G, GXY[,-ncol(GXY)]); colnames(GGX) <- paste0("G", colnames(GXY[,-ncol(GXY)]))
hasfriends    <- unlist(lapply(G, rowSums))
Ghasfriends   <- peer.avg(G, hasfriends)
J1            <- lapply(1:nsch, function(x) diag(sch.size[x]) - matrix(1, sch.size[x], sch.size[x])/sch.size[x])
J2            <- lapply(1:nsch, function(x) {
  onehat      <- rowSums(G[[x]])
  onecheck    <- 1 - onehat
  sonehat     <- sum(onehat); sonehat     <- ifelse(sonehat == 0, 1, sonehat)
  sonecheck   <- sum(onecheck); sonecheck <-  ifelse(sonecheck == 0, 1, sonecheck)
  diag(sch.size[x]) - onehat %*% t(onehat)/sonehat - onecheck %*% t(onecheck)/sonecheck
})
F1            <- fdataFs(J1)
F2            <- fdataFs(J2)

F1XY          <- peer.avg(F1, mydata[,va.names]); colnames(F1XY) <- paste0("F1_", va.names)
F2XY          <- peer.avg(F2, mydata[,va.names]); colnames(F2XY) <- paste0("F2_", va.names)
F1GXY         <- peer.avg(F1, GXY); colnames(F1GXY) <- paste0("F1", colnames(GXY))
F2GXY         <- peer.avg(F2, GXY); colnames(F2GXY) <- paste0("F2", colnames(GXY))
F1GGX         <- peer.avg(F1, GGX); colnames(F1GGX) <- paste0("F1", colnames(GGX))
F2GGX         <- peer.avg(F2, GGX); colnames(F2GGX) <- paste0("F2", colnames(GGX))
mydata        <- cbind(mydata, GXY, GGX, "hasfriends" = hasfriends, "G_hasfriends" = Ghasfriends)
mydataFE1     <- data.frame(F1XY, F1GXY, F1GGX)
mydataFE2     <- data.frame(F2XY, F2GXY, F2GGX)

rm(list = c("GXY", "GGX", "hasfriends", "Ghasfriends", "F1XY", "F2XY", "F1GXY", "F2GXY", "F1GGX", "F2GGX", "J1", "J2"))
gc()

### Formula and instruments
form.FE1      <- as.formula(paste("F1_gpa ~", paste(c(-1, paste0("F1_", va.exo), paste0("F1G_", va.exo)), collapse = "+"), "+ F1G_gpa"))
instr.FE1     <- as.formula(paste("~", paste(c(-1, paste0("F1_", va.exo), paste0("F1G_", va.exo), paste0("F1GG_", va.exo)), collapse = "+")))
form.FE2      <- as.formula(paste("F2_gpa ~", paste(c(-1, paste0("F2_", va.exo), paste0("F2G_", va.exo)), collapse = "+"), "+ F2G_gpa"))
instr.FE2     <- as.formula(paste("~", paste(c(-1, paste0("F2_", va.exo), paste0("F2G_", va.exo), paste0("F2GG_", va.exo)), collapse = "+")))

### Estimation
iviso.FE1     <- ivreg(formula = form.FE1, instruments = instr.FE1, data = mydataFE1)
(siviso.FE1   <- summary(iviso.FE1, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F1), function(x) rep(x, nrow(F1[[x]]))))))
saveRDS(siviso.FE1, file = "_output/siviso.FE1.RDS")
siv.FE1       <- readRDS(file = "_output/siviso.FE1.RDS")
write.csv(siv.FE1$coefficients, file = "_output/siviso.FE1.csv")

iviso.FE2     <- ivreg(formula = form.FE2, instruments = instr.FE2, data = mydataFE2)
(siviso.FE2   <- summary(iviso.FE2, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F2), function(x) rep(x, nrow(F2[[x]]))))))
saveRDS(iviso.FE2, file = "_output/iviso.FE2.RDS")
iviso.FE2     <- readRDS(file = "_output/iviso.FE2.RDS")
ml.iviso.FE2  <- foptim(iviso.FE2$residuals, iviso.FE2$coefficients["F2G_gpa"], 
                       G, fixed.effects = TRUE, F2, start = c(1, 0))
va.iviso.FE2  <- fvariance.iv(iviso.FE2, ml.iviso.FE2)  
va.iviso.FE2$output
va.iviso.FE2$sargan.stat
va.iviso.FE2$sargan.pvalue
va.iviso.FE2$sigma2eta
va.iviso.FE2$sigma2epsilon
va.iviso.FE2$rho
saveRDS(va.iviso.FE2, file = "_output/va.iviso.FE2.RDS")
va.iviso.FE2 <- readRDS(file = "_output/va.iviso.FE2.RDS")
write.csv(va.iviso.FE2$output, file = "_output/va.iviso.FE2.csv")

# ## Optimal GMM
# ### estimation of the fixed effects and E(Gy)
# EGy1         <- fEGy.iv(mydata, va.names, siv.FE1$coefficients[,"Estimate"], G, type = 1)
# EGy2         <- fEGy.iv(mydata, va.names, iv.FE2$coefficients, G, type = 2)
# mydata       <- mydata %>% select(!starts_with(c("EGy1", "EGy2"))) %>% mutate(EGy1 = EGy1$Egy, EGy2 = EGy2$Egy)
# mydataFE1    <- mydataFE1 %>% select(!starts_with("F1_EGy1")) %>% mutate(F1_EGy1 = peer.avg(F1, mydata$EGy1))
# mydataFE2    <- mydataFE2 %>% select(!starts_with("F2_EGy2")) %>% mutate(F2_EGy2 = peer.avg(F2, mydata$EGy2))
# 
# ### Model estimation
# ivopt.FE1    <- ivreg(formula = form.FE1, instruments = instropt.FE1, data = mydataFE1)
# (sivopt.FE1  <- summary(ivopt.FE1, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F1), function(x) rep(x, nrow(F1[[x]]))))))
# saveRDS(sivopt.FE1, file = "_output/sivopt.FE1.RDS")
# sivopt.FE1   <- readRDS(file = "_output/sivopt.FE1.RDS")
# 
# ivopt.FE2    <- ivreg(formula = form.FE2, instruments = instropt.FE2, data = mydataFE2)
# summary(ivopt.FE2, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F2), function(x) rep(x, nrow(F2[[x]]))))) 
# ml.ivopt.FE2 <- foptim(ivopt.FE2$residuals, ivopt.FE2$coefficients["F2G_gpa"], 
#                        G, fixed.effects = TRUE, F2, start = c(1, 0))
# va.ivopt.FE2 <- fvariance.iv(ivopt.FE2, ml.ivopt.FE2)  
# va.ivopt.FE2$output
# va.ivopt.FE2$sargan.stat
# va.ivopt.FE2$sargan.pvalue
# saveRDS(va.ivopt.FE2, file = "_output/va.ivopt.FE2.RDS")
# va.ivopt.FE2 <- readRDS(file = "_output/va.ivopt.FE2.RDS")

#' Identifying Peer Effects with Unobserved Effort and Isolated Students
#' Aristide Houndetoungan, Cristelle Kouame, and Michael Vlassopoulos
#' 
#' This file replicates the peer effect model estimation assuming that the network is exogenous
# set dvar
dvar  <- "gpa"
rm(list = ls()[ls() != "dvar"]) 
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

load(file = paste0("../../../Data/AHdata/PEEffort/AHD", dvar, ".rda"))
rm("Xlogit")
gc()

mydata        <- mydata %>% mutate(cst = 1)
va.names.cst  <- c("cst", va.names)
va.exo        <- va.names[-length(va.names)]
va.exo.cst    <- c("cst", va.exo)

# Econometric part
# data
G             <- norm.network(G)
XY            <- mydata[,va.names]
GXY           <- peer.avg(G, XY); colnames(GXY) <- paste0("G_", colnames(XY))
GGX           <- peer.avg(G, GXY[,-ncol(GXY)]); colnames(GGX) <- paste0("G", colnames(GXY[,-ncol(GXY)]))
hasfriends    <- unlist(lapply(G, rowSums))
Ghasfriends   <- peer.avg(G, hasfriends)

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

F1XY          <- peer.avg(F1, cbind(XY, hasfriends)); colnames(F1XY) <- paste0("F1_", c(colnames(XY), "hasfriends"))
F3XY          <- peer.avg(F3, XY); colnames(F3XY) <- paste0("F3_", colnames(XY))
F1GXY         <- peer.avg(F1, cbind(GXY, Ghasfriends)); colnames(F1GXY) <- paste0("F1", c(colnames(GXY), "G_hasfriends"))
F3GXY         <- peer.avg(F3, GXY); colnames(F3GXY) <- paste0("F3", colnames(GXY))
F1GGX         <- peer.avg(F1, GGX); colnames(F1GGX) <- paste0("F1", colnames(GGX))
F3GGX         <- peer.avg(F3, GGX); colnames(F3GGX) <- paste0("F3", colnames(GGX))
mydata        <- cbind(mydata, GXY, GGX, "hasfriends" = hasfriends, "G_hasfriends" = Ghasfriends)
mydataFE1     <- data.frame(F1XY, F1GXY, F1GGX)
mydataFE3     <- data.frame(F3XY, F3GXY, F3GGX)

rm(list = c("GXY", "GGX", "hasfriends", "Ghasfriends", "F1XY", "F3XY", "F1GXY", "F3GXY", "F1GGX", "F3GGX", "J1", "J3"))
gc()

# GMM
## formula and instruments
### Without fixed effects and without dummy for isolated students
form.noFE1    <- as.formula(paste0(dvar, "~", paste(c(-1, paste0("G_", dvar), va.exo.cst, paste0("G_", va.exo)), collapse = "+")))
instr.noFE1   <- as.formula(paste0("~", paste(c(-1, va.exo.cst, paste0("G_", va.exo), paste0("GG_", va.exo)), collapse = "+")))
### Without fixed effects and with dummy for isolated students
form.noFE2    <- as.formula(paste0(dvar, "~", paste(c(-1, paste0("G_", dvar), "hasfriends", va.exo.cst, paste0("G_", va.exo)), collapse = "+")))
instr.noFE2   <- as.formula(paste0("~", paste(c(-1, va.exo.cst, paste0("G_", va.exo), paste0("GG_", va.exo)), collapse = "+"), "+ hasfriends + G_hasfriends"))
### With fixed effects and without dummy for isolated students
form.FE1      <- as.formula(paste0("F1_", dvar, " ~", paste(c(-1, paste0("F1G_", dvar), paste0("F1_", va.exo), paste0("F1G_", va.exo)), collapse = "+")))
instr.FE1     <- as.formula(paste0("~", paste(c(-1, paste0("F1_", va.exo), paste0("F1G_", va.exo), paste0("F1GG_", va.exo)), collapse = "+")))
### With fixed effects and with dummy for isolated students
form.FE2      <- as.formula(paste0("F1_", dvar, " ~", paste(c(-1, paste0("F1G_", dvar), "F1_hasfriends", paste0("F1_", va.exo), paste0("F1G_", va.exo)), collapse = "+")))
instr.FE2     <- as.formula(paste0("~", paste(c(-1, paste0("F1_", va.exo), paste0("F1G_", va.exo), paste0("F1GG_", va.exo)), collapse = "+"), "+ F1_hasfriends + F1G_hasfriends"))
### Without fixed effects and with dummy for isolated students per school
form.FE3      <- as.formula(paste0("F3_", dvar, " ~", paste(c(-1, paste0("F3G_", dvar), paste0("F3_", va.exo), paste0("F3G_", va.exo)), collapse = "+")))
instr.FE3     <- as.formula(paste0("~", paste(c(-1, paste0("F3_", va.exo), paste0("F3G_", va.exo), paste0("F3GG_", va.exo)), collapse = "+")))

## IV estimation without fixed effects
### Without fixed effects and without dummy for isolated students
iv.noFE1     <- ivreg(formula = form.noFE1, instruments = instr.noFE1, data = mydata)
(siv.noFE1   <- summary(iv.noFE1, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(G), function(x) rep(x, nrow(G[[x]]))))))
saveRDS(siv.noFE1, file = paste0("_output/siv.noFE1", dvar, ".RDS"))
siv.noFE1    <- readRDS(file = paste0("_output/siv.noFE1", dvar, ".RDS"))
mean(siv.noFE1$residuals^2)
write.csv(siv.noFE1$coefficients, file = paste0("_output/siv.noFE1", dvar, ".csv"))

### Without fixed effects and with dummy for isolated students
iv.noFE2     <- ivreg(formula = form.noFE2, instruments = instr.noFE2, data = mydata)
(siv.noFE2   <- summary(iv.noFE2, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(G), function(x) rep(x, nrow(G[[x]]))))))
saveRDS(iv.noFE2, file = paste0("_output/iv.noFE2", dvar, ".RDS"))
iv.noFE2     <- readRDS(file = paste0("_output/iv.noFE2", dvar, ".RDS"))
ml.iv.noFE2  <- foptim(iv.noFE2$residuals, iv.noFE2$coefficients[paste0("G_", dvar)], 
                       G, fixed.effects = FALSE, F2, start = c(2.491167, 0.6053381))
va.iv.noFE2  <- fvariance.iv(iv.noFE2, ml.iv.noFE2, dvar)  
va.iv.noFE2$output
va.iv.noFE2$sargan.stat
va.iv.noFE2$sargan.pvalue
va.iv.noFE2$sigma2eta
va.iv.noFE2$sigma2epsilon
va.iv.noFE2$rho
saveRDS(va.iv.noFE2, file = paste0("_output/va.iv.noFE2", dvar, ".RDS"))
va.iv.noFE2  <- readRDS(file = paste0("_output/va.iv.noFE2", dvar, ".RDS"))
write.csv(va.iv.noFE2$output, file = paste0("_output/va.iv.noFE2", dvar, ".csv"))

## IV estimation with fixed effects
### With fixed effects and without dummy for isolated students
iv.FE1       <- ivreg(formula = form.FE1, instruments = instr.FE1, data = mydataFE1)
(siv.FE1     <- summary(iv.FE1, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F1), function(x) rep(x, nrow(F1[[x]]))))))
saveRDS(siv.FE1, file = paste0("_output/siv.FE1", dvar, ".RDS"))
siv.FE1      <- readRDS(file = paste0("_output/siv.FE1", dvar, ".RDS"))
mean(siv.FE1$residuals^2)
write.csv(siv.FE1$coefficients, file = paste0("_output/siv.FE1", dvar, ".csv"))

### With fixed effects and with dummy for isolated students
iv.FE2       <- ivreg(formula = form.FE2, instruments = instr.FE2, data = mydataFE1)
(siv.FE2     <- summary(iv.FE2, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F1), function(x) rep(x, nrow(F1[[x]]))))))
mean(siv.FE2$residuals^2)
saveRDS(iv.FE2, file = paste0("_output/iv.FE2", dvar, ".RDS"))
iv.FE2       <- readRDS(file = paste0("_output/iv.FE2", dvar, ".RDS"))
ml.iv.FE2    <- foptim(iv.FE2$residuals, iv.FE2$coefficients[paste0("F1G_", dvar)], 
                       G, fixed.effects = TRUE, F1, start = c(2.491167, 0.6053381))
va.iv.FE2    <- fvariance.iv(iv.FE2, ml.iv.FE2, dvar)  
va.iv.FE2$output
va.iv.FE2$sargan.stat
va.iv.FE2$sargan.pvalue
va.iv.FE2$sigma2eta
va.iv.FE2$sigma2epsilon
va.iv.FE2$rho
saveRDS(va.iv.FE2, file = paste0("_output/va.iv.FE2", dvar, ".RDS"))
va.iv.FE2    <- readRDS(file = paste0("_output/va.iv.FE2", dvar, ".RDS"))
write.csv(va.iv.FE2$output, file = paste0("_output/va.iv.FE2", dvar, ".csv"))

### With fixed effects and with dummy for isolated students per school
iv.FE3       <- ivreg(formula = form.FE3, instruments = instr.FE3, data = mydataFE3)
(siv.FE3     <- summary(iv.FE3, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F3), function(x) rep(x, nrow(F3[[x]]))))))
saveRDS(iv.FE3, file = paste0("_output/iv.FE3", dvar, ".RDS"))
iv.FE3       <- readRDS(file = paste0("_output/iv.FE3", dvar, ".RDS"))
ml.iv.FE3    <- foptim(iv.FE3$residuals, iv.FE3$coefficients[paste0("F3G_", dvar)], 
                       G, fixed.effects = TRUE, F3, start = c(2.491167, 0.6053381))
va.iv.FE3    <- fvariance.iv(iv.FE3, ml.iv.FE3, dvar)  
va.iv.FE3$output
va.iv.FE3$sargan.stat
va.iv.FE3$sargan.pvalue
va.iv.FE3$sigma2eta
va.iv.FE3$sigma2epsilon
va.iv.FE3$rho
saveRDS(va.iv.FE3, file = paste0("_output/va.iv.FE3", dvar, ".RDS"))
va.iv.FE3    <- readRDS(file = paste0("_output/va.iv.FE3", dvar, ".RDS"))
write.csv(va.iv.FE3$output, file = paste0("_output/va.iv.FE3", dvar, ".csv"))


#########################################################################################################
#####################################
## IV estimation with fixed effects removing completely isolated students
rm(list = ls()[ls() != "dvar"]) 
Rcpp::sourceCpp("codefiles/SourceCpp.cpp")
source("codefiles/SourceR.R")
load(file = paste0("../../../Data/AHdata/PEEffort/AHD", dvar, ".rda"))
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
XY            <- mydata[,va.names]
GXY           <- peer.avg(G, XY); colnames(GXY) <- paste0("G_", colnames(XY))
GGX           <- peer.avg(G, GXY[,-ncol(GXY)]); colnames(GGX) <- paste0("G", colnames(GXY[,-ncol(GXY)]))
hasfriends    <- unlist(lapply(G, rowSums))
Ghasfriends   <- peer.avg(G, hasfriends)

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

F1XY          <- peer.avg(F1, cbind(XY, hasfriends)); colnames(F1XY) <- paste0("F1_", c(colnames(XY), "hasfriends"))
F3XY          <- peer.avg(F3, XY); colnames(F3XY) <- paste0("F3_", colnames(XY))
F1GXY         <- peer.avg(F1, cbind(GXY, Ghasfriends)); colnames(F1GXY) <- paste0("F1", c(colnames(GXY), "G_hasfriends"))
F3GXY         <- peer.avg(F3, GXY); colnames(F3GXY) <- paste0("F3", colnames(GXY))
F1GGX         <- peer.avg(F1, GGX); colnames(F1GGX) <- paste0("F1", colnames(GGX))
F3GGX         <- peer.avg(F3, GGX); colnames(F3GGX) <- paste0("F3", colnames(GGX))
mydata        <- cbind(mydata, GXY, GGX, "hasfriends" = hasfriends, "G_hasfriends" = Ghasfriends)
mydataFE1     <- data.frame(F1XY, F1GXY, F1GGX)
mydataFE3     <- data.frame(F3XY, F3GXY, F3GGX)

rm(list = c("GXY", "GGX", "hasfriends", "Ghasfriends", "F1XY", "F3XY", "F1GXY", "F3GXY", "F1GGX", "F3GGX", "J1", "J3"))
gc()

# GMM
## formula and instruments
### With fixed effects and without dummy for isolated students
form.FE1      <- as.formula(paste0("F1_", dvar, " ~", paste(c(-1, paste0("F1G_", dvar), paste0("F1_", va.exo), paste0("F1G_", va.exo)), collapse = "+")))
instr.FE1     <- as.formula(paste0("~", paste(c(-1, paste0("F1_", va.exo), paste0("F1G_", va.exo), paste0("F1GG_", va.exo)), collapse = "+")))
### With fixed effects and with dummy for isolated students
form.FE2      <- as.formula(paste0("F1_", dvar, " ~", paste(c(-1, paste0("F1G_", dvar), "F1_hasfriends", paste0("F1_", va.exo), paste0("F1G_", va.exo)), collapse = "+")))
instr.FE2     <- as.formula(paste0("~", paste(c(-1, paste0("F1_", va.exo), paste0("F1G_", va.exo), paste0("F1GG_", va.exo)), collapse = "+"), "+ F1_hasfriends + F1G_hasfriends"))
### Without fixed effects and with dummy for isolated students per school
form.FE3      <- as.formula(paste0("F3_", dvar, " ~", paste(c(-1, paste0("F3G_", dvar), paste0("F3_", va.exo), paste0("F3G_", va.exo)), collapse = "+")))
instr.FE3     <- as.formula(paste0("~", paste(c(-1, paste0("F3_", va.exo), paste0("F3G_", va.exo), paste0("F3GG_", va.exo)), collapse = "+")))

## IV estimation with fixed effects
### With fixed effects and without dummy for isolated students
iv.FE1       <- ivreg(formula = form.FE1, instruments = instr.FE1, data = mydataFE1)
(siv.FE1     <- summary(iv.FE1, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F1), function(x) rep(x, nrow(F1[[x]]))))))
saveRDS(siv.FE1, file = paste0("_output/isosiv.FE1", dvar, ".RDS"))
siv.FE1      <- readRDS(file = paste0("_output/isosiv.FE1", dvar, ".RDS"))
mean(siv.FE1$residuals^2)
write.csv(siv.FE1$coefficients, file = paste0("_output/isosiv.FE1", dvar, ".csv"))

### With fixed effects and with dummy for isolated students
iv.FE2       <- ivreg(formula = form.FE2, instruments = instr.FE2, data = mydataFE1)
(siv.FE2     <- summary(iv.FE2, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F1), function(x) rep(x, nrow(F1[[x]]))))))
mean(siv.FE2$residuals^2)
saveRDS(iv.FE2, file = paste0("_output/isoiv.FE2", dvar, ".RDS"))
iv.FE2       <- readRDS(file = paste0("_output/isoiv.FE2", dvar, ".RDS"))
ml.iv.FE2    <- foptim(iv.FE2$residuals, iv.FE2$coefficients[paste0("F1G_", dvar)], 
                       G, fixed.effects = TRUE, F1, start = c(1.73344, 0.1682905))
va.iv.FE2    <- fvariance.iv(iv.FE2, ml.iv.FE2, dvar)  

va.iv.FE2$output
va.iv.FE2$sargan.stat
va.iv.FE2$sargan.pvalue
va.iv.FE2$sigma2eta
va.iv.FE2$sigma2epsilon
va.iv.FE2$rho
saveRDS(va.iv.FE2, file = paste0("_output/isova.iv.FE2", dvar, ".RDS"))
va.iv.FE2    <- readRDS(file = paste0("_output/isova.iv.FE2", dvar, ".RDS"))
write.csv(va.iv.FE2$output, file = paste0("_output/isova.iv.FE2", dvar, ".csv"))

### With fixed effects and with dummy for isolated students per school
iv.FE3       <- ivreg(formula = form.FE3, instruments = instr.FE3, data = mydataFE3)
(siv.FE3     <- summary(iv.FE3, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F3), function(x) rep(x, nrow(F3[[x]]))))))
saveRDS(iv.FE3, file = paste0("_output/isoiv.FE3", dvar, ".RDS"))
iv.FE3       <- readRDS(file = paste0("_output/isoiv.FE3", dvar, ".RDS"))
ml.iv.FE3    <- foptim(iv.FE3$residuals, iv.FE3$coefficients[paste0("F3G_", dvar)], 
                       G, fixed.effects = TRUE, F3, start = c(1.73344, 0.1682905))
va.iv.FE3    <- fvariance.iv(iv.FE3, ml.iv.FE3, dvar)  
va.iv.FE3$output
va.iv.FE3$sargan.stat
va.iv.FE3$sargan.pvalue
va.iv.FE3$sigma2eta
va.iv.FE3$sigma2epsilon
va.iv.FE3$rho
saveRDS(va.iv.FE3, file = paste0("_output/isova.iv.FE3", dvar, ".RDS"))
va.iv.FE3    <- readRDS(file = paste0("_output/isova.iv.FE3", dvar, ".RDS"))
write.csv(va.iv.FE3$output, file = paste0("_output/isova.iv.FE3", dvar, ".csv"))


#########################################################################################################
#####################################
## IV estimation with fixed effects removing any isolated students
rm(list = ls()[ls() != "dvar"]) 
Rcpp::sourceCpp("codefiles/SourceCpp.cpp")
source("codefiles/SourceR.R")
load(file = paste0("../../../Data/AHdata/PEEffort/AHD", dvar, ".rda"))
rm("Xlogit")
gc()

mydata        <- mydata %>% mutate(cst = 1) 
va.names.cst  <- c("cst", va.names)
va.exo        <- va.names[-length(va.names)]
va.exo.cst    <- c("cst", va.exo)

### We remove any isolated students
nhafr         <- lapply(G, rowSums) # number of friends
keep          <- lapply(1:nsch, function(x) nhafr[[x]] > 0)
mydata        <- mydata %>% filter(unlist(keep)) %>% group_by(sschlcde) %>% mutate(nstudent = n()) %>% ungroup() %>% filter(nstudent > 2)
G             <- lapply(1:nsch, function(x) G[[x]][keep[[x]], keep[[x]], drop = FALSE])
G             <- norm.network(G[sapply(G, is.matrix)])
Gnrow         <- sapply(G, nrow); G <- G[Gnrow > 2] # We remove school with completely isolated students
nsch          <- length(G)
sch.size      <- sapply(G, nrow)

### Data
XY            <- mydata[,va.names]
GXY           <- peer.avg(G, XY); colnames(GXY) <- paste0("G_", colnames(XY))
GGX           <- peer.avg(G, GXY[,-ncol(GXY)]); colnames(GGX) <- paste0("G", colnames(GXY[,-ncol(GXY)]))

J1            <- lapply(1:nsch, function(x) diag(sch.size[x]) - matrix(1, sch.size[x], sch.size[x])/sch.size[x])
F1            <- fdataFs(J1)

F1XY          <- peer.avg(F1, XY); colnames(F1XY) <- paste0("F1_", c(colnames(XY)))
F1GXY         <- peer.avg(F1, GXY); colnames(F1GXY) <- paste0("F1", c(colnames(GXY)))
F1GGX         <- peer.avg(F1, GGX); colnames(F1GGX) <- paste0("F1", colnames(GGX))
mydata        <- cbind(mydata, GXY, GGX)
mydataFE1     <- data.frame(F1XY, F1GXY, F1GGX)

rm(list = c("GXY", "GGX", "F1XY", "F1GXY", "F1GGX", "J1"))
gc()

# GMM
## formula and instruments
### With fixed effects and without dummy for isolated students
form.FE1      <- as.formula(paste0("F1_", dvar, " ~", paste(c(-1, paste0("F1G_", dvar), paste0("F1_", va.exo), paste0("F1G_", va.exo)), collapse = "+")))
instr.FE1     <- as.formula(paste0("~", paste(c(-1, paste0("F1_", va.exo), paste0("F1G_", va.exo), paste0("F1GG_", va.exo)), collapse = "+")))

## IV estimation with fixed effects
### With fixed effects and without dummy for isolated students
iv.FE1       <- ivreg(formula = form.FE1, instruments = instr.FE1, data = mydataFE1)
(siv.FE1     <- summary(iv.FE1, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F1), function(x) rep(x, nrow(F1[[x]]))))))
saveRDS(siv.FE1, file = paste0("_output/anyisosiv.FE1", dvar, ".RDS"))
siv.FE1      <- readRDS(file = paste0("_output/anyisosiv.FE1", dvar, ".RDS"))
mean(siv.FE1$residuals^2)
write.csv(siv.FE1$coefficients, file = paste0("_output/anyisosiv.FE1", dvar, ".csv"))

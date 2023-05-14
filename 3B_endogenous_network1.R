#' Identifying peer effects on academic achievements through students’ effort
#' Elysée Aristide Houndetoungan and Cristelle Kouame
#' 
#' This file replicates the peer effect model estimation controlling for network endogeneity.

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

# Network
G             <- norm.network(G)

# matrices J and F
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

rm(list = c("J1", "J2"))
gc()

festim        <- function(k, eff = "fe", interac = FALSE){
  data        <- readRDS(file = paste0("_output/endogeneity/data.np", k, ".RDS"))
  if((eff == "re") & interac){
    data      <- data*1e4 # has no influence on the result. We did this just because the values are low
  }
  data        <- cbind(mydata, data)
  
  newvar      <- c(paste0("muout.", eff, ".k", k, ".b.", 1:(k + 3)), paste0("muin.", eff, ".k", k, ".b.", 1:(k + 3)))
  if(interac){
    newvar    <- c(newvar, paste0("muinout.", eff, ".k", k, ".b.", sapply(1:(k + 3), function(x) sapply(1:(k + 3), function(y) paste0(x, ".", y)))))
  }
  va.exo      <- va.names[-length(va.names)]
  va.names.n  <- c(va.names, newvar)
  va.exo.n    <- c(va.exo, newvar)
  iInst       <- 1:length(va.exo)
  
  form.FE1    <- as.formula(paste("F1_gpa ~", paste(c(-1, paste0("F1_", va.exo), paste0("F1G_", va.exo), paste0("F1_", newvar)), collapse = "+"), "+ F1G_gpa"))
  instr.FE1   <- as.formula(paste("~", paste(c(-1, paste0("F1_", va.exo.n), paste0("F1G_", va.exo), paste0("F1GG_", va.exo)), collapse = "+")))
  form.FE2    <- as.formula(paste("F2_gpa ~", paste(c(-1, paste0("F2_", va.exo), paste0("F2G_", va.exo), paste0("F2_", newvar), paste0("F2G_", newvar)), collapse = "+"), "+ F2G_gpa"))
  instr.FE2   <- as.formula(paste("~", paste(c(-1, paste0("F2_", va.exo.n), paste0("F2G_", va.exo.n), paste0("F2GG_", va.exo)), collapse = "+")))
  
  GXY         <- peer.avg(G, data[,va.names.n]); colnames(GXY) <- paste0("G_", va.names.n)
  GGX         <- peer.avg(G, GXY[,iInst]); colnames(GGX) <- paste0("G", colnames(GXY[,iInst]))
  F1XY        <- peer.avg(F1, data[,va.names.n]); colnames(F1XY) <- paste0("F1_", va.names.n)
  F2XY        <- peer.avg(F2, data[,va.names.n]); colnames(F2XY) <- paste0("F2_", va.names.n)
  F1GXY       <- peer.avg(F1, GXY); colnames(F1GXY) <- paste0("F1", colnames(GXY))
  F2GXY       <- peer.avg(F2, GXY); colnames(F2GXY) <- paste0("F2", colnames(GXY))
  F1GGX       <- peer.avg(F1, GGX); colnames(F1GGX) <- paste0("F1", colnames(GGX))
  F2GGX       <- peer.avg(F2, GGX); colnames(F2GGX) <- paste0("F2", colnames(GGX))
  mydataFE1   <- data.frame(F1XY, F1GXY, F1GGX)
  mydataFE2   <- data.frame(F2XY, F2GXY, F2GGX)
  rm(list = c("GXY", "GGX", "F1XY", "F2XY", "F1GXY", "F2GXY", "F1GGX", "F2GGX", "data"))
  gc()
  
  # estimation
  iv.FE1      <- ivreg(formula = form.FE1, instruments = instr.FE1, data = mydataFE1)
  tryCatch({# May return error because of multicollinearity, so we tryCatch
    siv.FE1   <- summary(iv.FE1, diagnostics = TRUE, vcov = vcovCL, cluster = unlist(lapply(1:length(F1), function(x) rep(x, nrow(F1[[x]])))))
    endo.t1   <- ftestendo(iv.FE1$coefficients, siv.FE1$vcov)
    saveRDS(c(siv.FE1[c("coefficients", "sigma", "r.squared","adj.r.squared", "vcov", "diagnostics")], endo.t1), 
            file = paste0("_output/endogeneity/ed.siv.", eff, ifelse(interac, ".in", ".ad"), k,  ".FE1.RDS"))
  }, error = function(k) {print("summary(iv.FE2) failed")})
  
  iv.FE2      <- ivreg(formula = form.FE2, instruments = instr.FE2, data = mydataFE2)
  siv.FE2     <- list()
  tryCatch({# May return error because of multicollinearity, so we tryCatch
    siv.FE2   <- summary(iv.FE2, diagnostics = TRUE)
  }, error = function(k) {print("summary(iv.FE2) failed")})
  
  ml.iv.FE2   <- foptim(iv.FE2$residuals, iv.FE2$coefficients["F2G_gpa"], 
                        G, fixed.effects = TRUE, F2, start = c(2.26, 0.54))
  va.iv.FE2   <- fvariance.iv(iv.FE2, ml.iv.FE2); va.iv.FE2$Omega <- NULL
  va.iv.FE2$output
  va.iv.FE2$sargan.stat
  va.iv.FE2$sargan.pvalue
  va.iv.FE2$sigma2eta
  va.iv.FE2$sigma2epsilon
  va.iv.FE2$rho
  endo.t2     <- ftestendo(iv.FE2$coefficients, va.iv.FE2$var)
  saveRDS(c(siv.FE2[c("r.squared","adj.r.squared", "diagnostics")], va.iv.FE2, endo.t2), 
          file = paste0("_output/endogeneity/ed.va.iv", eff, ifelse(interac, ".in", ".ad"), k,  ".FE2.RDS"))
}


lapply(0:10, function(k) festim(k, eff = "fe", interac = FALSE))
lapply(0:10, function(k) festim(k, eff = "fe", interac = TRUE))
lapply(0:10, function(k) festim(k, eff = "re", interac = FALSE))
lapply(0:10, function(k) festim(k, eff = "re", interac = TRUE))

### subtract peer effect
siv.FE1      <- readRDS(file = "_output/siv.FE1.RDS")
va.iv.FE2    <- readRDS(file = "_output/va.iv.FE2.RDS")
PE0          <- data.frame(k       = rep(-1, 2),
                           coef    = c(siv.FE1$coefficients["F1G_gpa", "Estimate"], va.iv.FE2$output["F2G_gpa", "coef"]), 
                           SdtErr  = c(siv.FE1$coefficients["F1G_gpa", "Std. Error"], va.iv.FE2$output["F2G_gpa", "sd"]), 
                           effect  = rep("fixed effects", 2),
                           interac = "",
                           model   = c(1, 2))

PE           <- data.frame(k       = rep(NA, 88),
                           coef    = NA, 
                           SdtErr  = NA, 
                           effect  = NA,
                           interac = NA,
                           model   = NA)
j            <- 0
for (k in 0:10) {
  for (eff in c("fe", "re")) {
    for (interac in c(FALSE, TRUE)) {
      j               <- j + 1
      file            <- paste0("_output/endogeneity/ed.siv.", eff, ifelse(interac, ".in", ".ad"), k,  ".FE1.RDS")
      if(file.exists(file)){
        FE1           <- readRDS(file = file)
        PE$coef[j]    <- FE1$coefficients["F1G_gpa", "Estimate"]
        PE$SdtErr[j]  <- FE1$coefficients["F1G_gpa", "Std. Error"]
      }
      PE$model[j]     <- 3
      PE$k[j]         <- k
      PE$effect[j]    <- ifelse(eff == "fe", "fixed effects", "random effects")
      PE$interac[j]   <- ifelse(interac, 2, 1)
      j               <- j + 1
      file            <- paste0("_output/endogeneity/ed.va.iv", eff, ifelse(interac, ".in", ".ad"), k,  ".FE2.RDS")
      if(file.exists(file)){
        FE2           <- readRDS(file = file)
        PE$coef[j]    <- FE2$output["F2G_gpa", "coef"]
        PE$SdtErr[j]  <- FE2$output["F2G_gpa", "sd"]
      }
      PE$model[j]     <- 4
      PE$k[j]         <- k
      PE$effect[j]    <- ifelse(eff == "fe", "fixed effects", "random effects")
      PE$interac[j]   <- ifelse(interac, 2, 1)
    }
  }
}

PE                    <- rbind(PE0, PE)
PE$k                  <- factor(PE$k, labels = c("", 0:10))
PE$interac            <- factor(PE$interac, labels = c("", "B-splines", "Tensor products of B-splines"))
PE$model              <- factor(PE$model, labels = c("Standard model: Exo", "Proposed model: Exo", "Standard model: Endo", "Proposed model: Endo"))

library(ggplot2)
ggplot(PE %>% filter(effect == "fixed effects"), aes(x = k, colour = model)) + 
  geom_errorbar(width=.4, aes(ymin = coef - qnorm(0.975)*SdtErr, ymax = coef + qnorm(0.975)*SdtErr)) +
  geom_point(aes(y = coef, shape = model)) + theme_bw() +
  facet_grid(~ interac, scales = "free_x", shrink = TRUE, space='free') + 
  xlab("Number of knots") + ylab("Peer effect estimate") + 
  theme(legend.title = element_blank(), legend.position = "bottom")

# 8*3

# Export tables
# standard model
k            <- 10
eff          <- "re"
interac      <- FALSE
ED.siv.FE1   <- readRDS(file = paste0("_output/endogeneity/ed.siv.", eff, ifelse(interac, ".in", ".ad"), k,  ".FE1.RDS"))
write.csv(ED.siv.FE1$coefficients, file = paste0("_output/ed.siv.", eff, ifelse(interac, ".in", ".ad"), k,  ".FE1.csv"))
ED.siv.FE1$sigma^2
ED.siv.FE1$diagnostics
ED.siv.FE1$ed.wald.test
ED.siv.FE1$ed.wald.pvalue


# our approach
k            <- 10
eff          <- "re"
interac      <- FALSE
ED.iv.FE2    <- readRDS(file = paste0("_output/endogeneity/ed.va.iv", eff, ifelse(interac, ".in", ".ad"), k,  ".FE2.RDS"))
write.csv(ED.iv.FE2$output, file = paste0("_output/ed.va.iv", eff, ifelse(interac, ".in", ".ad"), k,  ".FE2.csv"))
ED.iv.FE2$output
ED.iv.FE2$sigma2eta
ED.iv.FE2$sigma2epsilon
ED.iv.FE2$rho
ED.iv.FE2$diagnostics
ED.iv.FE2$ed.wald.test
ED.iv.FE2$ed.wald.pvalue
ED.iv.FE2$sargan.stat
ED.iv.FE2$sargan.pvalue

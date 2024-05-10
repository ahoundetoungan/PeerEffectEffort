foptim <- function(resids, 
                   lambda, 
                   network, 
                   fixed.effects, 
                   Fs, 
                   start = c(1, 0), 
                   ctr = list(control = list(maxit = 1e6, abstol = 1e-9, reltol = 1e-9))){
  S    <- length(network)
  W    <- list()
  WW   <- list()
  I    <- list()
  nst  <- c()
  res  <- c()
  sumn <- NULL
  if(fixed.effects){
    tmp       <- ftoolsml(resids, Fs, network, lambda, S)
    nst       <- c(tmp$n)
    I         <- tmp$I
    W         <- tmp$W
    WW        <- tmp$WW
    res       <- unlist(tmp$lresid)
    sumn      <- sum(sapply(network, nrow))
  } else{
    nst       <- sapply(network, nrow)
    I         <- lapply(nst, diag)
    W         <- lapply(1:S, function(s) I[[s]] - lambda*network[[s]])
    WW        <- lapply(W, tcrossprod)
    res       <- resids
    sumn      <- sum(nst)
  }

  ctr1        <- c(ctr, list(par = c(log(start[1]), log((1 + start[2])/(1 - start[2]))), fn = floglike, 
                             resids = res, W = W, WW = WW, I = I, nst = nst, sumn = sumn, S = S))
  opt         <- do.call(optim, ctr1)
  
  tau         <- exp(opt$par[1])
  rho         <- (exp(opt$par[2]) - 1)/(exp(opt$par[2]) + 1)
  se2eps      <- fsigma2eps(tau, rho, res, W, WW, I, nst, sumn, S)
  Omega       <- se2eps$Om
  se2eps      <- se2eps$se2
  se2eta      <- tau*tau*se2eps
  c(opt, list(sigma2epsilon = se2eps, sigma2eta = se2eta, rho = rho, tau = tau, Omega = Omega))
}


fvariance.ols  <- function(ols.est, ml.est, dvar){
  mod          <- ols.est$model
  coef         <- ols.est$coefficients
  rem          <- names(coef)[which(is.na(coef))]
  coef         <- coef[!(names(coef) %in% rem)]
  mod          <- mod[, !(colnames(mod) %in% rem), drop = FALSE]
  K            <- length(coef)
  X            <- as.matrix(mod[,3:(K + 1)])
  R            <- as.matrix(cbind("Gy" = mod[,ifelse(paste0("G_", dvar) %in% colnames(mod), paste0("G_", dvar), 
                                                     ifelse(paste0("F1G_", dvar) %in% colnames(mod), paste0("F1G_", dvar), paste0("F3G_", dvar)))],
                                  X))
  Om           <- ml.est$Omega
  B0           <- crossprod(R)
  D0           <- crossprod(R, PartialNetwork::peer.avg(Om, R))
  vartheta     <- ml.est$sigma2epsilon*solve(B0, t(solve(t(B0), t(D0))))
  output       <- data.frame(variables = names(coef), coef = coef, sd = sqrt(diag(vartheta))) %>% mutate(t = coef/sd, prob = 2*(1 - pnorm(abs(t))))
  c(list(output = output, var = vartheta), ml.est)
}


fvariance.iv   <- function(iv.est, ml.est, dvar){
  mod          <- iv.est$model
  coef         <- iv.est$coefficients
  rem          <- names(coef)[which(is.na(coef))]
  coef         <- coef[!(names(coef) %in% rem)]
  mod          <- mod[, !(colnames(mod) %in% rem), drop = FALSE]
  K            <- length(coef)
  X            <- as.matrix(mod[,3:(K + 1)])
  R            <- as.matrix(cbind("Gy" = mod[,ifelse(paste0("G_", dvar) %in% colnames(mod), paste0("G_", dvar), 
                                                     ifelse(paste0("F1G_", dvar) %in% colnames(mod), paste0("F1G_", dvar), paste0("F3G_", dvar)))],
                                  X))
  Z            <- as.matrix(cbind(X, mod[, (K + 2):ncol(mod)]))
  Kz           <- ncol(Z)
  RZ           <- crossprod(R, Z)
  iZZRZ        <- solve(crossprod(Z), t(RZ), tol = 1e-30)
  B0           <- RZ %*% iZZRZ
  
  S            <- length(ml.est$Omega)
  nvec         <- sapply(ml.est$Omega, nrow)
  cumn         <- c(0, cumsum(nvec))
  Omch         <- matrix(0, Kz, Kz)
  for(s in 1:S){
    cat("Compute check(Omega) -- group: ", s, "/", S, "\n")
    Zs         <- Z[(cumn[s] + 1):cumn[s + 1], , drop = FALSE]
    Omch       <- Omch + crossprod(Zs, ml.est$Omega[[s]] %*% Zs)
  }
  
  D0           <- crossprod(iZZRZ, Omch %*% iZZRZ)
  
  ZOmZ         <- crossprod(Z, PartialNetwork::peer.avg(ml.est$Omega, Z))
  eZ           <- crossprod(iv.est$residuals, Z)
                                           
  sar.stat     <- c(eZ %*% solve(ZOmZ, t(eZ), tol = 1e-30)/ml.est$sigma2epsilon)
  sar.pva      <- 1 - pchisq(sar.stat, df = ncol(Z) - K)
    
  vartheta     <- ml.est$sigma2epsilon*solve(B0, t(solve(t(B0), t(D0), tol = 1e-30)), tol = 1e-30)
  output       <- data.frame(variables = names(coef), coef = coef, sd = sqrt(diag(vartheta))) %>% mutate(t = coef/sd, prob = 2*(1 - pnorm(abs(t))))
  c(list(output = output, var = vartheta, sargan.stat = sar.stat, sargan.pvalue = sar.pva), ml.est)
}


ftestendo      <- function(coef, var.cov){
  coef         <- coef[!is.na(coef)]
  su           <- grepl("_mu", names(coef))
  s.coef       <- coef[su]
  s.var.cov    <- var.cov[su, su]
  npar         <- length(s.coef)
  Test         <- sum(s.coef*solve(s.var.cov, s.coef, tol = 1e-30))
  Pvalue       <- 1 - pchisq(Test, df = npar)
  list("ed.wald.test" = Test, "ed.wald.pvalue" = Pvalue)
}

# fHausmanIV2    <- function(sexo, exo.lm, sendo){
#   mod          <- sexo$model
#   dmod         <- sendo$model
#   coef         <- sexo$coefficients
#   dcoef        <- sendo$coefficients
#   K            <- length(coef)
#   X            <- as.matrix(mod[, 2:K])
#   dX           <- as.matrix(dmod[, 2:K])
#   R            <- as.matrix(cbind(X, "G_gpa" = mod[, ifelse("G_gpa" %in% colnames(mod), "G_gpa", "F2G_gpa")]))
#   dR           <- as.matrix(cbind(dX, "G_gpa" = dmod[, ifelse("G_gpa" %in% colnames(dmod), "G_gpa", "F2G_gpa")]))
#   Z            <- as.matrix(cbind(X, mod[, (K + 2):ncol(mod)]))
#   dZ           <- as.matrix(cbind(dX, dmod[, (K + 2):ncol(dmod)]))
#   Kz           <- ncol(Z)
#   RZ           <- crossprod(R, Z)
#   dRZ          <- crossprod(dR, dZ)
#   iZZRZ        <- solve(crossprod(Z), t(RZ))
#   diZZRZ       <- solve(crossprod(dZ), t(dRZ))
#   B0           <- RZ %*% iZZRZ
#   dB0          <- dRZ %*% diZZRZ
#   
#   S            <- length(exo.lm$Omega)
#   nvec         <- sapply(exo.lm$Omega, nrow)
#   cumn         <- c(0, cumsum(nvec))
#   Omch         <- matrix(0, Kz, Kz)  #Z' x Omega x Z
#   DOmch        <- matrix(0, Kz, Kz)  #dZ'x Omega x dZ
#   dOmch        <- matrix(0, Kz, Kz)  #dZ'x Omega x  Z
#   for(s in 1:S){
#     cat("Compute check(Omega) -- group: ", s, "/", S, "\n")
#     Zs         <- as.matrix(mod[(cumn[s] + 1):cumn[s + 1], 2:K])
#     dZs        <- as.matrix(dmod[(cumn[s] + 1):cumn[s + 1], 2:K])
#     Zs         <- as.matrix(cbind(Zs, mod[(cumn[s] + 1):cumn[s + 1], (K + 2):ncol(mod)]))
#     dZs        <- as.matrix(cbind(dZs, dmod[(cumn[s] + 1):cumn[s + 1], (K + 2):ncol(dmod)]))
#     Omch       <- Omch + crossprod(Zs, exo.lm$Omega[[s]] %*% Zs)
#     dOmch      <- dOmch + crossprod(dZs, exo.lm$Omega[[s]] %*% Zs)
#     DOmch      <- DOmch + crossprod(dZs, exo.lm$Omega[[s]] %*% dZs)
#   }
#   
#   D0           <- crossprod(iZZRZ, Omch %*% iZZRZ)
#   dD0          <- crossprod(diZZRZ, dOmch %*% iZZRZ)
#   DD0          <- crossprod(diZZRZ, dOmch %*% diZZRZ)
#   
#   Bcom         <- as.matrix(bdiag(B0, dB0))
#   Dcomp        <- cbind(D0, t(dD0))
#   Dcomp        <- as.matrix(rbind(Dcomp, cbind(dD0, DD0)))
#   
#   vartheta     <- exo.lm$sigma2epsilon*solve(Bcom, t(solve(t(Bcom), t(Dcomp))))
#   Rmat         <- as.matrix(cbind(diag(K), -diag(K)))
#   RmatCoef     <- Rmat %*% c(coef, dcoef)
#   Haus         <- c(t(RmatCoef) %*% solve(Rmat %*% vartheta %*% t(Rmat), RmatCoef))
#   Prob         <- (1 - pchisq(Haus, K))
#   output       <- data.frame(variables = names(coef), coef = coef, sd = sqrt(diag(vartheta))) %>% mutate(t = coef/sd, prob = 2*(1 - pnorm(abs(t))))
#   list(variance = vartheta, test = Haus, p.value = Prob)
# }



# fEGy.iv        <- function(data, va.names, coefs, G, type = 2){
#   stopifnot(type %in% 1:2)
#   n.va.names   <- length(va.names)
#   va.exo       <- va.names[-n.va.names]
#   coefexo      <- head(coefs, length(coefs) - 1)
#   lambda       <- tail(coefs, 1)
#   Xb           <- c(as.matrix(data[,c(va.exo, paste0("G_", va.exo))]) %*% coefexo)
#   hattildv     <- data$gpa - lambda*data$G_gpa - Xb
#   nvec         <- sapply(G, nrow)
#   cumn         <- c(0, cumsum(nvec))
#   hasf         <- data$hasfriends
#   S            <- length(G)
#   FE           <- c()
#   
#   if(type == 1){
#     for (s in 1:S) {
#       vs       <- hattildv[(cumn[s] + 1):cumn[s + 1]]
#       FE       <- c(FE, rep(mean(vs), nvec[s]))
#     }
#   } else{
#     for (s in 1:S) {
#       vs       <- hattildv[(cumn[s] + 1):cumn[s + 1]]
#       hasfs    <- hasf[(cumn[s] + 1):cumn[s + 1]]
#       FE0      <- mean(vs[hasfs == 0])
#       FE1      <- mean(vs[hasfs == 1])
#       FE       <- c(FE, ifelse(hasfs == 0, FE0, FE1))
#     }
#   }
#   tmp          <- Xb #+ FE
#   Egy          <- unlist(lapply(1:S, function(s) solve(diag(nvec[s]) - lambda*G[[s]], tmp[(cumn[s] + 1):cumn[s + 1]])))
#   list(FE = FE, Egy = Egy)
# }
# foptimParallel <- function(resids, 
#                            lambda, 
#                            network, 
#                            fixed.effects, 
#                            J, 
#                            start  = c(0, 0), 
#                            ctr    = list(control = list(maxit = 1e6, abstol = 1e-13, reltol = 1e-13)),
#                            ncores = parallel::detectCores(all.tests = FALSE, logical = TRUE)){
#   S    <- length(network)
#   W    <- list()
#   WW   <- list()
#   I    <- list()
#   n    <- c()
#   res  <- list()
#   if(fixed.effects){
#     nvec       <- sapply(network, nrow)
#     ncum       <- c(0, cumsum(nvec))
#     for (s in 1:S) {
#       cat("Find Eigenvalues and Eigenvectors: ", s, "/", S, "\n")
#       es       <- resids[(ncum[s] + 1):ncum[s + 1]]
#       eig      <- eigen(J[[s]])
#       Fs       <- eig$vectors[,which(round(eig$values, 7) == 1)]
#       n[s]     <- ncol(Fs)
#       I[[s]]   <- diag(n[s])
#       W[[s]]   <- crossprod(Fs, I[[s]] - lambda*network[[s]])
#       WW[[s]]  <- tcrossprod(W[[s]])
#       W[[s]]   <-  W[[s]]%*%Fs
#       res[[s]] <- crossprod(Fs, es)
#     }
#   } else{
#     n          <- sapply(network, nrow)
#     ncum       <- c(0, cumsum(n))
#     I          <- lapply(n, diag)
#     W          <- lapply(1:S, function(s) I[[s]] - lambda*network[[s]])
#     WW         <- lapply(W, tcrossprod)
#     res        <- lapply(1:S, function(s) resids[(ncum[s] + 1):ncum[s + 1]])
#   }
#   sumn        <- sum(n)
#   ctr1        <- c(ctr, list(par = c(log(start[1]), log((1 + start[2])/(1 - start[2]))), fn = fllh, 
#                              lresids = res, W = W, WW = WW, I = I, n = n, sumn = sumn, S = S, ncores = ncores))
#   opt         <- do.call(optim, ctr1)
#   
#   tau         <- exp(opt$par[1])
#   rho         <- (exp(opt$par[2]) - 1)/(exp(opt$par[2]) + 1)
#   
#   se2eps      <- sum(unlist(parallel::mclapply(1:S, function(s) f1sig2(tau = tau, rho = rho, es = res[[s]], 
#                                                                        Ws = W[[s]], WWs = WW[[s]], Is = I[[s]]), mc.cores = ncores)))/sumn
#   se2eta      <- tau*tau*se2eps
#   c(opt, list(sigma2epsilon = se2eps, sigma2eta = se2eta, rho = rho, tau = tau))
# }
# 
# 
# fllh          <- function(x, lresids, W, WW, I, n, S, sumn, ncores){
#   tau         <- exp(x[1])
#   rho         <- (exp(x[2]) - 1)/(exp(x[2]) + 1)
#   cat("tau: ", tau, "\n")
#   cat("rho: ", rho, "\n")
#   
#   out1        <- parallel::mclapply(1:S, function(s) f1s(tau = tau, rho = rho, es = lresids[[s]], 
#                                                          Ws = W[[s]], WWs = WW[[s]], Is = I[[s]]), mc.cores = ncores)
#   sig2eps     <- sum(sapply(out1, function(o) o$se2s))/sumn
#   cat("sigma^2_epsilon: ", sig2eps, "\n")
#   lOm         <- lapply(out1, function(o) o$Oms)
#   
#   out2        <- parallel::mclapply(1:S, function(s) f2s(sig2eps, lOm[[s]], n[s]))
#   llh         <- sum(unlist(out2))
#   cat("log(likelihood): ", llh, "\n")
#   cat("*****************\n")
#   -llh
# }

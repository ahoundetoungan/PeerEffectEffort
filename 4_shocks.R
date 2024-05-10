# #' Identifying Peer Effects with Unobserved Effort and Isolated Students
#' Aristide Houndetoungan, Cristelle Kouame, and Michael Vlassopoulos
#'
#' This file simulates schocks to GPA
# set dvar
dvar  <- "gpa"
rm(list = ls()[ls() != "dvar"])
library(PartialNetwork)
library(AER)
library(dplyr)
library(ggplot2)

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

# Simulation functions
M       <- length(G)
nvec    <- sapply(G, nrow)
ncum    <- c(0, cumsum(nvec))
lNI     <- unlist(lapply(G, rowSums))
lI      <- 1 - lNI

fdya    <- function(lambda){
  dK    <- lI + (1 - lambda)*lNI
  unlist(lapply(1:M, function(x) solve(diag(nvec[x]) - lambda*G[[x]], dK[(ncum[x] + 1):ncum[x + 1]])))
}

fdyc    <- function(lambda){
  dK    <- lI + lNI
  unlist(lapply(1:M, function(x) solve(diag(nvec[x]) - lambda*G[[x]], dK[(ncum[x] + 1):ncum[x + 1]])))
}

# Ignoring network endogeneity
dy2     <- fdyc(0.507479404169283)
dy3     <- fdyc(0.751110035546685)
dy4a    <- fdya(0.855882248928115)
dy4c    <- fdyc(0.855882248928115)

# Controlling for network endogeneity
dy2e    <- fdyc(0.671899526930314)
dy3e    <- fdyc(0.728699298999939)
dy4ae   <- fdya(0.82788591602883)
dy4ce   <- fdyc(0.82788591602883)

dplot   <- data.frame(endo  = rep(0:1, each = 4*sum(nvec)),
                      shock = rep(rep(LETTERS[4:1], 2), each = sum(nvec)),
                      dy    = c(dy2, dy3, dy4a, dy4c, dy2e, dy3e, dy4ae, dy4ce)) %>%
  mutate(endo  = factor(endo, labels = c("Exogenous network", "Endogenous network")),
         dy1   = ifelse(dy == 1, 1, NA))

shockn  <- c(expression(paste("Model 4: Shock on ", c[s])),
             expression(paste("Model 4: Shock on ", alpha[s])),
             expression(paste("Model 3: Shock on ", c[s])),
             expression(paste("Model 2: Shock on ", c[s])))

(graph  <- ggplot(dplot %>% filter(endo == "Exogenous network"), aes(y = shock, x = dy, group = shock)) +
    geom_boxplot(outlier.shape = 1) +
    geom_point(aes(y = shock, x = dy1, group = shock), shape = 1) +
    ylab("") + xlab("GPA increase") + theme_bw() +
    # facet_wrap(~ endo, ncol = 2) +
    scale_y_discrete(labels = shockn) +
    scale_x_continuous(breaks = 1:7) +
    theme(strip.text = element_text(face = "italic"),
          text = element_text(size = 12, family = "Palatino"),
          axis.title = element_text(size = 12, family = "Palatino")))

ggsave("plot:shocks.pdf", path = "_output", plot = graph, device = "pdf", width = 7, height = 3)

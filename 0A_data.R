#' Identifying Peer Effects with Unobserved Effort and Isolated Students
#' Aristide Houndetoungan, Cristelle Kouame, and Michael Vlassopoulos
#' 
#' This file prepares the data set to be used.
#' The data set will be saved under the name "AHD.rda" and will be used in other files.
#' It also computes descriptive statistics

rm(list = ls())
library(readstata13)
library(ggplot2)
library(PartialNetwork)

proot <- c("~/GPAeffort",
           "~/Dropbox/Papers - In progress/EffortGPA/Code-EffortGPA")
root  <- sapply(proot, dir.exists)
root  <- proot[root][1]
setwd(root)


# the finale data set is save in the 'filname' path and can be loaded if saved before
# otherwise, the code will be ran to prepare the data
vdep           <- "gpa"
filname        <- paste0("../../../Data/AHdata/PEEffort/AHD", vdep, ".rda") 
if (file.exists(filname)) {
  load(file = filname)
} else {
  # Data
  mydata       <- read.dta13("../../../Data/AHdata/cleandta.dta")  # data from Stata
  mydata       <- mydata[order(mydata$sschlcde),] 
  mydata$club  <- ifelse(mydata$nclubs > 0, 1, 0)
  mydata$male  <- 1 - mydata$female
  
  mislist      <- c(55555555, 77777777, 88888888, 99999999, 99959995)
  mf_coln      <- paste0("mf", 1:5, "aid")
  ff_coln      <- paste0("ff", 1:5, "aid")
  f_coln       <- c(mf_coln, ff_coln)
  
  # mislist is an ID
  if (sum(mydata$aid %in% mislist) > 0) {
    stop("mislist is an ID")
  } else {
    cat("mislist is not an ID: OK", "\n")
  }
  
  # list of variable (excluding reference variables)
  va.all.names   <- c("male", "female", "age", "hispanic", "racewhite", "raceblack", "raceasian", "raceother", 
                      "withbothpar", "yearinschl", "club", "mehigh", "melhigh", "memhigh", "memiss", "mjprof", 
                      "mjhome", "mjother", "mjmiss", vdep)
  
  # list of variable (excluding reference variables for identification)
  va.names      <- c("female", "age", "hispanic", "raceblack", "raceasian", "raceother", 
                     "withbothpar", "yearinschl", "club", "melhigh", "memhigh", "memiss", "mjprof", 
                     "mjother", "mjmiss", vdep)
  
  # remove friend from different groups
  # remove self friendship
  # remove friend non found
  N       <- nrow(mydata)
  dscf    <- rep(0,N)
  sfre    <- rep(0,N)
  nffr    <- rep(0,N)
  for (i in 1:N) {
    for (j in f_coln) {
      k  <- which(mydata$aid == mydata[i, j, drop = TRUE])
      # remove if different school
      if (length(k) != 0) {
        if(mydata[i, "sschlcde", drop = TRUE] != mydata[k, "sschlcde", drop = TRUE]) {
          mydata[i, j]   <- -1
          dscf[i]        <- dscf[i] + 1
        }
        # remove if self frienship
        if(mydata[i, "aid", drop = TRUE] == mydata[k, "aid", drop = TRUE]) {
          mydata[i, j]   <- -2
          sfre[i]        <- sfre[i] + 1
        }
      }
      else {
        if (!((mydata[i, j, drop = TRUE] %in% mislist) | is.na(mydata[i, j, drop = TRUE]))) {
          mydata[i, j]   <- -3
          nffr[i]        <- nffr[i] + 1
        }
      }
    }
  }
  
  cat("remove", sum(dscf), "link(s) because students from different schools: their code are recode as -1", "\n")
  cat("remove", sum(sfre), "self-friendship(s): their code are recoded as -2", "\n")
  cat("remove", sum(nffr), "non-found friends: their code are recoded as -3", "\n")
  rm(list = c("i", "j", "k"))
  
  # Are there NAs?
  apply(mydata[,va.names], 2, function(w) sum(is.na(w)))
  
  # Keep row without NA
  keep1         <- as.logical((rowSums(is.na(mydata[,va.names])) == 0))
  mydata        <- mydata[keep1,]
  
  # schools
  schtable     <- table(mydata$sschlcde)
  school       <- as.numeric(names(schtable))
  sch.size     <- as.numeric(schtable)
  nsch         <- length(sch.size)
  
  # dependent variable
  mydata$y     <- mydata[,tail(va.names,1), drop = TRUE]
  
  mislistmis   <- c(55555555, 99999999, 99959995)
  
  # This function prepares the data and the network 
  gen.data  <- function(db) {
    G       <- vector("list", nsch)
    nmatch  <- vector("list", nsch)
    X       <- as.matrix(db[, va.names[-length(va.names)]])    
    Y       <- db$y
    
    
    # output for this loop
    # G, 
    # Gobs (the observed part of G)
    # X containing also number of missing links
    # ni is number of students in all schools
    for (i in 1:nsch) {
      cat("School :", i, "/", nsch, "\n")
      schi         <- school[i]
      dbi          <- db[db$sschlcde == schi,]
      Ni           <- nrow(dbi)
      Gi           <- matrix(0, Ni, Ni)
      
      nmatchi      <- numeric() #will contain missing links
      for (j in 1:Ni) {
        idx        <- which(dbi$aid %in% dbi[j, f_coln])
        Gi[j, idx] <- 1
        
        # missing links
        idx.miss   <- which(unlist(dbi[j, f_coln]) %in% c(mislistmis, -3))
        nmatchi[j] <- length(idx.miss) # Number of missing links
      }
      
      #  store G
      diag(Gi)     <- 0
      G[[i]]       <- Gi
      
      # unmatched
      nmatch[[i]]  <- nmatchi
    }
    
    list(G       = G,
         X       = X,
         nmatch  = nmatch,
         Y       = Y)
    
  }
  
  # Use the function to prepare the data
  tmp     <- gen.data(mydata)
  G       <- tmp$G
  X       <- tmp$X
  nmatch  <- tmp$nmatch
  Y       <- tmp$Y
  
  # dataset for the logistic model
  # # variable to include 0 (both male, both female, 0)
  # va.log0       <- "female"
  # variable to include 1 (indicator == 1 if same value)
  va.log1       <- c("female", "male", "hispanic", "racewhite", "raceblack", "raceasian", 
                     "club", "mehigh", "melhigh", "memhigh", "mjprof", "mjhome")
  # variable to include  2 (absolute value of the difference)
  va.log2       <- c("age", "yearinschl")
  
  # distance
  # dist0         <- function(x, y) ifelse(x == 1 & x == 1, 2, ifelse(x == 0 & x == 0, 1, 0))
  dist1         <- function(x, y) as.numeric(x == 1 & y == 1)
  dist2         <- function(x, y) abs(x - y)
  
  # data
  tmp           <- c(0, cumsum(sch.size))
  
  # for va.log1
  X1tmp         <- do.call("cbind", lapply(va.log1, function(z) {
    mat.to.vec(lapply(1:nsch, function(x) {
      va.zx     <- mydata[c(tmp[x] + 1):tmp[x+1],z]   
      matrix(kronecker(va.zx, va.zx, FUN = dist1), sch.size[x])}))
  }))
  # for va.log2
  X2tmp         <- do.call("cbind", lapply(va.log2, function(z) {
    mat.to.vec(lapply(1:nsch, function(x) {
      va.zx     <- mydata[c(tmp[x] + 1):tmp[x+1],z]   
      matrix(kronecker(va.zx, va.zx, FUN = dist2), sch.size[x])}))
  }))
  Xlogit            <- as.matrix(cbind(X1tmp, X2tmp))  
  colnames(Xlogit)  <- c(va.log1, paste0("diff.", va.log2))
  
  save(list = c("G", "dscf", "f_coln", "ff_coln", "filname", "mf_coln", "mislist", "mislistmis",
                "mydata", "nmatch", "nsch", "sch.size", "va.all.names", "va.log1", "va.log2", 
                "va.names", "Xlogit"), file = filname)
}
############################ Descriptive stat ##################################
rm("Xlogit")
gc()
# Descriptive statistic function 
library(dplyr)
my.stat.des    <- function(x) {
  out <- c(mean(x, na.rm = TRUE), 
           sd(x, na.rm = TRUE), 
           quantile(x, probs = c(0.005, 0.995, 0.025, 0.975, 0.05, 0.95, 0.25, 0.75), na.rm = TRUE))
  names(out)   <- c("Mean", "St. Dev.",   "Pctl(0.5)",   "Pctl(99.5)",   "Pctl(2.5)",   "Pctl(97.5)",   "Pctl(5)",   "Pctl(95)",   "Pctl(25)",   "Pctl(75)")
  out
}

# size
dim(mydata)
# 68430    54

# Number of schools and school's size
length(G)
quantile(sapply(G, nrow), prob = c(0, 0.25, 0.5, 0.75, 1))
# 0%  25%  50%  75% 100% 
# 18  216  368  605 2027
summary(sapply(G, nrow))
# Min.   1st Qu.  Median    Mean   3rd Qu.    Max. 
# 18.0   216.0     368.0   485.3   605.0    2027.0

# all the variables 
hnFr           <- which(unlist(lapply(G, rowSums)) == 0)
allVar         <- peer.avg(norm.network(G), mydata[,va.all.names])
allVar[hnFr,]  <- NA
allVar         <- cbind(mydata[,va.all.names], allVar)

# Descriptive stat
sdes           <- round(t(apply(allVar, 2, my.stat.des)), 3)[,c(1,2,9,10)]
sdes
print(sdes)
write.csv(sdes, file = paste0("_output/sdes", vdep, ".csv"))

# friends
cumsum         <- c(0, cumsum(sch.size))
nhafr          <- unlist(lapply(G, rowSums)) # number of friends
nisfr          <- unlist(lapply(G, colSums)) # number of times the student is a friend

nhafrmale      <- unlist(lapply(1:nsch, function(s) apply(G[[s]], 1, function(gi) sum(gi*mydata[(cumsum[s] + 1):cumsum[s + 1], "male"]))))
nhafrfemale    <- unlist(lapply(1:nsch, function(s) apply(G[[s]], 1, function(gi) sum(gi*mydata[(cumsum[s] + 1):cumsum[s + 1], "female"]))))
Friends        <- data.frame(gender      = ifelse(mydata$male == 1, "male", "female"),
                             hasfr       = nhafr,
                             isfr        = nisfr,     
                             hasfrmale   = nhafrmale,
                             hasfrfemale = nhafrfemale)
(sumFriends    <- Friends %>% group_by(gender) %>% summarise(across(c("hasfr", "hasfrmale", "hasfrfemale"), mean)) %>% bind_rows(
  Friends %>% summarise(across(c("hasfr", "hasfrmale", "hasfrfemale"), mean))))
write.csv(sumFriends, file = paste0("_output/FriendsByGender", vdep, ".csv"))
# gender   nfr  nffr  nmfr
# <chr>  <dbl> <dbl> <dbl>
# female  3.64  2.28  1.36
# male    3.24  1.41  1.83
# NA      3.44  1.85  1.59

# Descriptive stat for student who have no friends
sdes.nofriends <- round(t(apply(allVar[nhafr == 0,], 2, my.stat.des)), 3)[,c(1,2,9,10)]
sdes.nofriends
print(sdes.nofriends)
write.csv(sdes.nofriends, file = paste0("_output/sdes.nofriends", vdep, ".csv"))

# distribution of number of friends
disthasfr      <- Friends %>% group_by(hasfr) %>% summarise(count = length(hasfr))
print(as.data.frame(disthasfr), row.names = FALSE)
# hasfr count
# 0 14900
# 1  6402
# 2  7240
# 3  7889
# 4  7768
# 5  6886
# 6  5893
# 7  4920
# 8  3657
# 9  2185
# 10   690
write.csv(disthasfr, file = paste0("_output/distHasFriends", vdep, ".csv"))

# distribution of the number of times the student is a friends
distisfr        <- Friends %>% group_by(isfr) %>% summarise(count = length(isfr))
print(as.data.frame(distisfr), row.names = FALSE)
# isfr count
# 0 13161
# 1  9665
# 2  9777
# 3  8664
# 4  7041
# 5  5403
# 6  4132
# 7  2985
# 8  2238
# 9  1491
# 10  1112
# 11   727
# 12   598
# 13   419
# 14   311
# 15   206
# 16   156
# 17   109
# 18    75
# 19    48
# 20    31
# 21    23
# 22    19
# 23    18
# 24     6
# 25    10
# 27     2
# 28     1
# 29     1
# 35     1
write.csv(distisfr, file = paste0("_output/distIsFriends", vdep, ".csv"))

(tabFriends     <- table(Friends$isfr, Friends$hasfr))
# 0    1    2    3    4    5    6    7    8    9   10
# 0  7245 1604 1239 1001  735  546  330  240  139   65   17
# 1  2465 1639 1452 1297  975  704  490  327  209   96   11
# 2  1778 1211 1474 1463 1224  976  665  500  301  143   42
# 3  1195  774 1099 1282 1311 1001  764  606  355  222   55
# 4   769  461  720 1026 1009  992  774  616  407  210   57
# 5   473  294  475  660  824  757  675  554  416  214   61
# 6   361  141  283  419  575  592  569  489  421  215   67
# 7   199   96  194  266  391  402  467  427  289  187   67
# 8   136   72  128  183  250  307  364  305  263  173   57
# 9    95   40   56  109  151  187  231  237  203  140   42
# 10   59   15   41   69  118  139  177  173  164  120   37
# 11   32   21   30   35   67   83  107  115  118   80   39
# 12   31   12   10   28   49   71   79  104   94   85   35
# 13   20    9   17   22   20   49   68   53   69   67   25
# 14   17    5   11    7   27   26   48   53   66   34   17
# 15    9    2    2    6   13   20   25   30   41   43   15
# 16    5    1    5    5   15   12   20   24   32   17   20
# 17    5    3    1    3    6    5   13   17   27   21    8
# 18    0    0    1    2    3   10    5   15   18   14    7
# 19    1    1    0    1    2    0    4   14    7   13    5
# 20    1    0    0    1    1    2    6    7    4    7    2
# 21    2    0    0    0    2    1    6    3    4    5    0
# 22    0    1    2    2    0    1    1    5    1    5    1
# 23    0    0    0    0    0    1    4    3    4    4    2
# 24    0    0    0    1    0    1    0    0    1    2    1
# 25    2    0    0    1    0    0    0    2    4    1    0
# 27    0    0    0    0    0    0    0    0    0    2    0
# 28    0    0    0    0    0    1    0    0    0    0    0
# 29    0    0    0    0    0    0    0    1    0    0    0
# 35    0    0    0    0    0    0    1    0    0    0    0
write.csv(tabFriends, file = paste0("_output/tabFriends", vdep, ".csv"))

ggplot(data = data.frame(nfriends = factor(nhafr)), aes(x = nfriends)) +
  geom_bar(aes(y = (..count..)/sum(..count..)), color = "black", fill = "#eeeeee") + 
  theme_bw() + xlab("") + ylab("Frequency") + 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1L)) 
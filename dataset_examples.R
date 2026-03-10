library(netmeta)
library(poth)
library(stringr)
source("helperfuncs.R")

# Case study 1 -----------------------------

data("Linde2016")

# Conduct random effects network meta-analysis
#
net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
                data = Linde2016, ref = "placebo", sm = "OR", 
                common = FALSE, small.values = "undesirable")

# Additive model for treatment components
#
nc1 <- netcomb(net1)

# Design Matrix M

nc1$X.matrix

# Question 1: hierarchy of all observed treatments in terms of median
set.seed(2026)
linrank1 <- cnmaRank(nc1, set = nc1$trts, 
                     small.values = "undesirable", re = T,
                     iter = 1000,
                     metric = "medianrank")

linrank1

forest(nc1, rightcols = c("effect", "ci","value"), 
       rightlabs = c(NA, NA, "Median rank"),
       add.data = linrank1[order(row.names(linrank1)),],
       sortvar = -TE,
       just.addcols = "right",
       backtransf = FALSE)

# QUestion 2: Hierarchy of incremental effects in terms of probability of best value?

set.seed(2026)
linrank2 <- cnmaRank(nc1, set = c("Face-to-face PST",
                                  "Face-to-face interpsy",
                                  "Face-to-face CBT",
                                  "SSRI"), 
                     small.values = "undesirable", re = T,
                     iter = 1000,
                     metric = "Pbest")

linrank2
netcomparison(nc1, treat1 = c("Face-to-face PST",
                              "Face-to-face interpsy",
                              "Face-to-face CBT",
                              "SSRI"),
              treat2 = "Face-to-face CBT", backtransf = F)

# Case study 2 --------------------------------------------------
#Chronic Lymphocytic Leukemia data from Li et al supplementary info
Ex1 <- data.frame(id = c(2, 3, 1, 7, 4, 5, 6, 8, 9, 10),
                  n = c("261 (174/87)", "319 (160/159)", "389 (195/196)", "160(106/54)", "578 (289/289)", "389 (194/195)", "416 (207/209)", "220 (110/110)", "208 (104/104)", "117 (59/58)"),
                  t1 = c("Ide+Ofa", "Duv", "Ibr", "Ibr", "Ibr+Ben+Rit", "Ven+Rit", "Ide+Ben+Rit", "Ide+Rit", "Ibr", "Ubl+Ibr"),
                  t2 = c("Ofa", "Ofa", "Ofa", "Rit", "Ben+Rit", "Ben+Rit", "Ben+Rit", "Rit", "Ibr+Rit", "Ibr"),
                  lnpfshr = c(-1.34707364796661, -0.653926467406664, -2.26336437984076, -1.71479842809193, -1.59454929994035, -1.77195684193188, -1.10866262452161, -1.89711998488588, 0.150142658429719, -0.581605805827038),
                  selnpfshr = c(0.183812794578587, 0.149217754061151, 0.180294089952991, 0.274525365514299, 0.155552441740024, 0.197242318426909, 0.144212706390322, 0.319582389922288, 0.378985732075359, 0.484490089563701))

disEx1 <- discomb(lnpfshr, selnpfshr, t1, t2, id, reference.group = "Ofa",
                 sm="md", data = Ex1, small.values = "undesirable")

# Design matrix
M1 <- disEx1$X.matrix

# Some examples using checkIdentifiable function
checkIdentifiable(M1, v = c(1,0,0,0,0,0,0,-1))
checkIdentifiable(M1, v = c(0,0,0,0,1,0,0,-1))
checkIdentifiable(M1, v = c(1,0,0,0,0,-1,0,0))

# Question 1: Hierarchy of all observed treatments in terms of their point estimates
disrank1 <- cnmaRank(disEx1, set = disEx1$trts, verbose = T,
                  small.values = "undesirable", re = F,
                  metric = "pointestimate") # fails

# Question 2: Hierarchy of novel targeted agent components based on the expected rank
disrank2 <- cnmaRank(disEx1, set = c("Duv",
                                  "Ibr",
                                  "Ide",
                                  "Ubl",
                                  "Ven"), verbose = T,
                  small.values = "undesirable", re = F,
                  metric = "expectedrank")


# Revised question 2: Hierarchy of subset of novel targeted agent components based on E(rank)
set.seed(2026)
disrank2 <- cnmaRank(disEx1, set = c("Duv",
                                     "Ibr",
                                     "Ide",
                                     "Ubl"), verbose = F,
                     small.values = "undesirable", re = F,
                     iter = 1000,
                     metric = "expectedrank")

disrank2
netcomparison(disEx1, treat1 = c("Duv",
                                 "Ibr",
                                 "Ide",
                                 "Ubl"),
              treat2 = "Duv")


# Plotting ----------------------------------

othernet <- which(disEx1$trts %in% c("Ven+Rit",
                                     "Ide+Ben+Rit",
                                     "Ibr+Ben+Rit",
                                     "Ben+Rit"))

cols <- rep("navy", 12)
cols[othernet] <- "lightblue"

netgraph(disEx1, number.of.studies = F, points = T, cex.points = 3,
         col.points = cols, lwd = 1,
         main = "R/R CLL Network")

netgraph(net1, number.of.studies = F, 
         points = T, cex.points = 3,
         col.points = "cornflowerblue", 
         lwd = 1,
         main = "Depression Network")

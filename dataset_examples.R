library(netmeta)
library(poth)
library(stringr)
source("helperfuncs.R")

# Example 1 data from Li et al --------------------------------------------------
## Example 1: Chronic Lymphocytic Leukemia
Ex1 <- data.frame(id = c(2, 3, 1, 7, 4, 5, 6, 8, 9, 10),
                  n = c("261 (174/87)", "319 (160/159)", "389 (195/196)", "160(106/54)", "578 (289/289)", "389 (194/195)", "416 (207/209)", "220 (110/110)", "208 (104/104)", "117 (59/58)"),
                  t1 = c("Ide+Ofa", "Duv", "Ibr", "Ibr", "Ibr+Ben+Rit", "Ven+Rit", "Ide+Ben+Rit", "Ide+Rit", "Ibr", "Ubl+Ibr"),
                  t2 = c("Ofa", "Ofa", "Ofa", "Rit", "Ben+Rit", "Ben+Rit", "Ben+Rit", "Rit", "Ibr+Rit", "Ibr"),
                  lnpfshr = c(-1.34707364796661, -0.653926467406664, -2.26336437984076, -1.71479842809193, -1.59454929994035, -1.77195684193188, -1.10866262452161, -1.89711998488588, 0.150142658429719, -0.581605805827038),
                  selnpfshr = c(0.183812794578587, 0.149217754061151, 0.180294089952991, 0.274525365514299, 0.155552441740024, 0.197242318426909, 0.144212706390322, 0.319582389922288, 0.378985732075359, 0.484490089563701))
## Create the M matrix for different examples
disEx1 = discomb(lnpfshr, selnpfshr, t1, t2, id, reference.group = "Ofa",
                 sm="md", data = Ex1)
M1 <- disEx1$X.matrix

netgraph(disEx1, number.of.studies = F, points = T, cex.points = 3)

# Some examples using checkIdentifiable function
checkIdentifiable(M1, v = c(1,0,0,0,0,0,0,-1))
checkIdentifiable(M1, v = c(0,0,0,0,1,0,0,-1))
checkIdentifiable(M1, v = c(1,0,0,0,0,-1,0,0))

# Example: Ranking novel targeted agent components 
rank1 <- cnmaRank(disEx1, set = c("Duv",
                                  "Ibr",
                                  "Ide",
                                  "Ubl",
                                  "Ven"), verbose = F,small.values = "undesirable", re = F) # fails

rank2 <- cnmaRank(disEx1, set = c("Duv",
                                  "Ibr",
                                  "Ide",
                                  "Ubl"), verbose = F,small.values = "undesirable", re = F)

poth(rank2)

# Example: Ranking observed treatments
rank3 <- cnmaRank(disEx1, set = c("Ofa",
                                  "Ide+Ofa",
                                  "Ide+Rit",
                                  "Rit",
                                  "Ubl+Ibr",
                                  "Ven+Rit",
                                  "Ben+Rit",
                                  "Duv",
                                  "Ibr",
                                  "Ibr+Ben+Rit",
                                  "Ibr+Rit",
                                  "Ide+Ben+Rit"), 
                  verbose = T,small.values = "undesirable", re = F) # can't do

# run on only treatments from connected sub-network 1
rank4 <- cnmaRank(disEx1, set = c("Ide+Ofa",
                                  "Ide+Rit",
                                  "Ofa",
                                  "Rit",
                                  "Ubl+Ibr",
                                  "Duv",
                                  "Ibr",
                                  "Ibr+Rit"), 
                  verbose = T,small.values = "undesirable", re = F) 
poth(rank4)

rank5 <- cnmaRank(disEx1, set = c("Ibr+Ben+Rit",
                                  "Ide+Ben+Rit",
                                  "Ven+Rit",
                                  "Ben+Rit"), 
                  verbose = T,small.values = "undesirable", re = F) # can't do

# Other possible examples?? -----------------------------

data("Linde2016")
face <- subset(Linde2016, id %in% c(16, 24, 49, 118))

# Conduct random effects network meta-analysis
#
net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
                data = face, ref = "placebo", sm = "OR", common = FALSE)
net1

# Additive model for treatment components (with placebo as inactive
# treatment)
#
nc1 <- netcomb(net1)
dim(nc1$X.matrix)


# other data

stay <- read.csv("length_of_stay.csv")

staypw <- pairwise(treat = treatment,
                   n = N, studlab = ID, mean = Mean, sd = SD, data = stay)

staynw <- netmeta(staypw)

staycnma <- netcomb(staynw)

staycnma$X.matrix

checkIdentifiable(staycnma$X.matrix, v = c(0,0,0,0,0,0,0))

data(Linde2016)

net1 <- netmeta(lnOR, selnOR, treat1, treat2, id,
                data = Linde2016, ref = "placebo", sm = "OR", common = FALSE)

nc1 <- netcomb(net1)
nc1

allrank <-cnmaRank(nc1, set = nc1$trts, 
         verbose = F,small.values = "undesirable", re = F)
allrank
poth(allrank)



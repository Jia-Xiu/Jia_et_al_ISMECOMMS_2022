# Using sample-specific rarity cutoffs to define the rare and common biopsheres
# Author: Xiu Jia
# Date: 09-04-2019

# load library
library(vegan)

# load rarefied feature table
com <- read.csv(choose.files(), sep=",",  header=1, row.names=1)
com <- t(com[, 1:60]) # remove columns with taxa information
str(com)

# calculate Chao1 for estimation of the sequencing depth
# to be notice, Chao1 could only represent the lower bundary of richness estimation
Chao <- as.data.frame(t(estimateR(com)))
Chao$slope <- Chao$S.obs/Chao$S.chao1
head(Chao)


# built a empty matrix to store the sample-specific rarity cutoffs
cutoffs <- matrix(NA, nrow(com), 3)
row.names(cutoffs) <- row.names(com)
cutoffs[,2] <- Chao$slope
cutoffs[,3] <- Chao$S.obs
colnames(cutoffs) <- c("Rarity.cutoffs", "Slopes", "S.obs")

# using method calculate H-index to calculate rarity cutoff per sample
for (j in 1:nrow(com)) {
  com_j = sort(as.numeric(com[j,]), decreasing = TRUE)
  com_j <- com_j[com_j!=0]
  slope <- cutoffs[j,2]
  for (i in 1:length(com_j)){
    if (com_j[i]>=i*slope){
      H=i
    }
  }
  cutoffs[j, 1] <- H
}
cutoffs <- as.data.frame(cutoffs)
head(cutoffs)


# generate the dataset of the rare biosphere
df <- com
for (j in 1:nrow(df)) {
  for (i in 1:ncol(df)) {
    if (df[j, i] > cutoffs[j,1]) {
      df[j, i] <- NA
    }
  }
}

df[is.na(df)] <- 0
df <- df[, colSums(df)!=0] 
df <- df[, rowSums(df)!=0] 
write.csv(t(df), "rare_biosphere_specitic_cutoffs.csv")

# generate the dataset of the common biosphere
df <- com
for (j in 1:nrow(df)) {
  for (i in 1:ncol(df))
    if (df[j, i] <= cutoffs[j,1]) {
      df[j,i] <- NA
    }
}

df[is.na(df)] <- 0
df <- df[, colSums(df)!=0] 
df <- df[, rowSums(df)!=0] 
write.csv(t(df), "common_biosphere_specitic_cutoffs.csv")

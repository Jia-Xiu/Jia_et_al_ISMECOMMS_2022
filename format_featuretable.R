# Format feature table
# Author: Xiu Jia
# Date: 29-10-2018

rm(list=ls())

# load the directory
directory = '~/Dropbox'
subfolder = 'dada2'

setwd(paste(directory, subfolder, sep="/"))
getwd()

# load packages
library(vegan)

##########################################################################################################

# No need to run this section anymore!!!

##########################################################################################################
df1 <- read.csv("feature_table-silva-filtered-nontaxa.csv", row.names=1, header=1, sep=";")
df2 <- read.csv("taxonomy-silva.csv", row.names=1, header=1, sep=";")
df <- merge(df1, df2,  by="row.names")  # merge by feature id
# write.csv(de, "OTU-table-dada2-silva-raw.csv")
# remove redundant details, e.g. "D_0__", "cloroplast", "mitochondria", etc.

df <- read.csv("feature_table-dada2-silva.csv", row.name=1, header=1, sep=";")
df2 <- df[ , -which(names(df) %in% c("Kindom","Phylum","Class", "Order", "Family", "Genus", "Species"))]
str(df2)

# rarefy OTU table
df2 <- df2[rowSums(df2)!=0,] # remove 0 abundance taxa
df2 <- t(df2)
(raremax <- min(rowSums(df2)))
set.seed(10000)
rare <- rrarefy(df2, raremax)
rare <- rare[, colSums(rare)!=0]
rare <- t(rare)
dim(rare)
colSums(rare)
rare[1:5, 1:6]

# write.csv(rare,"otutable_dada2_silva_filtered_rarefied.csv")
# get taxonomy info
taxa <- read.csv("otutable-dada2-silva.csv", row.name=1, header=1, sep=";")
taxa <- taxa[, -c(1:60)]
str(taxa)

# combine taxonomy info with the rarefied table
com <- merge(rare, taxa, by="row.names")  # merge by feature id
row.names(com) <- com$Row.names
com <- com[,-1]
str(com)
#write.csv(com, "schier_cdna_otutable_dada2_silva_filtered_rarefied_taxonomy.csv")

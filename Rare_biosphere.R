# split the rare and abundant biospheres, and analysis for those two biosphere
# Author: Jia Xiu 
# Date: 17-04-2019 (built on 09-11-2018)

# load packages 
rm(list=ls())

library(VennDiagram)
library(ggplot2)
library(RColorBrewer); display.brewer.all()
library(ggpubr)
library(vegan)
library(ape) 
library(reshape2)
library(doBy) 
library(plyr)
library(GUniFrac)


# load directory 
setwd()


# plot theme
mytheme <- theme_bw()+
  theme(text = element_text(size=12),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.box.background = element_rect(),
        legend.box.margin = margin(1,1,1,1),
        legend.position = "right", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(1, 0.5, 1, 0.5),"cm")) 

# load the rarefied feature table ----------------------------------------------------------------------------
feature_table <- read.csv("schier_cdna_feature_table_dada2_silva_rarefied_taxonomy.csv", sep=",", 
                          header=1, row.names=1)
feature_table$Phylum <- gsub("Candidate division ", "", feature_table$Phylum)
feature_table[feature_table==""] <- NA
levels(factor(feature_table$Phylum))
wholeDS <- feature_table[, c(1:60)]
wholeDS <- t(wholeDS)
str(wholeDS)
cat("\nrarefied to", rarefactiondepth <- mean(rowSums(wholeDS)))
cat("\nthe number of samples is:", nrow(wholeDS), "\nthe number of species/ASVs is:", ncol(wholeDS))
cat("\nthe range of sequence number among samples is:", range(rowSums(wholeDS)))


# Calculate the lowest relative abundance (%) of the maximum ESV occurance in all samples
colMax <- function(data) sapply(data, max, na.rm = TRUE)
cat("rarity cutoff no more than", round((min(colMax(feature_table[,c(1:60)]))/median(rowSums(wholeDS)))*100,2), "%\n")

# Set the cutoff for rarity
cutoff = 0.1/100

# source the trucate function
source("TruncateTable.r") #https://github.com/Jia-Xiu/rare_biosphere_assembly_2020/blob/master/TruncateTable.R

# The truncated datasets can be stored as follows: 
truncated_ds_dominant <-TruncateTable(wholeDS, cutoff, typem="dominant") 
str(truncated_ds_dominant)
#write.csv(t(truncated_ds_dominant), paste("truncated_ds_dominant", cutoff, "cutoff.csv", sep="_"))

truncated_ds_rare_without_dominant <-TruncateTable(wholeDS, cutoff, typem="rare") 
str(truncated_ds_rare_without_dominant)
#write.csv(t(truncated_ds_rare_without_dominant),  paste("truncated_ds_rare_without_dominant", cutoff, "cutoff.csv", sep="_"))

name.to.keep <- row.names(t(truncated_ds_rare_without_dominant))
truncated_ds_rare_taxonomy <- subset(feature_table, row.names(feature_table) %in% name.to.keep)
#write.csv(truncated_ds_rare_taxonomy, paste("truncated_ds_rare_taxonomy", cutoff, "cutoff.csv", sep="_"))
truncated_ds_rare <- truncated_ds_rare_taxonomy[,c(1:60)]
truncated_ds_rare <- t(truncated_ds_rare)


# define different types of rarity and commonness----------------------------------------------------------------------------------------------
df <- as.data.frame(t(truncated_ds_rare))

warning("Check whether rarefactiondepth (", rarefactiondepth, ") is the value your mean (i.e. 31500)?")

# source rarity type function
source("rarity_type.R")

truncated_ds_rare_classified <- rarity_type(df, cutoff, rarefactiondepth)
dim(truncated_ds_rare_classified)
truncated_ds_rare_classified[1:5, 59:61]
cat("Is the number of categoried species (", sum(table(truncated_ds_rare_classified$vector_category)), 
    ")\nsame as the number of rare species (", ncol(truncated_ds_rare_without_dominant), ") ?\n",
    sum(table(truncated_ds_rare_classified$vector_category)) == ncol(truncated_ds_rare_without_dominant))

# seperate different types of rarity
# permanently rare
A <- subset(truncated_ds_rare_classified, vector_category == "A")
# transiently rare
B <- truncated_ds_rare_classified[truncated_ds_rare_classified$vector_category == "B", ]
# conditionally rare/dominant
C <- truncated_ds_rare_classified[truncated_ds_rare_classified$vector_category == "C", ]
# permanently dominant
D <- as.data.frame(subset(t(truncated_ds_dominant), 
                          !(row.names(t(truncated_ds_dominant)) %in% row.names(t(truncated_ds_rare_without_dominant)))))
D$vector_category <- "D"
str(D)

cat("Is the sum of ASVs in each rarity type (", nrow(A)+nrow(B)+nrow(C)+nrow(D), 
    ")\nsame as the total number of ASVs in the whole data set (",  ncol(wholeDS), ") ?\n",
    nrow(A)+nrow(B)+nrow(C)+nrow(D) == ncol(wholeDS))
#write.csv(A[, 1:60],  paste("truncated_ds_type_A", cutoff, "cutoff.csv", sep="_"))





# plot types of rarity and commonness -----------------------------------------------

# Summerizing each type of rarity 
source("SummarizeRarityTypes.R")

df <- SummarizeRarityTypes(truncated_ds_rare_without_dominant, truncated_ds_rare_classified, rarefactiondepth, binary)
head(df)
colMeans(df);colSds(df)/sqrt(nrow(df))

group_info <- data.frame(row.names=rownames(df),  t(as.data.frame(strsplit(as.character(row.names(df)), "_"))))

df1 <- data.frame(Year=as.factor(group_info[,2]), Month=as.factor(group_info[,3]), df)
df1$Month <- factor(df1$Month, levels=c("5", "7", "9", "11"), labels=c("M", "J", "S", "N"))
df1$Year <- factor(df1$Year, levels=c("0", "10", "40", "70", "110"))

df2 <- melt(df1, id=c("Year","Month"))
df2$variable <- factor(df2$variable, levels=c("PermanentlyRare", "TransientlyRare", "ConditionallyRare"),
                       labels=c("Permanently Rare", "Transiently Rare", "Conditionally Rare"))

dstats<-function(x)(c(n=length(x), mean=mean(x), sd=sd(x), se=sd(x)/sqrt(length(x))))
data<-summaryBy(value~Year+Month+variable, data=df2, FUN=dstats)
head(data)

# custom colors
# my_palette = c(brewer.pal(8, "Pastel1")[c(1:3)])
my_palette = c(brewer.pal(8, "Greys")[c(1:3)])

# stacked-bar plot
(f2 <- ggplot(data, aes(x=Month, y=value.mean, fill=variable)) + 
    geom_bar(stat="identity", width=0.8, colour="black") +
    scale_fill_manual(values=my_palette) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) + #2500
    facet_grid(~Year, switch = "x", scales = "free_x") +
    guides(fill=guide_legend(title="Types of rarity"))+
    xlab("Stage of succession (Years)") +
    ylab("Relative abundance (%)") +
    #ylab("Number of ASVs") +
    ggtitle("Rare biosphere") +
    mytheme)

(f <- ggarrange(f1, f2, labels = c("A", "B"), common.legend = TRUE, legend = "right", ncol = 2))

ggsave(paste("Rarity_types_sn_ab", cutoff, "greys.png", sep = "_"), 
       width = 17, height = 8, units = "cm", f, scale = 2, dpi = 300)


# summarize dominant types ---------------------------------
source("SummarizeCommonnessTypes.R")
df <- SummarizeCommonnessTypes(truncated_ds_dominant, truncated_ds_rare_without_dominant, rarefactiondepth, binary)
head(df)
colMeans(df);colSds(df)/sqrt(nrow(df))

group_info <- data.frame(row.names=rownames(df), t(as.data.frame(strsplit(as.character(row.names(df)), "_"))))

df1 <- data.frame(Year=as.factor(group_info[,2]), Month=as.factor(group_info[,3]), df)
df1$Month <- factor(df1$Month, levels=c("5", "7", "9", "11"), labels=c("M", "J", "S", "N"))
df1$Year <- factor(df1$Year, levels=c("0", "10", "40", "70", "110"))

df2 <- melt(df1, id=c("Year","Month"))
df2$variable <- factor(df2$variable, levels=c("ConditionallyDominant", "PermanentlyDominant"),
                       labels=c("Conditionally Common", "Permanently Common"))

data<-summaryBy(value~Year+Month+variable, data=df2, FUN=dstats)
head(data)

#my_palette = c(brewer.pal(8, "Pastel2")[c(5:6)])
my_palette = c(brewer.pal(8, "Greys")[c(4, 6)])
# stacked-bar plot
(p2 <- ggplot(data, aes(x=Month, y=value.mean, fill=variable)) + 
    geom_bar(stat="identity", width=0.8, colour="black") +
    scale_fill_manual(values=my_palette) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 2500)) + #2500
    facet_grid(~Year, switch = "x", scales = "free_x") +
    guides(fill=guide_legend(title="Commonness/Rarity types"))+
    xlab("Stage of succession (Years)") +
    #ylab("Relative abundance (%)") +
    ylab("Number of ASVs") +
    ggtitle("Common biosphere") +
    mytheme)

(f4 <- ggplot(data, aes(x=Month, y=value.mean, fill=variable)) + 
    geom_bar(stat="identity", width=0.8, colour="black") +
    scale_fill_manual(values=my_palette) + # 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 2600)) + #1800
    facet_grid(~Year, switch = "x", scales = "free_x") +
    guides(fill=guide_legend(title="Rarity/Commonness types"))+
    xlab("Stage of succession (Years)") +
    #ylab("Relative abundance (%)") +
    ylab("Number of ASVs") +
    mytheme)

# combine plots
(fp <- ggarrange(p1, f2, p2, f1, labels = c("A", "B", "C", "D"), common.legend = TRUE, legend = "right", ncol = 2, nrow = 2))

ggsave(paste("Rarity_Common_types_sn_ab", cutoff, "greys.pdf", sep = "_"), 
       width = 13, height = 11, units = "cm", fp, scale = 2)




# Beta diversity analysis -------------------------------------------------------------------------------------
# getPalette = colorRampPalette(brewer.pal(8,"Set2")), getPalette(8)
mytheme <- theme_bw()+
  theme(text = element_text(size=15),
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

# combine the rare and dominant biosphere together
row.names(truncated_ds_dominant) <- gsub("cDNA", "dominant", row.names(truncated_ds_dominant))
dominant <- as.data.frame(t(truncated_ds_dominant))
dominant[1:5, 1:3] 

row.names(truncated_ds_rare_without_dominant) <- gsub("cDNA", "rare", row.names(truncated_ds_rare_without_dominant))
rare <- as.data.frame(t(truncated_ds_rare_without_dominant))
rare[1:5, 1:3] 

com <- transform(merge(dominant, rare,
                       by="row.names", all=TRUE), row.names=Row.names, Row.names=NULL)  
com[is.na(com)] <- 0
com <- t(com)
com[1:5, 1:2]
dim(com)

# Unifrac
tree <- read.tree("wholeDS_tree_matched_to_com_table.tre")
str(tree)
str(com)
unifracs <- GUniFrac(com, tree, alpha=c(0, 0.5, 1))$unifracs
str(unifracs)

dw <- unifracs[, , "d_1"] # Weighted UniFrac
dw[56:65, c(1,61)]
du <- unifracs[, , "d_UW"] # Unweighted UniFrac
dv <- unifracs[, , "d_VAW"] # Variance adjusted weighted UniFrac
d0 <- unifracs[, , "d_0"] # GUniFrac with alpha 0
d5 <- unifracs[, , "d_0.5"] # GUniFrac with alpha 0.5

dist <- as.dist(as.matrix(dw))

# Bray Curtis
dist <- vegdist(com, method="bray", binary=FALSE, diag=1) 
str(dist)

# Some distance measures may result in negative eigenvalues. In that case, add a correction:
re <- pcoa(dist, correction="none", rn=NULL) # ?correction = "cailliez")
str(re)

group_info <- data.frame(row.names=row.names(re$vectors), 
                         t(as.data.frame(strsplit(as.character(row.names(re$vectors)), "_"))))
head(group_info)

df <- data.frame(x=re$vectors[,1],y=re$vectors[,2],
                 Biosphere=as.factor(group_info[,1]),
                 Year=as.factor(group_info[,2]),
                 Month=as.factor(group_info[,3]),
                 replicates=as.factor(group_info[,4]))

df$Year <- factor(df$Year, levels=c("0", "10", "40", "70", "110"))

df$Month <- factor(df$Month, levels=c("5", "7", "9", "11"), 
                   labels=c("May", "July", "September", "November"))

df$Biosphere <- factor(df$Biosphere, levels=c("rare", "dominant"), labels=c("Rare ", "Common")) # , #, "cDNA" "Whole"

str(df)

(f1 <- ggplot(df, aes(x, y, shape=Year, fill=Biosphere))+
    geom_point(size=6, alpha=0.7)+ 
    labs(x=paste("PCoA1 (", round(re$values$Rel_corr_eig[1]*100, 2), "%)", sep=""), 
         y=paste("PCoA2 (", round(re$values$Rel_corr_eig[2]*100, 2), "%)", sep=""),
         title = "")+
    scale_fill_manual(values = c("white", "gray"), guide=guide_legend(override.aes = list(shape=21)))+
    scale_shape_manual(values=c(23, 22, 21, 24, 25))+ 
    mytheme
)

# combine plots
f <- ggarrange(f1, f2, f3, labels = c("A 0.2%", "B  0.1%", "C  0.05%"), 
               common.legend = TRUE, legend = "right", ncol = 3, nrow = 1)
f

ggsave("PCoA_bray_rare_dominant_three_cutoffs.pdf", 
       width = 17, height = 5, units = "cm", f, scale = 2)
ggsave("PCoA_bray_rare_dominant_three_cutoffs.png",
       width = 17, height = 5, units = "cm", f, scale = 2, dpi = 300)

ggsave(paste("PCoA_weighted_unifrac_rare_dominant", cutoff, "cutoff.pdf", sep="_"), 
       width = 7.5, height = 6, units = "cm", f1, scale = 2)
ggsave(paste("PCoA_weighted_unifrac_rare_dominant", cutoff, "cutoff.png", sep="_"),
       width = 7.5, height = 6, units = "cm", f1, scale = 2, dpi = 300)





# PERMANOVA ------------------------------------------------------------------------------------------
com[1:6,1:5]

df <- data.frame(row.names=rownames(com),t(as.data.frame(strsplit(rownames(com),"_"))))

df <- rename(df, c("X1"="Dataset", "X2"="Year", "X3"="Month", "X4"="replicates"))
colnames(df) <- c("Dataset", "Year", "Month", "Replicates")
df$Year <- factor(df$Year, levels=c("0", "10", "40", "70", "110"))

df$Month <- factor(df$Month, levels=c("5", "7", "9", "11"), 
                   labels=c("May", "July", "September", "November"))
head(df)

# two way permanova (successional stages vs. season)
result <- adonis(com ~ Dataset*Year, data=df, method="bray", permutation=9999) # two way permanova
result
write.csv(result$aov.tab, paste("PERMANOVA", name, cutoff, "result1.csv", sep="_"))



# Phyla composition ------------------------------------------------------------------------------------
# Rare & dominant species composition in Phyla level (relative abundance) 
# add taxonomy info
row.names(truncated_ds_dominant) <- gsub("cDNA", "dominant", row.names(truncated_ds_dominant))
dominant <- as.data.frame(t(truncated_ds_dominant))
dominant[1:5, 1:3] 

row.names(truncated_ds_rare_without_dominant) <- gsub("cDNA", "rare", row.names(truncated_ds_rare_without_dominant))
rare <- as.data.frame(t(truncated_ds_rare_without_dominant))
rare[1:5, 1:3] 

com <- transform(merge(dominant, rare,
                       by="row.names", all=TRUE), row.names=Row.names, Row.names=NULL)  

com <- transform(merge(feature_table[, -c(1:61)], com, by="row.names"), row.names=Row.names, Row.names=NULL)  

com$Phylum <- gsub("Candidate division ", "", com$Phylum)

com[is.na(com)] <- 0
com <- com[, -c(1, 3:8)]
com <- com[!is.na(com$Phylum), ]
com$Phylum
com[1:45, 1:2]

phyla <- aggregate(com[,-1], list(com$Phylum), sum)
row.names(phyla) <- phyla$Group.1
phyla <- phyla[, -1]
phyla <- t(phyla)
phyla[1:5, 1:3]
str(phyla)

# relative abundance
df <- 100*phyla/rarefactiondepth
df[is.na(df)] <- 0
df <- t(df)
df <- as.data.frame(df)
df$mean <- rowMeans(df)
df <- df[with(df, order(-mean)), ] #for stacked-bar plot order(mean), for line plot order(-mean)
df$mean <- NULL
df <- t(df)
colnames(df)
str(df)
df[1:5, 1:3]

group_info <- data.frame(row.names=rownames(df), 
                         t(as.data.frame(strsplit(as.character(row.names(df)), "_"))))

df1 <- data.frame(Year=as.factor(group_info[,2]),
                  Month=as.factor(group_info[,3]),
                  Group=as.factor(group_info[,1]),
                  df)
df1$Month <- factor(df1$Month, levels=c("5", "7", "9", "11"), 
                    labels=c("M", "J", "S", "N"))

df1$Year <- factor(df1$Year, levels=c("0", "10", "40", "70", "110"))

df1$Group <- factor(df1$Group, levels=c("rare", "dominant"), 
                    labels=c("Rare", "Common"))

df2 <- melt(df1, id=c("Year","Month","Group"))
names(table(df2$variable))
str(df2)

dstats<-function(x)(c(n=length(x), mean=mean(x), sd=sd(x),  se=sd(x)/sqrt(length(x))))
data<-summaryBy(value~Year+Month+Group+variable, data=df2, FUN=dstats)
head(data)

# stacked-bar plot
f <- ggplot(data, aes(x=interaction(Month,Year), y=value.mean, fill=Group)) + 
  geom_bar(stat="identity", width=0.7) +
  scale_fill_manual(values = c("#66c2a5", "#fc8d62"))+ 
  facet_wrap(~variable, ncol=5, scales="free_y") +
  guides(fill = guide_legend(title = NULL))+
  labs(x=" ",y="Relative abundance (%)", title=" ") +
  theme_classic()+
  theme(axis.text.x=element_text(size=rel(0.7), angle=45, hjust=1, vjust=1),
        axis.title.y = element_text(size = rel(1.2), angle = 90),
        legend.position = c(0.92, 0.04), # c(0,0) bottom left, c(1,1) top-right.
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.title = element_text(face = "bold"))
f

pdf(paste("Phyla_rare_vs_dominate", cutoff, "cutoff.pdf", sep="_"), width=9, height=14)
print(f)
dev.off()


# taxonomy distribution of all types of rarity and commonness -----------------------------------------------

types <- c('Permanently_Rare', 'Transiently_Rare', 'Conditionally_Rare', 
           'Conditionally_Common', 'Permanently_Common')

datalist = list()

for (type in types) {
  if (type == 'Permanently_Rare') {
    name.to.keep <- row.names(A)
    df <- subset(t(truncated_ds_rare_without_dominant), 
                 row.names(t(truncated_ds_rare_without_dominant)) %in% name.to.keep)
  } else if (type  == 'Transiently_Rare') {
    name.to.keep <- row.names(B)
    df <- subset(t(truncated_ds_rare_without_dominant), 
                 row.names(t(truncated_ds_rare_without_dominant)) %in% name.to.keep)
  } else if (type  == 'Conditionally_Rare') {
    name.to.keep <- row.names(C)
    df <- subset(t(truncated_ds_rare_without_dominant), 
                 row.names(t(truncated_ds_rare_without_dominant)) %in% name.to.keep)
  } else if (type  == 'Conditionally_Common') {
    name.to.keep <- row.names(C)
    df <- subset(t(truncated_ds_dominant), 
                 row.names(t(truncated_ds_dominant)) %in% name.to.keep)
  } else if (type  == 'Permanently_Common') {
    name.to.keep <- row.names(D)
    df <- subset(t(truncated_ds_dominant), 
                 row.names(t(truncated_ds_dominant)) %in% name.to.keep)
  }
  
  df <- transform(merge(feature_table[63:67], df, by = "row.names"), row.names=Row.names, Row.names=NULL)
  df <- df[, -c(2:5)]
  
  df <- aggregate(df[,-1], list(df$Phylum), sum)
  row.names(df) <- df$Group.1; df <- df[, -1]
  df <- t(df)
  
  # relative abundance
  df <- 100*df/rarefactiondepth
  df[is.na(df)] <- 0
  df <- as.data.frame(t(df))
  df <- df[with(df, order(df$cDNA_110_11_A)), ]
  df <- t(df)
  colnames(df)
  str(df)
  df[1:5, 1:3]
  
  group_info <- data.frame(row.names=rownames(df), 
                           t(as.data.frame(strsplit(as.character(row.names(df)), "_"))))
  
  df1 <- data.frame(Year=as.factor(group_info[,2]),
                    Month=as.factor(group_info[,3]),
                    Type=type,
                    df)
  
  df1$Month <- factor(df1$Month, levels=c("5", "7", "9", "11"), 
                      labels=c("M", "J", "S", "N"))
  
  df1$Year <- factor(df1$Year, levels=c("0", "10", "40", "70", "110"))
  
  df1 <- melt(df1, id=c("Year","Month", "Type"))
  names(table(df1$variable))
  str(df1)
  head(df1)
  datalist[[type]] <- df1
}

str(datalist)
big_data = do.call(rbind, datalist)
big_data$Type <- factor(big_data$Type, levels = c('Permanently_Rare', 'Transiently_Rare', 'Conditionally_Rare', 
                                                  'Conditionally_Common', 'Permanently_Common'),
                        labels = c('Permanently Rare', 'Transiently Rare', 'Conditionally Rare', 
                                   'Conditionally Common', 'Permanently Common'))
head(big_data)
str(big_data)


dstats<-function(x)(c(n=length(x), mean=mean(x), sd=sd(x),  se=sd(x)/sqrt(length(x))))
data<-summaryBy(value~Year+Month+Type+variable, data=big_data, FUN=dstats)
data <- data[with(data, order(data$value.mean)), ]
head(data)

(colourCount = length(unique(data$variable)))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col=sample(col_vector, colourCount)
pie(rep(1,colourCount), col=sample(col_vector, colourCount))

# stacked-bar plot
f <- ggplot(data, aes(x= interaction(Month,Year), y=value.mean, fill=variable)) + 
  geom_bar(stat="identity", width=0.8) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 75))+
  scale_fill_manual(values=col)+
  facet_wrap(~Type, ncol=5, scales="free_y") +
  labs(x=" ",y="Relative abundance (%)", title=" ") +
  guides(fill=guide_legend(title="Phyla"))+
  theme_bw()+
  theme(text = element_text(size=15),
        axis.text.x=element_text(size=9, angle=45, hjust=1, vjust=1),
        legend.text=element_text(size=11),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.box.background = element_rect(),
        legend.box.margin = margin(1,1,1,1),
        legend.position = "bottom", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(1, 0.5, 1, 0.5),"cm")) 
f

pdf(paste("Phyla_Types", cutoff, "cutoff.pdf", sep="_"), width=14, height=9)
print(f)
dev.off()


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
library(ape) # for pcoa 
library(reshape2)
library(doBy) # se function
library(plyr) # rbind.fill # apply function
library(GUniFrac) # UniFrac
#library(rmarkdown)

# load directory 
directory = '~/Dropbox/'
subfolder = 'Schier/RNA'

setwd(paste(directory, subfolder, sep="/"))
getwd()

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


# Correlation between occurance and abundance for each typ of rarity/commonness --------------------
# permanently rare - A 
# transiently rare - B
# conditionally rare/dominant - C 
# permanently dominant - D
Ac <- data.frame(ASVs = row.names(A), catergory = A[,61], row.names = row.names(A))
Bc <- data.frame(ASVs = row.names(B), catergory = B[,61], row.names = row.names(B))
Cc <- data.frame(ASVs = row.names(C), catergory = C[,61], row.names = row.names(C))
Dc <- data.frame(ASVs = row.names(D), catergory = D[,61], row.names = row.names(D))
com_cat <- rbind(Ac, Bc, Cc, Dc)

# calculate occurance of each ASVs across all samples
occurance <- data.frame(occurance = rowSums(t(wholeDS) != 0), 
                        t(wholeDS), 
                        row.names = row.names(t(wholeDS)))

warning("Becareful of put species as rows! Put it in a wrong way, my crappy computer almost dead!!!")
com_cat <- transform(merge(com_cat, occurance, by="row.names", all=TRUE), row.names=Row.names, Row.names=NULL)  
com_cat[1:5, 1:4]

com_cat <- melt(com_cat, id=c("ASVs", "catergory", "occurance"))
com_cat <- com_cat[com_cat$value != 0,]
com_cat$variable <- as.character(com_cat$variable)
str(com_cat)

group_info <- data.frame(t(as.data.frame(strsplit(com_cat$variable,"_"))))

df <- data.frame(com_cat,
                 Year=as.factor(group_info[,2]),
                 Month=as.factor(group_info[,3]),
                 replicates=as.factor(group_info[,4]))
head(df)

# loess fit plot ----------------------------------------------------------------------------
(my_palette = c(brewer.pal(9, "Set1")[c(1:4)]))
(p <- ggplot(df, aes(x = occurance, y = value, color = catergory, shape = catergory)) + 
    geom_point(size = 2, alpha = 0.9)+
    facet_grid(Year ~ Month) +
    geom_smooth(aes(fill = catergory), method = loess,  linetype="dashed") +
    scale_color_manual(values = my_palette) + 
    scale_size_manual(values = c(16, 15, 17, 18)) +
    labs(x="Occurance", y="Abundance (reads)")+
    mytheme)
ggsave("Occurance_abundance.png", width = 17, height = 15, units = "cm", p, scale = 1.5, dpi = 300)


# CCA for conditionally and permanently rare -------------------------------------------------------
# inertia explained by axes in CCA plots
#' @title Get percent of total inertia explained by principal or constrained axes
#' @param mod cca.object
#' @param axes A numeric vector indicating the axes for which percent explained inertia should be returned for
#' @return Returns a named vector containing explained inertia as percentages. 
#' Returns constrained inertia fo CCA and unconstrained for CA
axis.expl <- function(mod, axes = 1:2) {
  
  if(is.null(mod$CCA)) {
    sapply(axes, function(i) {
      100*mod$CA$eig[i]/mod$tot.chi
    })
  } else {
    sapply(axes, function(i) {
      100*mod$CCA$eig[i]/mod$tot.chi
    })
  }
  
}


# Environmental factors
env <- read.csv("SoilParameters_2017.csv", sep=",", row.names=1, header=TRUE)
row.names(env) <- paste("cDNA", row.names(env), sep = "_")
env <- scale(env)
env <- as.data.frame(env)
head(env)

subcom <- t(D[, -61]) # A or C
# change to present/absent data
subcom[subcom > 0] <- 1
subcom <- subcom[rownames(env),]
str(subcom)

# Use adonis to find significant environmental variables
com.adonis <- adonis(subcom ~ ., data=env)
com.adonis$aov.tab
#write.csv(com.adonis$aov.tab, "adonis_permanently_common.csv")

#Extract the best variables, remove NA entries
bestEnvVariables <- rownames(com.adonis$aov.tab)[com.adonis$aov.tab$"Pr(>F)" <= 0.01]
bestEnvVariables <- bestEnvVariables[!is.na(bestEnvVariables)]
bestEnvVariables

# The length of first DCA axis > 4 S.D. indicates heterogeneous dataset on which unimodal methods should be used
# while the length < 3 S.D. indicates homogeneous dataset for which linear methods are suitable 
# when the length between 3 and 4 S.D., both linear and unimodal methods are OK. 
dca <- decorana(subcom)
dca

# We are now going to use only those environmental variables in cca that were found significant
# eval - Evaluate an R expression in a specified environment.
im.cca <- eval(parse(text=paste("sol <- cca(subcom ~ ", do.call(paste, c(as.list(bestEnvVariables), sep=" + ")),", data=env)", sep="")))
sum.cca <- summary(im.cca)
# tests the significance of the variation in species composition explained by explanatory variables
anova.cca(im.cca)
coef(im.cca)
plot(im.cca)

#generalized variance-inflation factors
(vif<-vif.cca(im.cca))

# rda-R2
(R2 <- RsquareAdj(im.cca)$r.squared)
(R2adj <- RsquareAdj(im.cca)$adj.r.squared)

#site scores
scrs <- scores(sol, display=c("wa","lc","bp","cn"))

#Check the attributes
attributes(scrs)

#Extract site data first
df <- data.frame(scrs$sites, t(as.data.frame(strsplit(rownames(scrs$sites),"_"))))
colnames(df) <- c("CCA1","CCA2", "Datasets", "Year", "Month", "replicates")

df$Year <- factor(df$Year, levels=c("0", "10", "40", "70", "110"))
df$Month <- factor(df$Month, levels=c("5", "7", "9", "11"), 
                   labels=c("May", "July", "September", "November"))
head(df)

labs <- axis.expl(im.cca )

#Draw biplots
multiplier <- ordiArrowMul(scrs$biplot) 
#Support functions to assist with drawing of vectors (arrows) on ordination plots. 
#ordiArrowMul finds the multiplier for the coordinates of the head of the vector such that they accupy fill proportion of the plot region. 
#ordiArrowTextXY finds coordinates for the locations of labels to be drawn just beyond the head of the vector.

df_arrows <- scrs$biplot*multiplier
colnames(df_arrows) <- c("x","y")
df_arrows <- as.data.frame(df_arrows)
#rownames(df_arrows) <- factor(rownames(df_arrows), labels = c("SWC", "Nitrate", "SOM", "Sodium", "pH"))
#labels = c("Soil water content", "Nitrate/Nitrate", "Soil organic carbon", "Sodium", "Total Nitrogen", "pH"))
df_arrows

# scatter plot
(p3 <- ggplot() + geom_point(data=df, aes(CCA1, CCA2, fill=Year, shape=Month), size=7, alpha=0.8)+
    scale_fill_brewer(palette="Set3", guide=guide_legend(override.aes = list(shape=21)))+
    scale_shape_manual(values=c(24, 22, 23, 21))+
    geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y, colour=Year),
                 arrow = arrow(length = unit(0.3, "cm")), color="#990000", alpha=0.8)+
    geom_text(data=as.data.frame(df_arrows*1.1), aes(x, y, label = rownames(df_arrows)), color="#990000", alpha=0.8)+
    labs(x=paste0(names(labs[1]), " (", sprintf("%.1f", labs[1]), "%)"), 
         y=paste0(names(labs[2]), " (", sprintf("%.1f", labs[2]), "%)"), 
         title="Permanently Common") + #Permanently Conditionally 
    annotate('text', label = paste("R2adj =", round(R2adj, 3)), 
             x = Inf, y = -Inf, size = 4, hjust = 1.2, vjust = -.9, size = 4) +
    mytheme)

# combine plots
p <- ggarrange(p3, p1, p2, labels = c("A", "B", "C"), common.legend = TRUE, legend = "right", ncol = 1, nrow = 3)
p

ggsave("CCA_weighted_permanently_common.pdf", width = 9.5, height = 8, units = "cm", p3, scale = 1.5)
ggsave("CCA_weighted_permenantly_conditionally_rare.png", width = 17, height = 7, units = "cm", p, scale = 1.5, dpi = 300)


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




# NMDS ----------------------------------------------------------------------------------------------------

# First step is to calculate a distance matrix. See PCOA for more information about the distance measures
# Here we use bray-curtis distance, which is recommended for abundance data
dist <- vegdist(com,  method = "bray")

# In this part, we define a function NMDS.scree() that automatically 
# performs a NMDS for 1-10 dimensions and plots the nr of dimensions vs the stress
NMDS.scree <- function(x) { #where x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

# Use the function that we just defined to choose the optimal nr of dimensions
NMDS.scree(dist)

# Because the final result depends on the initial 
# random placement of the points 
# we`ll set a seed to make the results reproducible
set.seed(123)

# Here, we perform the final analysis and check the result
# If you don`t provide a dissimilarity matrix, metaMDS automatically applies Bray-Curtis. 
NMDS1 <- metaMDS(dist, k = 2, trymax = 100, trace = F, autotransform = F)
# Do you know what the trymax = 100 and trace = F means?
# Let's check the results
NMDS1

# check the results of NMDS1 with a stressplot
stressplot(NMDS1)

# There is a good non-metric fit between observed dissimilarities (in our distance matrix) and 
# the distances in ordination space. Also the stress of our final result was ok 
# (do you know how much the stress is?). So we can go further and plot the results:
plot(NMDS1, type = "t")

str(NMDS1)

group_info <- data.frame(row.names=row.names(NMDS1$points), 
                         t(as.data.frame(strsplit(as.character(row.names(NMDS1$points)), "_"))))

df <- data.frame(x = NMDS1$points[,1], y = NMDS1$points[,2],
                 Biosphere = as.factor(group_info[,1]),
                 Year = as.factor(group_info[,2]),
                 Month = as.factor(group_info[,3]),
                 replicates = as.factor(group_info[,4]))

df$Year <- factor(df$Year, levels=c("0", "10", "40", "70", "110"))

df$Month <- factor(df$Month, levels=c("5", "7", "9", "11"), 
                   labels=c("May", "July", "September", "November"))
df$Biosphere <- factor(df$Biosphere, levels=c("rare", "dominant"), labels=c("Rare ", "Common")) # , #, "cDNA" "Whole"

str(df)


library(wesanderson)
scale_color_manual(values=wes_palette(n=5, name="Darjeeling"), guide=guide_legend(override.aes = list(shape=21)))

(f <- ggplot(df, aes(x, y, shape=Year, fill=Biosphere))+
    geom_point(size=6, alpha=0.7)+ 
    labs(x= "NMDS1", y = "NMDS2", title = "")+
    scale_fill_manual(values = c("white", "gray"), guide=guide_legend(override.aes = list(shape=21)))+
    scale_shape_manual(values=c(23, 22, 21, 24, 25))+ 
    #scale_fill_grey(start=.9, end=0, guide=guide_legend(override.aes = list(shape=21))) +
    #scale_fill_brewer(palette="Set3", guide=guide_legend(override.aes = list(shape=21)))+
    #scale_shape_manual(values=c(21, 23))+ 
    annotate('text', label = paste("Stress =", round(NMDS1$stress, 3)), 
             x = Inf, y = -Inf, size = 4, hjust = 1.2, vjust = -.8, size = 4) +
    mytheme
)

ggsave(paste("NMDS_bray_rare_dominant", cutoff, "cutoff.png", sep="_"),
       width = 7.5, height = 6, units = "cm", f, scale = 2, dpi = 300)

#PERMANOVA ------------------------------------------------------------------------------------------------
str(com)
df <- data.frame(row.names=rownames(com), t(as.data.frame(strsplit(rownames(com),"_"))))

names(df) <- c("Dataset", "Year", "Month", "replicates")
df$Year <- factor(df$Year, levels=c("0", "10", "40", "70", "110"))
df$Month <- factor(df$Month, levels=c("5", "7", "9", "11"), 
                   labels=c("May", "July", "September", "November"))
head(df)

# two way permanova (successional stages vs. season)
set.seed(123)
result <- adonis(com ~ Dataset*Month*Year, data=df, method="bray", permutation=9999) # two way permanova
result
write.csv(result$aov.tab, paste("PERMANOVA_result", cutoff, "rare_dominant.csv", sep = "_"))


# draw PCoA for rare and dominant together -------------------------------------------------------------------
com <- t(rare)
#com <- t(dominant)
com[1:5, 1:2]; dim(com)

group_info <- data.frame(row.names=rownames(com),t(as.data.frame(strsplit(rownames(com),"_"))))
head(group_info)

#PCoA analysis weighted bray curtis
bray <- vegdist(com, method="bray", binary=FALSE, diag=1) 
re <- pcoa(bray, correction="none", rn=NULL)
str(re)

df <- data.frame(x=re$vectors[,1],y=re$vectors[,2],
                 Year=as.factor(group_info[,2]),
                 Month=as.factor(group_info[,3]),
                 replicates=as.factor(group_info[,4]))

df$Year <- factor(df$Year, levels=c("0", "10", "40", "70", "110"))

df$Month <- factor(df$Month, levels=c("5", "7", "9", "11"), 
                   labels=c("May", "July", "September", "November"))
head(df)

f1 <- ggplot(df, aes(x, y, shape=Month, fill=Year))+
  geom_point(size=6, alpha=0.7)+ 
  labs(x=paste("PCoA1 (", round(re$values$Relative_eig[1]*100, 2), "%)", sep=""), 
       y=paste("PCoA2 (", round(re$values$Relative_eig[2]*100, 2), "%)", sep=""), title = " ")+
  scale_x_continuous(limits = c(-0.6, .4)) +
  scale_y_continuous(limits = c(-0.5, .4)) +
  scale_shape_manual(values=c(24, 22, 23, 21))+
  scale_fill_brewer(palette="Set3", guide=guide_legend(override.aes = list(shape=21)))+
  mytheme
f1

# combine plots
fp <- ggarrange(f1, f2, labels = c("Rare biosphere", "Common biosphere"), 
                common.legend = TRUE, legend = "right", ncol = 2, nrow = 1)
fp

ppi<-300
png(paste("PCoA_rare_dominant_biosphere", cutoff, "cutoff.png", sep="_"), width=11*ppi, height=5*ppi, res=ppi)
print(fp)
dev.off()

pdf(paste("PCoA_rare_dominant_biosphere", cutoff, "cutoff.pdf", sep="_"), width=11*ppi, height=5*ppi)
print(fp)
dev.off()

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




# heatmap ---------------------------------------------------------------------------------------------------------
df <- D[, 1:60]
colnames(df) <- sub("cDNA_", "", colnames(df))
df <- as.matrix(df)
str(df)

# conditionally rare heatmap
#df <- df[1:30, ]
df.m <- melt(df, id.vars="Names")
df.m$rescale <- log(round(df.m$value/315, 2))
df.m$Var2 <- factor(df.m$Var2, levels = c("0_5_A", "0_5_B", "0_5_C", "0_7_A", "0_7_B", "0_7_C", 
                                          "0_9_A", "0_9_B", "0_9_C", "0_11_A", "0_11_B", "0_11_C",
                                          "10_5_A", "10_5_B", "10_5_C", "10_7_A", "10_7_B", "10_7_C", 
                                          "10_9_A", "10_9_B", "10_9_C", "10_11_A", "10_11_B", "10_11_C",
                                          "40_5_A", "40_5_B", "40_5_C", "40_7_A", "40_7_B", "40_7_C", 
                                          "40_9_A", "40_9_B", "40_9_C", "40_11_A", "40_11_B", "40_11_C", 
                                          "70_5_A", "70_5_B", "70_5_C", "70_7_A", "70_7_B", "70_7_C", 
                                          "70_9_A", "70_9_B", "70_9_C", "70_11_A", "70_11_B", "70_11_C",
                                          "110_5_A", "110_5_B", "110_5_C", "110_7_A", "110_7_B", "110_7_C", 
                                          "110_9_A", "110_9_B", "110_9_C", "110_11_A", "110_11_B", "110_11_C"))
head(df.m)
base_size = 9
p <- ggplot(df.m, aes(Var2, Var1)) + 
  geom_tile(aes(fill = rescale), colour = "white") + 
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_grey(base_size = base_size) + labs(x = "", y = "") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  theme(legend.position = "none", 
        plot.margin=unit(c(0, 1, .2, .2),"cm"), #t, r, b, l
        axis.ticks = element_blank(), 
        axis.text.x = element_text(size = base_size * 0.8, angle = 330, hjust = 0, colour = "grey50"))
p

ggsave(paste("Heat_map_permanently_common", cutoff, "cutoff.pdf", sep="_"), width=14, height=14)



# The ratio of rare species -------------------------------------------------------------------------
binary <- readline("binary, TRUE or FALSE: ")

if (binary == "TRUE") {
  cat("take species number into consideration \n")
  truncated_ds_dominant.bi <- as.matrix((truncated_ds_dominant > 0) + 0)
  truncated_ds_rare_without_dominant.bi <- as.matrix((truncated_ds_rare > 0) + 0)
  t(truncated_ds_dominant.bi)[1:25, 1:3]
  t(truncated_ds_rare_without_dominant.bi)[1:25, 1:3]
  df <- data.frame(cbind(rowSums(truncated_ds_dominant.bi), rowSums(truncated_ds_rare_without_dominant.bi), rowSums(as.matrix((wholeDS > 0) + 0))))
  df$rareratio <- df[,2]/df[,3]
  cat("The number of rare species of each community on average is", round(mean(df$rareratio*100), 2), "+-",
      round(sd(df$rareratio)/sqrt(length(df$rareratio)), 2),"%")
} else {
  cat("take species abundance into consideration \n")
  t(truncated_ds_dominant)[1:25, 1:3]
  t(truncated_ds_rare_without_dominant)[1:25, 1:3]
  df <- data.frame(cbind(rowSums(truncated_ds_rare_without_dominant), rowSums(truncated_ds_rare_without_dominant), rowSums(wholeDS)))
  df$rareratio <- df[,2]/df[,3]
  cat("The abundance of rare species of each community on average is", round(mean(df$rareratio*100), 2), "+-",
      round(sd(df$rareratio)/sqrt(length(df$rareratio)), 2),"%")
}

group_info<-data.frame(row.names=rownames(df), t(as.data.frame(strsplit(rownames(df),"_"))))
head(group_info)

df2 <- data.frame(Year = as.factor(group_info[,2]),
                  Month = as.factor(group_info[,3]),
                  rareratio = df$rareratio)

row.names(df2) <- row.names(df)

df2$Year <- factor(df2$Year, levels=c("0", "10", "40", "70", "110"))

df2$Month <- factor(df2$Month, levels=c("5", "7", "9", "11"), 
                    labels=c("May", "July", "September", "November"))
str(df2); mean(df2$rareratio)

dstats<-function(x)(c(n=length(x), mean=mean(x), sd=sd(x)))
data<-summaryBy(rareratio~Year+Month, data=df2, FUN=dstats)
data$se <- data$rareratio.sd / sqrt(data$rareratio.n)
head(data)
str(data)

# bar plot of the ratio of rare species in the whole community
p <- ggplot(data, aes(x=Year, y=round(rareratio.mean*100, 2), fill=Month)) + 
  geom_bar(stat="identity", position="dodge", colour="black") +
  geom_errorbar(aes(ymin=round((rareratio.mean-se)*100, 2), ymax=round((rareratio.mean+se)*100, 2)), 
                colour="black", width=.3, position=position_dodge(0.9)) +
  scale_fill_manual(values=brewer.pal(n = 4, name = "RdBu")) +
  guides(fill = guide_legend(title = NULL))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 80)) +
  labs(title="",
       x = "Stage of succession (Years)", 
       y = "Proportion of Rare Biosphere (%)")+
  mytheme
p

ppi=300
png(paste("Rare_ratio", cutoff, "cutoff.png", sep="_"), width=6*ppi, height=4*ppi,res=ppi)
print(p)
dev.off()

# Venn Plot for the overlap of rare (or dominate) species between successional stages --------------------
# df <- wholeDS
# df <- truncated_ds_dominant
df <- truncated_ds_rare_without_dominant

group<-data.frame(row.names=rownames(df),
                  t(as.data.frame(strsplit(rownames(df),"_"))))
head(group)

df2 <- data.frame(df, Year=as.factor(group[,2]))
df2$Year <- factor(df2$Year, levels=c("0","10","40","70","110"), ordered=TRUE)

yr0 <- df2[df2$Year == '0', ]
yr0$Year <- NULL
yr0 <- yr0[, colSums(yr0)!=0] 
yr10 <- df2[df2$Year == '10', ]
yr10$Year <- NULL
yr10 <- yr10[, colSums(yr10)!=0] 
yr40 <- df2[df2$Year == '40', ]
yr40$Year <- NULL
yr40 <- yr40[, colSums(yr40)!=0] 
yr70 <- df2[df2$Year == '70', ]
yr70$Year <- NULL
yr70 <- yr70[, colSums(yr70)!=0] 
yr110 <- df2[df2$Year == '110', ]
yr110$Year <- NULL
yr110 <- yr110[, colSums(yr110)!=0] 

#colorRampPalette(brewer.pal(8, "Set1"))

venn.plot <- venn.diagram(x= list(yr0 = colnames(yr0), yr10 = colnames(yr10), 
                                  yr40 = colnames(yr40), yr70 = colnames(yr70), 
                                  yr110 = colnames(yr110)), 
                          filename = paste("venn_diagram_dominant", cutoff, "cutoff.png", sep="_"), 
                          imagetype="png", 
                          height = 2000, width = 2000, resolution = 300,
                          col="white", lwd=0.6, fill=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),
                          alpha = 0.50, 
                          fontfamily = "sans",
                          cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
                                  1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
                          cat.cex = 1.5,
                          cat.fontface = "bold",
                          cat.fontfamily = "sans",
                          #cat.dist = c(0.06, 0.1, 0.06, 0.06, 0.09),
                          category.names = c("0 year", "10 years", "40 years", "70 years", "110years"),
                          cat.col=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),
                          margin = 0.1)



#####################################
#Correlation and Procrustes calculations

res <- matrix(NA,3,3)	#create a matrix to store mantel and procrustes data
colnames(res)<-c("dataset","corrcoeff","R Procrustes")

# change truncated matrix by hand
matrix <- truncated_ds_dominant
matrix <- truncated_ds_rare_without_dominant
matrix <- truncated_ds_rare

ori_dist <- vegdist(wholeDS ,method="bray")	#distance matrix of the original dataset
tru_dist <- vegdist(matrix, distance="bray")	#distance matrix of the truncated dataset

cat('correlation coefficient = spearman or others?\n')
corcoef = readline("corcoef=? (spearman or others)...\t\n")

ori_tru_cor <- cor.test(ori_dist, tru_dist, method=corcoef) #correlation between matrices

require(MASS)
ori_NMDS <- monoMDS(ori_dist)	#NMDS for original dataset ?isoMDS
tru_NMDS <- monoMDS(tru_dist)		#NMDS for truncated dataset
ori_tru_procrustes <- protest(ori_NMDS,tru_NMDS) #procrustes
summary(ori_tru_procrustes)
str(ori_tru_procrustes)
plot(ori_tru_procrustes)


# change row numer by hand
res[3,1] <- 'truncated_ds_rare'
res[3,2]<-cbind(ori_tru_cor$estimate)
res[3,3]<-cbind(ori_tru_procrustes$t0)
res <- as.data.frame(res)
write.csv(res, "correlation_bray_each_fraction.csv")
str(res)



#### Do not run this part any more, 
# rare species composition in Phyla level
taxa <- feature_table[, -c(1:60)]

name ='truncated_ds_rare_without_dominant'
#com <- cond_rare
#com <- t(truncated_ds_dominant)
com <- t(truncated_ds_rare_without_dominant)
#com <- truncated_ds_rare
com[1:4, 1:8]

# add taxonomy info
com <- merge(taxa, com, by="row.names")  # merge by feature id
row.names(com) <- com$Row.names
com <- com[, -c(1, 2, 4:8)]
str(com)
levels(com$Phylum)

phyla <- aggregate(com[,-1], list(com$Phylum), sum)
row.names(phyla) <- phyla$Group.1
phyla <- phyla[, -1]
phyla <- 100*phyla/rowSums(wholeDS)
phyla$mean <- rowMeans(phyla)
phyla <- subset(phyla, mean > 0.01) 
phyla <- phyla[with(phyla, order(mean)), ] #for stacked-bar plot order(mean), for line plot order(-mean)
phyla$mean <- NULL
phyla <- t(phyla)
colnames(phyla)

group_info<-data.frame(row.names=rownames(phyla), t(as.data.frame(strsplit(rownames(phyla),"_"))))
head(group_info)

phyla2 <- data.frame(Year=as.factor(group_info[,2]),
                     Month=as.factor(group_info[,3]),
                     replicates=as.factor(group_info[,4]),
                     phyla)
row.names(phyla2) <- row.names(phyla)
head(phyla2)

phyla3 <- melt(phyla2,id=c("Year","Month","replicates"), variable = 'Phyla')

phyla3$Year <- factor(phyla3$Year, levels=c("0", "10", "40", "70", "110"))

phyla3$Month <- factor(phyla3$Month, levels=c("5", "7", "9", "11"), 
                       labels=c("May", "July", "September", "November"))

names(table(phyla3$Phyla))
str(phyla3)

dstats<-function(x)(c(n=length(x), mean=mean(x), sd=sd(x)))
data<-summaryBy(value~Year+Month+Phyla, data=phyla3, FUN=dstats)
data$se <- data$value.sd / sqrt(data$value.n)
head(data)

# stacked-bar plot
f <- ggplot(data, aes(x=interaction(Month, Year), y=value.mean, fill=Phyla)) + 
  geom_bar(stat="identity", width=0.8) +
  scale_fill_manual(values = c("#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4", "#fed9a6",
                               "#d9d9d9", "#b3de69", "#fdb462", "#80b1d3", "#fb8072", 
                               "#bc80bd", "#ccebc5", "#ffed6f", "#bebada", "#ffffb3", 
                               "#8dd3c7", "#b15928", "#ffff99", "#6a3d9a", "#cab2d6", 
                               "#ff7f00", "#fdbf6f", "#33a02c", "#b2df8a", "#e31a1c", 
                               "#fb9a99", "#1f78b4", "#a6cee3",   "#ff7f00", "#fdbf6f", "#33a02c", "#b2df8a"))+ 
  # 
  labs(x=" ",y="Relative abundance (%)", title=" ") +
  theme_classic()+
  theme(axis.text.x=element_text(size=rel(1.2), angle=45, hjust=1, vjust=1),
        axis.title.y = element_text(size = rel(1.2), angle = 90),
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.title = element_text(face = "bold"))
f

ppi=300
png(paste("Phylum_composition", name, cutoff, "cutoff.png", sep="_"), width=9*ppi, height=6*ppi,res=ppi)
print(f)
dev.off()


# line plot
# move move errorbars .05 to the left and right to avoid overlap
head(data)
pd <- position_dodge(0.2)  

mytheme <- theme_classic()+
  theme(legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.title = element_text(face = "bold")) 

p <- ggplot(data, aes(x=Year, y=value.mean, group=Month, colour=Month, shape=Month)) + 
  geom_errorbar(aes(ymin=value.mean-se, ymax=value.mean+se), colour="gray", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.5) +
  scale_shape_manual(values=c(15,16,17,23)) +
  scale_color_brewer(palette="Dark2") +
  facet_wrap(~Phyla, nrow=5, scales="free_y") +
  labs(title=" ",
       x = "Stage of succession (Year)", 
       y = "Relative abundance (%)")+
  mytheme
p

ppi=300
png(paste("Phylum_abundance", name, cutoff, "cutoff.png", sep="_"), width=9*ppi, height=9*ppi,res=ppi)
print(p)
dev.off()

pdf(paste("Phylum_abundance", name, cutoff, "cutoff.pdf", sep="_"), width=10, height=10)
print(p)
dev.off()

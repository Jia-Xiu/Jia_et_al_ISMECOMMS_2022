# taxa summarize for the entire community
# Author: Jia Xiu 
# Date: 17-04-2019 

rm(list=ls())

# Load the directory
directory = '~/Dropbox/' 
subfolder = 'Schier/RNA'

setwd(paste(directory, subfolder, sep="/"))
getwd()

library(reshape2)
library(vegan)
library(ape) # for pcoa 
library(ggplot2)
library(RColorBrewer)
library(doBy)
#display.brewer.all()


# taxonomy summarize --------------------------------------------------------------------------------------

# Phyla
com <- read.csv("schier_cdna_feature_table_dada2_silva_rarefied_taxonomy.csv", sep=",", 
                header=1, row.names=1, na.strings=c(""," ","NA"))
com <- com[, c(63, 1:60)]
com$Phylum <- gsub("Candidate division ", "", com$Phylum)
com$Phylum <- factor(com$Phylum)
com[1:4,1:9]
str(com)
levels(com$Phylum)

(rarefactiondepth <- mean(colSums(com[,2:60])))

phyla <- aggregate(com[,-1], list(com$Phylum), sum)
row.names(phyla) <- phyla$Group.1
phyla <- phyla[, -1]
phyla <- 100*phyla/rarefactiondepth
phyla$mean <- rowMeans(phyla)
#phyla <- subset(phyla, mean > 0.2) 
#phyla <- phyla[with(phyla, order(phyla$cDNA_110_11_A)), ]
phyla <- phyla[with(phyla, order(mean)), ] #for stacked-bar plot order(mean), for line plot order(-mean)
phyla$mean <- NULL
phyla <- t(phyla)
phyla <- as.data.frame(phyla)
#phyla$Others <- 100-rowSums(phyla)
#phyla <- phyla[, c(17, 1:16)]
colnames(phyla)
str(phyla)
class(phyla)

# split treatment info
group_info<-data.frame(row.names=rownames(phyla),
                       t(as.data.frame(strsplit(rownames(phyla),"_"))))
head(group_info)

# combine treatment info with phyla relative abundance
phyla2 <- data.frame(Year=as.factor(group_info[,2]),
                  Month=as.factor(group_info[,3]),
                  replicates=as.factor(group_info[,4]),
                  phyla)
row.names(phyla2) <- row.names(phyla)
head(phyla2)

phyla3 <- melt(phyla2,id=c("Year","Month","replicates"), variable = 'Phyla')

phyla3$Year <- factor(phyla3$Year, levels=c("0", "10", "40", "70", "110"))

phyla3$Month <- factor(phyla3$Month, levels=c("5", "7", "9", "11"), 
                   #labels=c("M", "J", "S", "N"))
                   labels=c("May", "Jul", "Sep", "Nov"))

names(table(phyla3$Phyla))
#phyla3$Phyla <- factor(phyla3$Phyla, levels(phyla3$Phyla)[c(1:18, 20, 19)])
str(phyla3)


dstats<-function(x)(c(n=length(x), mean=mean(x), sd=sd(x), se=sd(x)/sqrt(length(x))))
data<-summaryBy(value~Year+Month+Phyla, data=phyla3, FUN=dstats)
head(data)

(colourCount = length(unique(data$Phyla)))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
set.seed(111); col=sample(col_vector, colourCount); pie(rep(1,colourCount), col=sample(col_vector, colourCount))

# stacked-bar plot
f <- ggplot(data, aes(x=Month, y=value.mean, fill=Phyla)) + 
  geom_bar(stat="identity", colour="black", width=0.8) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,100))+
  facet_grid(~Year, switch = "x", scales = "free_x") +
  scale_fill_manual(values = col)+
  guides(fill=guide_legend(reverse = TRUE, title = NULL)) +
  labs(x=" ",y="Relative abundance (%)", title=" ") +
  theme_classic()+
  theme(text = element_text(size=15),
        legend.text=element_text(size=11),
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.box.background = element_rect(),
        legend.box.margin = margin(6, 6, 6, 6),
        legend.position = "bottom", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(1, 0.5, 1, 0.5),"cm")) 
f






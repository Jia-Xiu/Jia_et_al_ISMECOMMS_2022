# alpha- & beta-diversity analysis 
# Date: 15-10-2019
# Author: Jia Xiu 

rm(list=ls())

# load the directory
directory = '~/Dropbox/'
subfolder = 'Schier/cdna'

setwd(paste(directory, subfolder, sep="/"))
getwd()

# load packages
library(vegan)
library(ggplot2)
library(RColorBrewer) 
library(ggpubr)
library(reshape2)
library(ape)

mytheme<-theme_bw()+
  theme(text = element_text(size=15),
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

# load the rarefied feature table ----------------------------------------------------------------------------
com <- read.csv("schier_cdna_feature_table_dada2_silva_rarefied_taxonomy.csv", header=1, row.names=1)
com <- com[, c(1:60)]
com <- t(com)
str(com)


# alpha-diversity -----------------------------------------------------------------------------------------
Shannon <- diversity(com)
Simpson <- diversity(com, "simpson")
Richness <- specnumber(com) 
Chao <- estimateR(com)
Chao1 <- as.data.frame(t(Chao))

group_info<-data.frame(row.names=rownames(com), t(as.data.frame(strsplit(rownames(com),"_"))))

# read phylogenetic diversity, which calculate by picante
data.set.name = 'wholeDS'
phylo.div <- read.csv(paste(data.set.name, "_PD_SR.csv", sep=""), row.names=1)
head(phylo.div)

df <- data.frame(cbind(Richness, Chao1 = Chao1$S.chao1, Shannon, PD = phylo.div$PD, Simpson),
                 Year=as.factor(group_info[,2]),
                 Month=as.factor(group_info[,3]))

df$Year <- factor(df$Year, levels=c("0", "10", "40", "70", "110"))
#df$Year <- as.numeric(as.character(df$Year))

df$Month <- factor(df$Month, levels=c("5", "7", "9", "11"), 
                   labels=c("May", "July", "September", "November"))
str(df)
#write.csv(df, paste("alpha_diversity.csv", sep="_"))

# Two-way ANOVA -----------------------------------------------------------------------------------------
# Richness, Shannon, Simpson and PD
aov <- aov(Richness ~ Year * Month, data = df) # NO interaction effect
summary(aov)

# check the homogeneity of variances
bartlett.test(Richness ~ interaction(Year, Month), data=df)
# p > 0.05, suggesting that the variance across groups is not statistically significantly different.
plot(aov, 1)

# check normality
plot(aov, 2) # draw the correlation between a given sample and the normal distribution.

# Shapiro-Wilk test of normality for univariate
shapiro.test(df$Richness)
# p-value > 0.05, implying that the distribution of the data are not significantly different from normal distribution. 

# post hoc test
library(agricolae)
#LSD method
LSD.test(df$Richness, df$Year) 
TukeyHSD(aov, "Year")



# violin plot
df1 <- df[,c(1:2,6:7)]
df1 <- melt(df1, id=c("Year", "Month"))
str(df1)

f1 <- ggplot(df1, aes(x=Year, y=value, fill=Month))+
  geom_violin(trim = FALSE,  position=position_dodge(0.8)) +
  stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1),
               geom = "pointrange", color = "#990000", size = 0.4, position=position_dodge(0.8))+
  facet_wrap(. ~ variable, nrow=1)+
  scale_fill_brewer(palette="Set3")+
  labs(x="Stage of succession (Years)", y=" ", title=" ")+
  mytheme
f1

df1 <- df[,c(3:4, 6:7)]
df1 <- melt(df1, id=c("Year", "Month"))
df1$variable <- factor(df1$variable, levels = c("Shannon", "PD") ,labels = c("Shannon", "Phylogenetic diversity"))
str(df1)

f2 <- ggplot(df1, aes(x=Year, y=value, fill=Month))+
  geom_violin(trim = FALSE,  position=position_dodge(0.8)) +
  stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1),
               geom = "pointrange", color = "#990000", size = 0.4, position=position_dodge(0.8))+
  facet_wrap(. ~ variable, nrow=1, scale="free_y")+
  scale_fill_brewer(palette="Set3")+
  labs(x="Stage of succession (Years)", y=" ", title=" ")+
  mytheme 
f2

f <- ggarrange(f1, f2, common.legend = TRUE, legend = "right", ncol = 1, nrow = 2)
f

ppi=300
png(paste("alpha_diversity", data.set.name, "all.png", sep="_"), width=9*ppi, height=7*ppi, res=ppi)
print(f)
dev.off()


# rarefaction curve -------------------------------------------------------------------------------
df <- read.csv("schier_cdna_rarefaction.csv", header=1, row.names=1, sep=";")
df <- df[, c(1:100)]
df <- t(df)

group <- data.frame(row.names=rownames(df),
                       t(as.data.frame(strsplit(rownames(df),"_"))))
head(group)

df <- data.frame(df, depth=group[,1])

df2 <- melt(df, id=c("depth"))

library(doBy)
dstats <- function(x)(c(n=length(x), mean=mean(x), sd=sd(x), se=sd(x)/sqrt(length(x))))
df3 <- summaryBy(value~depth+variable, data=df2, FUN=dstats)
df3$se <- df3$value.sd / sqrt(df3$value.n)
group2 <- data.frame(row.names=paste(df3$variable, seq(from=1, to=length(df3$variable)), sep="_"),
                    t(as.data.frame(strsplit(as.character(df3$variable), "_"))))
head(group2)

depth <-  data.frame(t(as.data.frame(strsplit(as.character(df3$depth), ".", fixed = TRUE))))
depth[,2] <- as.numeric(as.character(depth[,2]))
head(depth)
str(depth)

df4 <- data.frame(df3[, c(1,4,6)],
                  Depth=depth$X2,
                  Year=as.factor(group2[,2]),
                  Month=as.factor(group2[,3]),
                  Replicates=as.factor(group2[,4]),
                  Sample=df3$variable)

df4$Year <- factor(df4$Year, levels=c("0", "10", "40", "70", "110"))

df4$Month <- factor(df4$Month, levels=c("5", "7", "9", "11"), 
                       labels=c("May", "July", "September", "November"))

df4$Sample <- factor(gsub('cDNA_', '', df4$Sample))

df4$Sample <- factor(df4$Sample, levels = c("0_5_A", "0_5_B", "0_5_C", "0_7_A", "0_7_B", "0_7_C", 
                                            "0_9_A", "0_9_B", "0_9_C", "0_11_A", "0_11_B", "0_11_C",
                                            "10_5_A", "10_5_B", "10_5_C", "10_7_A", "10_7_B", "10_7_C", 
                                            "10_9_A", "10_9_B", "10_9_C", "10_11_A", "10_11_B", "10_11_C",
                                            "40_5_A", "40_5_B", "40_5_C", "40_7_A", "40_7_B", "40_7_C", 
                                            "40_9_A", "40_9_B", "40_9_C", "40_11_A", "40_11_B", "40_11_C", 
                                            "70_5_A", "70_5_B", "70_5_C", "70_7_A", "70_7_B", "70_7_C", 
                                            "70_9_A", "70_9_B", "70_9_C", "70_11_A", "70_11_B", "70_11_C",
                                            "110_5_A", "110_5_B", "110_5_C", "110_7_A", "110_7_B", "110_7_C", 
                                            "110_9_A", "110_9_B", "110_9_C", "110_11_A", "110_11_B", "110_11_C"))
str(df4)

pd <- position_dodge(0.2)  

p <- ggplot(df4, aes(x=Depth, y=value.mean, colour=Sample)) + 
  geom_errorbar(aes(ymin=value.mean-value.se, ymax=value.mean+value.se), colour="gray", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=1.5) +
  labs(title=" ",
       x = "Number of sequence reads", 
       y = "Number of ASVs")+
  theme_bw()+
  theme(text = element_text(size=15),
        legend.position = "right",
        legend.box.background = element_rect(),
        legend.box.margin = margin(1, 1, 1, 1),
        legend.title = element_text(face = "bold", size = 10)) 
p

ggsave("rarefaction.pdf", width = 12, height = 7, units = "cm", p, scale = 2)


# beta-diversity -----------------------------------------------------------------------------------------

# Bray-Curtis or Jaccard
distance_method <- "bray"
#distance_method <- "jaccard"
dist <- vegdist(com, method=distance_method, diag=1) 

# PCoA
re <- pcoa(dist, correction="none", rn=NULL)
str(re)

group_info <- data.frame(row.names=rownames(re$vectors),t(as.data.frame(strsplit(rownames(re$vectors),"_"))))
head(group_info)

df <- data.frame(x=re$vectors[,1],y=re$vectors[,2],
                 Year=as.factor(group_info[,2]),
                 Month=as.factor(group_info[,3]),
                 replicates=as.factor(group_info[,4]))

df$Year <- factor(df$Year, levels=c("0", "10", "40", "70", "110"))

df$Month <- factor(df$Month, levels=c("5", "7", "9", "11"), 
                   labels=c("May", "July", "September", "November"))

(p1 <- ggplot(df, aes(x, y, shape=Month, fill=Year))+
  geom_point(size=5, alpha=0.7)+ #, shape=21
  labs(x=paste("PCoA1 (", round(re$values$Relative_eig[1] * 100, 2), "%)", sep=""), 
       y=paste("PCoA2 (", round(re$values$Relative_eig[2] * 100, 2), "%)", sep=""), 
       title="Jaccard")+
  scale_fill_brewer(palette="Set3", guide=guide_legend(override.aes = list(shape=21)))+
  scale_shape_manual(values=c(24, 22, 23, 21))+ 
  mytheme)

(p <- ggarrange(p1, p2, labels = c("A", "B"), 
               ncol = 2, nrow = 1 , common.legend = TRUE, legend = "right"))

ggsave("PCoA_jaccard_bray.png", width = 14, height = 6, units = "cm", p, scale = 1.5, dpi = 300)



#PERMANOVA ------------------------------------------------------------------------------------------------
str(dist)
row.names(as.matrix(dist))

df <- data.frame(row.names=row.names(as.matrix(dist)), t(as.data.frame(strsplit(row.names(as.matrix(dist)),"_"))))

df <- rename(df, c("X1"="Dataset", "X2"="Year", "X3"="Month", "X4"="replicates"))
df$Year <- factor(df$Year, levels=c("0", "10", "40", "70", "110"))

df$Month <- factor(df$Month, levels=c("5", "7", "9", "11"), 
                   labels=c("May", "July", "September", "November"))
head(df)

# two way permanova (successional stages vs. season)
result <- adonis(dist ~ Year*Month, data=df, permutation=999) # two way permanova
result
write.csv(result$aov.tab, paste("PERMANOVA_result_wholeDS_", distance_method,".csv", sep = ""))



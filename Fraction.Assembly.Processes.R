# Quantify the contribution of each assembly processes
# Created on 29-01-2018 by Jia Xiu
# Updated on 30-04-2019 by Jia Xiu

rm(list=ls())

# Load the libraries
library(vegan)
library(ggplot2)
library(RColorBrewer) 
library(reshape2) 
library(scales)
library(ggforce)
library(dplyr)


# load directory --------------------------------------------------------------------------------------------
directory = '~/Dropbox/'
subfolder = 'Schier/cDNA'

setwd(paste(directory, subfolder, sep="/"))
getwd()

# make a matrix to store the number of pairwise samples of each assembly process ----------------------------
df <- matrix(NA, nrow = 5, ncol = 3)
row.names(df) <- c("Variable.selection", "Homogeneous.selection", 
                   "Dispersal.limitation", "Homogenizing.dispersal", "Undominated.processes")
colnames(df) <- c("Whole", "Dominant", "Rare")
df

# set parameters --------------------------------------------------------------------------------------------
cutoff = 0.1/100
iteration = 999

# all three data sets
data.set.names = c('wholeDS', 'truncated_ds_dominant', 'truncated_ds_rare_without_dominant')

# a loop to calculte for the assembly processes of each rarity/commonnness type
for (xx in data.set.names) {
  data.set.name = xx
  
  # Quantify each assembly process ------------------------------------------------------------------------
  if (data.set.name == 'wholeDS') {
    nti <- read.csv(paste(data.set.name, "weighted_bNTI.csv", sep="_"), 
                    header=1, row.names=1, check.names=FALSE) 
    rc <- read.csv(paste("RC-bray", data.set.name, iteration, ".csv", sep="_"), 
                   header=1, row.names=1, check.names=FALSE)
  } 
  else if (data.set.name == 'truncated_ds_dominant') { 
    nti <- read.csv(paste(data.set.name, cutoff, "weighted_bNTI.csv", sep="_"), 
                    header=1, row.names=1, check.names=FALSE)
    rc <- read.csv(paste("RC-bray", data.set.name, cutoff, "999.csv", sep="_"), 
                   header=1, row.names=1, check.names=FALSE)
  } 
  else if (data.set.name == 'truncated_ds_rare_without_dominant') {
    nti <- read.csv(paste(data.set.name, cutoff, "weighted_bNTI.csv", sep="_"), 
                    header=1, row.names=1, check.names=FALSE)
    rc <- read.csv(paste("RC-bray", data.set.name, cutoff, "999.csv", sep="_"), 
                   header=1, row.names=1, check.names=FALSE)
  }
  
  # read weighted beta NTI -----------------------------------------------------------------------------------
  colnames(nti) <- sub("cDNA_", "", colnames(nti))
  row.names(nti) <- sub("cDNA_", "", row.names(nti))
  nti <- as.matrix(nti)
  # Function to extract pairwise value from a n*n lower trianglar matrix
  nti <- data.frame(as.table(nti))[lower.tri(nti, diag = FALSE), ]
  cat("should got:", (60*60-60)/2, "pair-wise distance\n"); 
  cat("Actually we got:", length(nti$Freq), "pair-wise distance\n")
  cat("the mean beta-NTI is:", round(mean(na.omit(nti$Freq)),2), "\n")
  row.names(nti) <- paste(nti$Var1, nti$Var2, sep = "_")
  head(nti)
  str(nti)
   
  # RC-bray -------------------------------------------------------------------------------------------------
  colnames(rc) <- sub("cDNA_", "", colnames(rc))
  row.names(rc) <- sub("cDNA_", "", row.names(rc))
  rc <- as.matrix(rc)
  rc <- data.frame(as.table(rc))[lower.tri(rc, diag = FALSE), ]
  cat("should got:", (60*60-60)/2, "pair-wise distance\n"); 
  cat("Actually we got:", length(rc$Freq), "pair-wise distance\n")
  cat("the mean RC-bray is:", round(mean(na.omit(rc$Freq)),2), "\n")
  row.names(rc) <- paste(rc$Var1, rc$Var2, sep = "_")
  head(rc)
  str(rc)
  
  
  # Combine the beta-NTI values with RC-bray
  nti.rc <- merge(nti, rc, by=0, all=TRUE)  # merge by row names (by=0 or by="row.names")
  nti.rc <- data.frame(nti = nti.rc$Freq.x, rc = nti.rc$Freq.y, row.names = nti.rc$Row.names)
  
  # Invalid the value of RC-bray in which the beta-NTI larger than +2 or less than -2
  for (i in 1:nrow(nti.rc)) {
    if (nti.rc[i,1] > 2 | nti.rc[i,1] < -2) {
      nti.rc[i, 2] <- NA
    }
  }
  
  head(nti.rc)
  str(nti.rc)
  
  # Quantify each assembly process ------------------------------------------------------------------------
  if (data.set.name == 'wholeDS') {
    i = 1 } else if (data.set.name == 'truncated_ds_dominant') { 
      i = 2 } else if (data.set.name == 'truncated_ds_rare_without_dominant') {
        i = 3
      }
  
  # Variable selection
  Variable.selection <- nti.rc$nti > 2
  cat('Number of variable  selection:', table(Variable.selection)['TRUE'], 'within', nrow(nti.rc), 'pairwise samples\n')
  df[1, i] <- length(Variable.selection[Variable.selection == TRUE])
  
  # Homogenous selction
  Homogeneous.selection <- nti.rc$nti < -2
  cat('Number of homogenous selection:', table(Homogeneous.selection)['TRUE'], 'within', nrow(nti.rc), 'pairwise samples\n')
  df[2, i] <- table(Homogeneous.selection)['TRUE']#length(c[Homogeneous.selection == TRUE])
  
  # Dispersal limilation
  Dispersal.limilation <- na.omit(nti.rc$rc) > 0.95
  cat('Number of dispersal limitation:', table(Dispersal.limilation)['TRUE'], 'within', nrow(nti.rc), 'pairwise samples\n')
  df[3, i] <- length(Dispersal.limilation[Dispersal.limilation == TRUE])
  
  # Homogenizing dispersal
  Homogenizing.dispersal <- na.omit(nti.rc$rc) < -0.95
  cat('Number of homogenizing dispersal:', table(Homogenizing.dispersal)['TRUE'], 'within', nrow(nti.rc), 'pairwise samples\n')
  df[4, i] <- length(Homogenizing.dispersal[Homogenizing.dispersal == TRUE])
  
  # Undominated processes
  Undominated.processes <- na.omit(nti.rc$rc) <= 0.95 & na.omit(nti.rc$rc) >= -0.95
  cat('Number of Undominated processes:', table(Undominated.processes)['TRUE'], 'within', nrow(nti.rc), 'pairwise samples\n')
  df[5, i] <- length(Undominated.processes[Undominated.processes == TRUE])
}

# calculate relatice impacts of each process --------------------------------------------------------------------------
df1 <- melt(df)
df1$value <- round(df1$value*100/nrow(nti.rc), 2)
colnames(df1) <- c("Processes", "groups", "value")
df1[df1 == 0] <- NA
df1 <- na.omit(df1)

df1$Processes <- factor(df1$Processes, levels = c('Variable.selection', 'Homogeneous.selection', 'Dispersal.limitation', 
                                                'Homogenizing.dispersal', 'Undominated.processes'),
                       labels = c('Variable selection', 'Homogeneous selection', 'Dispersal limilation', 
                                  'Homogenizing dispersal', 'Undominated processes'))
df1$groups <- factor(df1$groups, levels = c('Whole', 'Dominant', 'Rare'), 
                       labels = c('Whole community', 'Common biosphere', 'Rare biosphere'))

df1

# first way to generate a pie plot ---------------------------------------------------------------------------------
# clockwise

# calculate the start and end angles for each pie
dat_pies <- left_join(df1,
                      df1 %>% 
                        group_by(groups) %>%
                        summarize(value_total = sum(value))) %>%
  group_by(groups) %>%
  mutate(end_angle = 2*pi*cumsum(value)/value_total,      # ending angle for each pie slice
         start_angle = lag(end_angle, default = 0),   # starting angle for each pie slice
         mid_angle = 0.5*(start_angle + end_angle))   # middle of each pie slice, for the text label

rpie = 1 # pie radius
rlabel = 0.6 * rpie # radius of the labels; a number slightly larger than 0.5 seems to work better, 0.5 would place it exactly in the middle as the question asks for.

# draw the pies
pie <- ggplot(dat_pies) + 
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = rpie, start = start_angle, end = end_angle, fill = Processes)) +
  geom_text(aes(x = rlabel*sin(mid_angle), y = rlabel*cos(mid_angle), label = paste(round(value,2), "%")), 
            hjust = 0.5, vjust = 0.5, size=4) +
  coord_fixed() +
  scale_x_continuous(limits = c(-1, 1), name = "", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-1, 1), name = "", breaks = NULL, labels = NULL) +
  facet_grid(.~groups)+
  scale_fill_brewer(palette="Pastel2", direction = -1)+
  theme_minimal()+
  theme(legend.title=element_text(size=13),
        text = element_text(size=15),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank())
pie

pdf(paste("Fraction_assembly_processes", cutoff, iteration, 'pie.pdf', sep = "_"),  width=10, height=5)
print(pie)
dev.off()


# the second way to generate a pie plot
# Anti-clockwise
pie <- ggplot(df1, aes(x="", y=value, fill=Processes))+
  geom_bar(stat="identity", width=1, colour = "black") +
  facet_grid(facets=. ~ groups)+
  coord_polar("y", start=0)+
  scale_fill_brewer(palette="Pastel2", direction = -1)+
  geom_text(aes(label = paste0(round(value, 2), "%")), position = position_stack(vjust = 0.5))+
  theme_minimal()+
  theme(text = element_text(size=14),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.ticks = element_blank(),
        plot.title=element_text(size=15, face="bold"))
pie




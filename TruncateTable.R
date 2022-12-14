# Based on Angelique Gobet & Alban Ramette, MultiCoLA Manual - Quick and Easy version (2011) 
# Author: https://www.mpi-bremen.de/Binaries/Binary1660/MultiCoLA.1.4.-.zip
# Modified by Xiu Jia 
# Date: 01-09-2018

# see new updates here https://github.com/Jia-Xiu/collaborations/blob/main/Araujo_et_al_2022/TruncateTable.R

# sample-based cutoff (Gobet at al., 2010)
TruncateTable<-function(dataset,cutoff,typem){
        # remove columns in the matrix for which the sum of the line is 0
        CLrow<-function(m) {
          m <- m[, colSums(m)!=0] 
          return(m)
         }

        CLcol<-function(m) {
          # remove all the rows in the matrix for which the sum of the line is 0
          # the case happend only if you increased the cutoff untill it reaches the lowest number of the maximum OTU occurance in all samples
          m <- m[, rowSums(m)!=0]
          return (m)
         }

			# Application of a percentage cut-off to the original dataset to obtain abundant/rare dataset
				Q1<-dataset
				if(typem == "dominant") {Q1[Q1/rowSums(Q1) <= cutoff] <- 0}	# all species presents no less than j times =0
				if(typem == "rare") {Q1[Q1/rowSums(Q1) > cutoff] <- 0}	# all species presents more than j times =0
				Q3 <- CLcol(CLrow(Q1))	#remove rows and columns whose sum=0
return(Q3)
} #end TruncateTable


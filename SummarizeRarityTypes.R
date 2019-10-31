# Function for summarizing types of rarity
# Author: Xiu Jia
# Date: 09-04-2019

# set parameters
cat("Set the species*site metrix as binary data?
    e.g. take considering of the abundance of each rarity types? Type 'FALSE', set 'rarefactiondepth' if binary==FALSE
    Take considering of the species bumber of each rarity type? Type 'TRUE'")
binary <- as.logical(readline("Set the species*site metrix as binary data, TRUE or FALSE\t"))

SummarizeRarityTypes <- function(dataset, rare_classified, rarefactiondepth, binary) {
  if (binary == FALSE) {
    print("binary == FALSE, take abundance into consideration")
    rare_classified <- transform(merge(as.data.frame(t(dataset)),
                                       as.data.frame(rare_classified$vector_category, row.names = row.names(rare_classified)), 
                                       by='row.names'), row.names=Row.names, Row.names=NULL)  
  } else {
    print("binary = TRUE, take present-absent into consideration")
    dataset <- as.data.frame(decostand(t(dataset), "pa"))
    rare_classified <- transform(merge(dataset,
                                       as.data.frame(rare_classified$vector_category, row.names = row.names(rare_classified)), 
                                       by='row.names'), row.names=Row.names, Row.names=NULL)  
  }
  colnames(rare_classified)[61] <- "rarity"
  rare_classified$rarity <- factor(rare_classified$rarity, levels=c("A", "B", "C"),
                                   labels=c("PermanentlyRare", "TransientlyRare", "ConditionallyRare"))
  rare_type_sum <- aggregate(rare_classified[,-61], list(rare_classified$rarity), sum)
  row.names(rare_type_sum) <- rare_type_sum$Group.1
  rare_type_sum <- rare_type_sum[, -1]
  rare_type_sum <- t(rare_type_sum)
  print("structure of summarized rarity types")
  str(rare_type_sum)
  
  # relative abundance of each type of rarity or number of ASVs of each type of rarity
  if  (binary == FALSE) {
    print("Summerizing the relative abundance of each type of rarity")
    name <- "rare_type_sum"
    df <- 100*rare_type_sum/rarefactiondepth
    cat("the average relative abundance of each type of rarity:\n", 
        "Permanently Rare is", round(mean(df[, 1]/rowSums(df))*100,2), "+-", 
        round(sd(df[, 1]*100/rowSums(df))/sqrt(nrow(df)), 2), "%\n", 
        "Transiently Rare is", round(mean(df[, 2]/rowSums(df))*100,2), "+-", 
        round(sd(df[, 2]*100/rowSums(df))/sqrt(nrow(df)), 2),  "%\n", 
        "Conditionally Rare is", round(mean(df[, 3]/rowSums(df))*100,2), "+-", 
        round(sd(df[, 3]*100/rowSums(df))/sqrt(nrow(df)), 2),  "%\n")
    cat("average relative abundance of rare biosphere", round(mean(rowSums(df)),2), "%")
  } else {
    print("Summerizing the number of ASVs of each type of rarity")
    name <- "rare_type_sum_binary"
    df <- rare_type_sum
  }
  return(df)
}

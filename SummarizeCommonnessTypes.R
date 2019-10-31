# Function for summarizing types of commonness
# Author: Xiu Jia
# Date: 09-04-2019

# set parameters

warning(cat("Check whether rarefactiondepth", rarefactiondepth, "is the value your mean?"))

cat("Set the species*site metrix as binary data?
    e.g. take considering of the abundance of each commonness types? Type 'FALSE', set 'rarefactiondepth' if binary==FALSE
    Take considering of the species bumber of each commonness type? Type 'TRUE'")
binary <- as.logical(readline("Set the species*site metrix as binary data, TRUE or FALSE\t"))

dominant_classified <- as.data.frame(t(truncated_ds_dominant))
#dominant_classified <- truncated_ds_dominant_classified


SummarizeCommonnessTypes <- function(truncated_ds_dominant, truncated_ds_rare_without_dominant, rarefactiondepth, binary) {
  if (binary == FALSE) {
    print("binary == FALSE, take abundance into consideration")
    dominant_classified <- as.data.frame(t(truncated_ds_dominant))
    dominant_classified$DominantType <- NA
    dominant_classified$DominantType <- ifelse(row.names(dominant_classified) %in% colnames(truncated_ds_rare_without_dominant), 
                                               "ConditionallyDominant", "PermanentlyDominant")
    dominant_classified <- dominant_classified
  } else {
    print("binary = TRUE, take present-absent into consideration")
    dominant_classified <- as.data.frame(t(truncated_ds_dominant))
    dominant_classified <- as.data.frame(decostand(dominant_classified, "pa"))
    dominant_classified$DominantType <- NA
    dominant_classified$DominantType <- ifelse(row.names(dominant_classified) %in% colnames(truncated_ds_rare_without_dominant), 
                                               "ConditionallyDominant", "PermanentlyDominant")
  }
  dominant_type_sum <- aggregate(dominant_classified[,-61], list(dominant_classified[,61]), sum)
  row.names(dominant_type_sum) <- dominant_type_sum$Group.1
  dominant_type_sum <- dominant_type_sum[, -1]
  dominant_type_sum <- t(dominant_type_sum)
  str(dominant_type_sum)

# Using SummarizeDominantTypes function to summarize dominant types
if  (binary == FALSE) {
  print("take considering of the abundance of each type of commonness")
  df <- 100*dominant_type_sum/rarefactiondepth
  cat("average abundance of dominant biosphere", round(mean(rowSums(df)), 2), "%")
  cat("the average relative abundance of each dominant types\n", 
      "conditionally dominant is", round(mean(df[, 1]/rowSums(df))*100,2), "+-", 
      round(sd(df[, 1]*100/rowSums(df))/sqrt(nrow(df)), 2), "%\n", 
      "permanently dominant is", round(mean(df[, 2]/rowSums(df))*100,2), "+-", 
      round(sd(df[, 2]*100/rowSums(df))/sqrt(nrow(df)), 2), "%\n")
} else {
  print("Take considering of the ASVs numer of each type of commonness")
  df <- dominant_type_sum
}
 return(df)
}

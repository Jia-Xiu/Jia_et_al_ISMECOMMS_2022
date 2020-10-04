# Classify different types of rarity or commonness
# Author: Xiu Jia
# E-mail: xibeihenai@gmail.com
# Builted on 26-11-2018
# Updated on 04-10-2020
# Acknowledge: Leonel Herrera Alsina

warning("This function only applies to rarified feature table")

# @my_data is the feature/OTU table with ASVs/OTUs (~species) in rows, and samples in columns
# @cutoff is a value of relative abundance below which species defined as rare, otherwise common
# @total_abundance is the value used to rarefy feature/OTU table

rarity_type <- function(my_data, cutoff, total_abundance) {
  
  my_data <- as.data.frame(my_data)
  new_table <- my_data
  new_table$vector_category <- NA
  
  for (i in 1:nrow(my_data)) {
    focal_species <- as.numeric(my_data[i, 1:ncol(my_data)])
    
    if (sum(focal_species) == 0) {
      cat("The ASV", i, "has abundance of zero \n")
      stop()
    }
    
    if (min(focal_species[focal_species > 0]) > cutoff * total_abundance) {
      # permanently common
      category <- 'PC'
      
    } else {
      
      if (max(focal_species) > cutoff * total_abundance &
          min(focal_species[focal_species > 0]) <= cutoff * total_abundance) {
        # conditionally rare or common
        category <- 'CRC'
        
      } else {
        
        if (sum(focal_species != 0) == 1) {
          # transiently rare
          category <- 'TR'
          
        } else {
          # permanently rare
          category <- 'PR'
        }
      }
    }
    
    new_table[i, 'vector_category'] <- category
  }
  
  return(new_table)
}

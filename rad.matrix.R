# generate a rank abundance matrix for plotting rank abundance curve
# built by Jia, Xiu on 28-01-2019

# mock data which you can use to test this function
# x <- matrix(sample(rep(1:9, each = 2, len = 50)), nrow = 10, dimnames = list(c(letters[1:10]), c(LETTERS[1:5])))
# x[x>7] <- 0; x<- as.data.frame(x); x

# for a feature table that rows as species, columns as samples
rad.matrix <- function (wholeDS) {
  
  # find sample that has the maximum species number
  max.nrow <- max(colSums(wholeDS != 0))
  
  # create a matrix with rows based on sample that has the maximum number of species, and columns as species 
  rad.df <- matrix(nrow = max.nrow, ncol = 1 + ncol(wholeDS))
  colnames(rad.df) <- c("rank", colnames(wholeDS))
  
  # reorder species by their abundance per sample
  for(i in 1:ncol(wholeDS)){
    subset <- wholeDS[ ,i]
    names(subset) <- row.names(wholeDS)
    take <- subset > 0
    nm <- names(subset) # names of each species, but not important
    rad.vector <- subset[take]
    names(rad.vector) <- nm[take]
    
    # order taxa
    rad.vector <- sort(rad.vector, decreasing = TRUE) 
    length(rad.vector) <- max.nrow
    rad.df[, 1] <- rep(1:max.nrow)
    rad.df[, i+1] <- rad.vector
    rad.df[is.na(rad.df)] <- 0
    rad.df <- as.data.frame(rad.df)
    
  }
  
  return(rad.df)
}

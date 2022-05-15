## parameter num.dim for number of singular components
## parameter weighted for kurtosis weights in continuous variables
## parameter ntrees for isolation forest algorithm 

FAMDAD <- function(data, weighted = TRUE, num.dim = NULL, cutoff = 0.5,  ntrees = 500) {
  require(e1071)
  require(isotree)
  require(FactoMineR)
  
  ## additional function, input: factor data
  OHE <- function(data) {
    data <- as.data.frame(data)
    new.data <- data.frame(matrix(0, nrow = nrow(data), ncol = length(unlist(sapply(data, levels)))))
    len <- function(i, cat.data) {
      return(length(levels(cat.data[, i])))
    }
    colnames(new.data) <- paste(rep(names(data), sapply(1:ncol(data), len, cat.data = data)), (unlist(sapply(data, levels))), sep="_") 
    colN = 0
    for (i in 1:ncol(data)) {
      m =  length(levels(data[, i]))
      for (j in 1:m) {
        new.data[which(data[, i] ==levels(data[, i])[j]), (colN+j)] <- 1  
      }
      colN = colN + m 
    }
    return(new.data)
  }

  Y <- OHE(data[, sapply(data, is.factor)])
  p <- colSums(Y)/nrow(data)
  Zd <- as.data.frame(t(t(Y) / p) - 1)
  Zc <- scale(data[, sapply(data, is.numeric)])
  Z <- cbind(Zd, Zc)
  
  if (weighted) {
    kapa <- apply(data[, sapply(data, is.numeric)], 2, kurtosis, na.rm = TRUE)
    W <- c(p,ifelse(kapa <=0, 1, kapa))  
  } else {
    W <- c(p, rep(1,  sum(sapply(data, is.numeric))))
  }
  B <- rep(1/nrow(data), nrow(data))
  SVD <- svd.triplet(as.matrix(Z), row.w = B, col.w = W)
  if(is.null(num.dim)) {
    num.dim = sum(SVD$vs >= 1)
  }
  FA <- data.frame(as.matrix(Z)%*%(diag(sqrt(W)))%*%SVD$V[, 1:num.dim]) 
  iso <- isolation.forest(FA,  output_score = TRUE, output_dist = FALSE, output_imputations = FALSE, ntrees = ntrees)
  

  outliers <- rep(0, nrow(data))
  outliers[which(iso$scores > cutoff)] <- 1
  return(list(scores = iso$scores, outliers = outliers))
}

##########################################################################################################################################################


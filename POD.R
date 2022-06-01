## Outlier detection method for Mixed-attribute Data 
## Proposed in "An Effective Pattern Based Outlier Detection Approach for Mixed Attribute Data" (DOI 10.1007/978-3-642-17432-2_13)
#############################################################################################################################################################
## Inputs: data frame data with observations in rows and attributes in columns 
## parameter n: k_nn paremeter for k-nearest neighbor

POD <- function(data, n = round(0.01*nrow(data)), prior = FALSE) {  
  require(spatstat.geom)
  
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
  
  ## additional function
  pi <- function(vec, model) {
    return(1/(1+ exp(model$coefficients%*%vec)))
  }
  outliers <- rep(0, nrow(data))
  
  cont.data <- data[, sapply(data, is.numeric)]
  cat.data <- OHE(data[, sapply(data, is.factor)])
  COF <- data.frame(matrix(0, nrow = nrow(data), ncol = ncol(cat.data)))
  names(COF) <- names(cat.data)
  
  for (i in 1:ncol(cat.data)) {
    glm.data = data.frame(y= cat.data[, i], cont.data)
    if (prior) {
      prior.weights = ifelse(cat.data[, i]==0, mean(cat.data[, i]), 1-mean(cat.data[, i]))
      model <- glm(y~., family = binomial(link = "logit"), data = glm.data, weights = prior.weights)
    } else {
      model <- glm(y~., family = binomial(link = "logit"), data = glm.data)
    }
    pis <- apply(cbind(1, cont.data), 1, pi, model = model)
   
    COF[which(cat.data[, i]==1), i] <- 1 - 0.5*pis[which(cat.data[, i]==1)]
    COF[which(cat.data[, i]==0), i] <- 0.5 + 0.5*pis[which(cat.data[, i]==0)]
  }

  COF$NN <- exp(1+ nndist(scale(cont.data), k = round(0.01*nrow(data))))
  if(sum(is.infinite(COF$NN)) > 0) {
    COF$NN[which(is.infinite(COF$NN))] <- exp(1)
  }
  MADOF <- sqrt(rowSums(COF^2))
  scores <- (MADOF - min(MADOF))/(max(MADOF) - min(MADOF))
  outliers[order(MADOF, decreasing =TRUE)[1:n]] <- 1
  return(list(MADOF= MADOF, outliers= outliers, scores = scores))
}

##############################################################################################################################################

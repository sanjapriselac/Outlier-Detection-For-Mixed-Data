## paramters k for the number of discretization interval
## parameter method = {"ZDisc", "KMeans"}

ZDisc_KMeans <- function(data, method = "ZDisc", k = round(log2(nrow(data)+1)), cutoff =0.5) { 
  ZDisc <- function(cont.data, k) {
    cont.data <- as.data.frame(scale(cont.data))
    disc.data <- cont.data
    for (j in 1:ncol(cont.data)) { 
      a = min(cont.data[, j],  na.rm = TRUE)
      b = max(cont.data[, j],  na.rm = TRUE)
      d = (b-a)/k
      for (i in 1:(k-1)) {
        disc.data[which((cont.data[, j] >= (a+(i-1)*d)) & (cont.data[, j] < (a+i*d))), j] <- paste0("[", (a+(i-1)*d), ",",(a+i*d), ")") 
      }
      disc.data[which((cont.data[, j] >= (a+(k-1)*d)) & (cont.data[, j] <= b)), j] <- paste0("[", (a+(k-1)*d), ",",b, ")") 
      disc.data[, j] <- as.factor(disc.data[, j])
    }
    return(disc.data)
  }
  
  K_means <- function(cont.data, k) {
    for (j in 1:ncol(cont.data)) {
      k.means <- kmeans(cont.data[, j], k)
      cont.data[, j] <- factor(k.means$cluster)
    }
    return(cont.data)
  }
  
  AVF <- function(data) {
    scores <- rep(0, nrow(data))
    for (at in 1:ncol(data)) {
      for(level in levels(data[, at])) {
        observations <- which(data[, at] == level)
        scores[observations] <- scores[observations] + length(observations)
      } 
    }
    scores <- scores/ncol(data)
    return(scores)
  }
  
  cat.data <- data[, sapply(data, is.factor)]
  cont.data <- data[, sapply(data, is.numeric)]
  if (method == "ZDisc") {
    ZScore_data <- cbind(cat.data, ZDisc(cont.data, k=k))
    scores <- AVF(ZScore_data)
    scores <- 1 -(scores - min(scores, na.rm = TRUE))/(max(scores, na.rm = TRUE) - min(scores, na.rm = TRUE))
    outliers <- rep(0, nrow(data))
    outliers[which(scores > cutoff)] <- 1
  } 
  if (method == "KMeans"){
    KMeans_data <- cbind(cat.data, K_means(cont.data, k=k))
    scores <- AVF(KMeans_data)
    scores <- 1 -(scores - min(scores, na.rm = TRUE))/(max(scores, na.rm = TRUE) - min(scores, na.rm = TRUE))
    outliers <- rep(0, nrow(data))
    outliers[which(scores > cutoff)] <- 1
  }
  return(list(scores = scores, outliers = outliers))
}

################################################################################################################################################################


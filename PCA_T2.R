## parameter k for the number of principal components

PCAmix_T2 <- function(data, alpha = 0.05, k = NULL) {
  require(PCAmixdata)
  
  ndim = sum(sapply(data, is.numeric)) + length(unlist(sapply(data[, sapply(data, is.factor)], levels)))
  pca <- PCAmix(data[, sapply(data, is.numeric)], data[, sapply(data, is.factor)], rename.level = TRUE, ndim = ndim, graph = FALSE)
  if (is.null(k)) {
    k = sum((pca$eig[, "Cumulative"] <= 80))    ## number of components that explains more than 80%
  }
  
  t2 <- function(vec, pca, mu, k) {
    return(t(vec-mu)%*%diag(1/pca$eig[1:k, "Eigenvalue"])%*%(vec-mu))
  }
  T2 <- apply(pca$scores[, 1:k], 1, t2, pca = pca, mu = colMeans(pca$scores)[1:k], k=k)
  n = nrow(data)
  Conv.Control <- (k*(n+1)*(n-1)/(n^2 - n*k))*qf(1 - alpha, k, n-k)
  outliers <- rep(0, n)
  outliers[which(T2 > Conv.Control)] <- 1
  scores <- (T2 - min(T2))/(max(T2) - min(T2))
  return(list(T2 = T2, outliers = outliers, scores = scores))
}

######################################################################################################################################################################


## parameters min_sup, max_len for Apriori algorithm
## parameters mc for minimum count, pbar = {"uniform", "one"} 

ABOD <- function(data, min_sup = 0.1, mc = round(0.4*nrow(data)), pbar = "one", k = round(log2(nrow(data)+1)), max_len = 7) { 
  require(arules)
  
  ZDisc <- function(cont.data, k) {
    cont.data <- as.data.frame(scale(cont.data))
    disc.data <- cont.data
    for (j in 1:ncol(cont.data)) { 
      a = min(cont.data[, j],  na.rm = TRUE)
      b = max(cont.data[, j],  na.rm = TRUE)
      d = (b-a)/k
      for (i in 1:(k-1)) {
        disc.data[which((cont.data[, j] >= (a+(i-1)*d)) & (cont.data[, j] < (a+i*d))), j] <- paste0("[", round((a+(i-1)*d), 3), ",",round((a+i*d), 3), ")") 
      }
      disc.data[which((cont.data[, j] >= (a+(k-1)*d)) & (cont.data[, j] <= b)), j] <- paste0("[", round((a+(k-1)*d), 3), ",",round(b, 3), ")") 
      disc.data[, j] <- as.factor(disc.data[, j])
    }
    return(disc.data)
  }
  
  get.associated <- function(data, min_sup, mc, max_len) {
    tic()
    freq.itemset <- apriori(data,  parameter = list(target = "frequent itemsets", support = min_sup, maxlen = max_len), control = list(verbose = FALSE))
    toc()
    freq.itemset <- DATAFRAME(freq.itemset, setStart = '', itemSep = ', ', setEnd = '')
    freq.itemset$items <- as.character(freq.itemset$items)
    ## remove itemsets with len 1 (1 attribute)
    freq.itemset <- freq.itemset[sapply(freq.itemset$items, grepl, pattern= " "),]
  
    if (pbar == "uniform") {
      associated.attributes <-  data.frame(matrix(ncol = 3, nrow = 0))
      names(associated.attributes) <- c("attributes", "count", "combinations")
      
      for (i in 1:nrow(freq.itemset)) {
        attr.split <- unlist(strsplit(freq.itemset$items[i], split = ", "))
        attributes <-  c()
        combinations = 1
        for (j in 1:length(attr.split)) {
          attr.value <- unlist(strsplit(attr.split[j], split = "="))
          attributes <- paste(attributes, attr.value[1], sep = " ")
          combinations = combinations * length(levels(data[, attr.value[1]]))
        }
        attributes <- substring(attributes, 2)
        if (!(attributes %in% associated.attributes$attributes)) {
          j = nrow(associated.attributes) + 1
          associated.attributes[j,] = c(attributes, 1, combinations)
          associated.attributes$count <- as.numeric(associated.attributes$count)
          associated.attributes$combinations <- as.numeric(associated.attributes$combinations)
        } else {
          associated.attributes$count[which(associated.attributes$attributes == attributes)] <- associated.attributes$count[which(associated.attributes$attributes == attributes)] + 1
        }
      } 
    } else {
      associated.attributes <-  data.frame(matrix(ncol = 2, nrow = 0))
      names(associated.attributes) <- c("attributes", "count")
      
      for (i in 1:nrow(freq.itemset)) {
        attr.split <- unlist(strsplit(freq.itemset$items[i], split = ", "))
        attributes <-  paste(unlist(strsplit(attr.split, split = "="))[seq(1, 2*length(attr.split), 2)], collapse = ' ')
        
        if (!(attributes %in% associated.attributes$attributes)) {
          j = nrow(associated.attributes) + 1
          associated.attributes[j, "attributes"] = attributes
          associated.attributes[j, "count"] = 1
        } else {
          associated.attributes[which(associated.attributes$attributes == attributes), "count"] <- associated.attributes[which(associated.attributes$attributes == attributes), "count"] + 1
        }
      }
    }
    print(nrow(associated.attributes))
    if (sum(associated.attributes$count > mc) > 0) {
      associated.attributes <- associated.attributes[which(associated.attributes$count > mc), ]
    }
    print(nrow(associated.attributes))
    return(associated.attributes)
  }
  
  if (sum(sapply(data, is.numeric)) != 0) {
    if (sum(sapply(data, is.factor)) != 1) {
      data = cbind(data[, sapply(data, is.factor)], ZDisc(data[, sapply(data, is.numeric)], k = k))
    } else {
      data = cbind(data[, sapply(data, is.factor)], ZDisc(data[, sapply(data, is.numeric)], k = k))
      names(data)[-grep("X\\d",names(data))] <- "X1"
    }
    
  }
  ass.att <- get.associated(data, min_sup = min_sup, mc = mc, max_len = max_len)
 
  sds <- function(i, ass.att, data, n) {
    fun <- function(j, ass.att, n) {
      atts <- unlist(strsplit(ass.att$attributes[j], split = " "))
      p <- sum(sapply(data[, atts], as.integer) == as.integer(data[i, atts]))/n
      return(-log2(p)/p)
    }
    return(sum(sapply(1:nrow(ass.att), fun, ass.att = ass.att, n = n)))
  }
  
  n = nrow(data)
  tic()
  SD <- sapply(1:n, sds, ass.att = ass.att, data = data, n = n)
  toc()
  
  SD <- (SD - min(SD))/(max(SD) - min(SD)) 
  return(list(scores = SD))
}


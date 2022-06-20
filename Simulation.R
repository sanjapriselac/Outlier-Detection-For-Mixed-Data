## Sanja Priselac
## last edit on 20/06/2022
## Syhtetic Data Generation 

library(MASS)
library(dplyr)
library("scales")
library("reticulate")
library(pROC)
library(tictoc)
library(extraDistr)
library('reshape2')
library(gridExtra)
library('latex2exp')

#############################################################################################################################################################
## produce data set with the characteristics
#############################################################################################################################################################

## n: # observations, 
## p: # variables, 
## outn: outlier proportion, 
## cat: # of categorical variables, 
## rho: correlation coeff , 
## sd = sd for the outlier group, 
## k: # intervals for categorization 
## out.tyoe: "cat"/"cont"/"both", are the outlying points only in categorical/continious space or in both
## ZDisc for categorization 

produce.data <- function(n, p, cat, outn, rho, mu, out.type = "both", k = round(log2(n+1)), seed = 100) {
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
  
  set.seed(seed)
  sig <- matrix(NA, p, p)
  sig <- rho^abs(row(sig) - col(sig))
  
  x <- data.frame(mvrnorm(n=n, mu = rep(0, p), Sigma = sig))

  if (out.type == "cat") {
    x[(round(n*(1-outn))+1):nrow(x), 1:cat] <-  x[(round(n*(1-outn))+1):nrow(x), 1:cat]  + mu
  } 
  if (out.type == "cont") {
    x[(round(n*(1-outn))+1):nrow(x), (cat+1):ncol(x)] <-   x[(round(n*(1-outn))+1):nrow(x), (cat+1):ncol(x)] + mu
  } 
  if (out.type == "both") {
    x[(round(n*(1-outn))+1):nrow(x), ] <-  x[(round(n*(1-outn))+1):nrow(x), ] + mu
  }
 
  data = cbind(ZDisc(x[, 1:cat], k = k), x[,(cat+1):ncol(x)])
  labels <- rep(c(0, 1), c(round(n*(1-outn)), round(n*outn)))
  return(list(data = data, outliers = labels))
}

########################################################################
## load the methods 
########################################################################

################################################################
## set the working directory 
################################################################
#setwd("") 


source("Methods/POD.R")
source("Methods/FAMDAD.R")
source("Methods/ZDisc_Kmeans.R")
source("Methods/PCA_T2.R")
source("Methods/SECODA_run.R")
source("Methods/ABOD.R")
source("Methods/FRGOD.R")


#############################################################################################################################################################
## experiments 
## n = c(100, 500, 1000)
## p = c(10, 50, 100)

## firstly: only n= 500, and p=100
## cat continuously from 1 to p; 
## nout continuously from 0.01 to 0.5
## rho continuously from  0 to 0.9

#############################################################################################################################################################

##################################
## !!! Methods MIX: create the file data in the working directory!!!
##################################

simulation <- function(parameter, out.type = 'both', n = 500, p = 30, repeats = 100,  methods = c("POD", "FAMDAD", "SECODA_out", "ZDisc", "KMeans", "PCAmix_T2", "ABOD", "MIX")) {
  if (parameter == "nout") {
    outn = seq(0.01, 0.5,  0.02)
    len = length(outn)
    AUC <- data.frame(matrix(NA, nrow = repeats*len, ncol = (length(methods)+2)))
    names(AUC) = c("Sim", "outn", methods)
    AUC$outn = rep(outn, repeats)
    Time <- AUC
    strparam <- paste0("n=",n,", p=",p, ", cat=", round(p/2), ", outn=outn[i], rho=0.5, out.type='",out.type,"', mu=1.5, seed = rep")
    params <- as.list(parse(text=paste0("f(", strparam , ")"))[[1]])[-1]
  } else {
    if (parameter == "cat") {
      cat = seq(1, 27, 2)
      len = length(cat)
      AUC <- data.frame(matrix(NA, nrow = repeats*len, ncol = (length(methods)+2)))
      names(AUC) = c("Sim", "cat", methods)
      AUC$cat = rep(cat, repeats)
      Time <- AUC
      strparam <- paste0("n=",n,", p=",p, ", cat =cat[i], outn=0.1, rho=0.5, out.type='",out.type,"', mu=1.5, seed = rep")
      params <- as.list(parse(text=paste0("f(", strparam , ")"))[[1]])[-1]
    } else {
      if(parameter == "rho") {
        rho = seq(0, 0.9, 0.05)
        len = length(rho)
        AUC <- data.frame(matrix(NA, nrow = repeats*len, ncol = (length(methods)+2)))
        names(AUC) = c("Sim", "rho", methods)
        AUC$rho = rep(rho, repeats)
        Time <- AUC
        strparam <- paste0("n=500, p=30, cat =15, outn=0.1, rho=rho[i], out.type='",out.type,"', mu=1.5, seed = rep")
        params <- as.list(parse(text=paste0("f(", strparam , ")"))[[1]])[-1]
      } else { 
        if (parameter == "mu") {
          mu = seq(0.3, 5, 0.2)
          len = length(mu)
          AUC <- data.frame(matrix(NA, nrow = repeats*len, ncol = (length(methods)+2)))
          names(AUC) = c("Sim", "mu", methods)
          AUC$mu = rep(mu, repeats)
          Time <- AUC
          strparam <- paste0("n=",n,", p=",p, ", cat=", round(p/2), ", outn=0.1, rho=0.5, out.type='",out.type,"', mu=mu[i], seed = rep")
          params <- as.list(parse(text=paste0("f(", strparam , ")"))[[1]])[-1]
        } else {
          if (parameter == "n") {
            n =  seq(100, 1000, 100)
            len = length(n)
            AUC <- data.frame(matrix(NA, nrow = repeats*len, ncol = (length(methods)+2)))
            names(AUC) = c("Sim", "n", methods)
            AUC$n = rep(n, repeats)
            Time <- AUC
            strparam <- paste0("n= n[i], p = ",p, ", cat=", round(p/2), ", outn=0.1, rho=0.5, out.type='",out.type,"', mu=1.5, seed = rep")
            params <- as.list(parse(text=paste0("f(", strparam , ")"))[[1]])[-1]
          } else {
            if (parameter == "p") {
              p =  seq(10, 100, 10)
              len = length(p)
              AUC <- data.frame(matrix(NA, nrow = repeats*len, ncol = (length(methods)+2)))
              names(AUC) = c("Sim", "p", methods)
              AUC$p = rep(p, repeats)
              Time <- AUC
              strparam <- paste0("n=", n,", p=p[i], cat= round(p[i]/2) , outn=0.1, rho=0.5, out.type='",out.type,"', mu=1.5, seed = rep")
              params <- as.list(parse(text=paste0("f(", strparam , ")"))[[1]])[-1]
            }
           else {
             if (parameter == "class") {
               class = 2:12
               len = length(class)
               AUC <- data.frame(matrix(NA, nrow = repeats*len, ncol = (length(methods)+2)))
               names(AUC) = c("Sim", "class", methods)
               AUC$class = rep(class, repeats)
               Time <- AUC
               strparam <- paste0("n=",n,", p=",p, ", cat =round(p/2), outn=0.1, rho=0.5, out.type='",out.type,"', mu=1.5, seed = rep, k = class[i]")
               params <- as.list(parse(text=paste0("f(", strparam , ")"))[[1]])[-1]
             } else {
               cat("Invalid parameter input \n")
             }
          }
        }
      }
    }
  }
 }

  
  for(rep in 1:repeats) {
    for (i in 1:len) {
      cat("Simulation run", rep, "from", repeats, ", num", i, "from", len, "\n")
      p.data = do.call(produce.data, params)
      AUC[(rep-1)*len + i, "Sim"] <- rep
      Time[(rep-1)*len + i, "Sim"] <- rep
      for(method in methods) {
        skip_to_next <- FALSE
        if(method == "ZDisc") {
          tic()
          tryCatch(res <- ZDisc_KMeans(p.data$data, k = 9),
                   error = function(e) { skip_to_next <<- TRUE})
          t <- toc()
          mydata <- data.frame(scores = res$scores, outliers = p.data$outliers)
          AUC[(rep-1)*len + i, method] <- auc(outliers~scores, data=mydata, plot = FALSE)
          Time[(rep-1)*len + i, method] <- unname(t$toc - t$tic)
        } else {
          if (method == "KMeans") {
            tic()
            tryCatch(res <- ZDisc_KMeans(p.data$data,  method = "KMeans", k = 9),
                     error = function(e) { skip_to_next <<- TRUE})
            t <- toc()
            mydata <- data.frame(scores = res$scores, outliers = p.data$outliers)
            AUC[(rep-1)*len + i, method] <- auc(outliers~scores, data=mydata, plot = FALSE)
            Time[(rep-1)*len + i, method] <- unname(t$toc - t$tic)
          } else {
            if (method == "MIX") {
              data.mix <- data.frame(cbind(p.data$data, p.data$outliers))
              names(data.mix) <- c(paste0("A", 1:sum(sapply(p.data$data, is.factor))), paste0("B", 1:sum(sapply(p.data$data, is.numeric))), "class")
              write.csv(sapply(data.mix, as.numeric), paste0("/data/", i, ".csv"), row.names = FALSE)
              } else { 
              if (method == "FRGOD") {
                f <- get(method)
                tryCatch(res <- f(p.data$data),
                         error = function(e) { skip_to_next <<- TRUE})
                mydata <- data.frame(scores = res$scores, outliers = p.data$outliers)
                AUC[(rep-1)*len+ i, method] <- auc(outliers~scores, data=mydata, plot = FALSE)
                Time[(rep-1)*len + i, method] <- res$time
              } else {
                f <- get(method)
                tic()
                tryCatch( res <- f(p.data$data),
                         error = function(e) { skip_to_next <<- TRUE})
                t <- toc()
                if (! skip_to_next) {
                  mydata <- data.frame(scores = res$scores, outliers = p.data$outliers)
                  AUC[(rep-1)*len+ i, method] <- auc(outliers~scores, data=mydata, plot = FALSE)
                  Time[(rep-1)*len + i, method] <- unname(t$toc - t$tic)
                }
              }
            }
          }
        }
      }
    }
    if ("MIX" %in% methods) {
      unlink('C:/Users/Sanja/Desktop/Faks/DataScience/MasterThesis/Code/repo/Methods-/Methods/MIX-master/out.txt')
      system("python C:/Users/Sanja/Desktop/Faks/DataScience/MasterThesis/Code/repo/Methods-/Methods/MIX-master/main.py")
      mix.auc <- read.table('C:/Users/Sanja/Desktop/Faks/DataScience/MasterThesis/Code/repo/Methods-/Methods/MIX-master/out.txt', sep = ",", col.names = c("num", "AUC", "rep", "time"))
      AUC[((rep-1)*len+1):(rep*len), "MIX"] <- mix.auc[order(mix.auc$num), "AUC"]
      as.time <- function(i, times) {
        return(as.numeric(sub("s", "", times[i])))
      }
      Time[((rep-1)*len+1):(rep*len), "MIX"] <-  sapply(1:nrow(mix.auc), as.time, mix.auc[order(mix.auc$num), "time"])
      unlink("C:/Users/Sanja/Desktop/Faks/DataScience/MasterThesis/Code/repo/Methods-/Methods/MIX-master/data/*")
     print(sapply(1:nrow(mix.auc), as.time, mix.auc[order(mix.auc$num), "time"]))
    }
  }

  write.csv(AUC, file = paste0("C:/Users/Sanja/Desktop/Faks/DataScience/MasterThesis/Code/repo/Methods-/Results/AUC_const", parameter,"_rep_", repeats,"_outtype_", out.type, ".csv"), row.names = FALSE)
  write.csv(Time, file = paste0("C:/Users/Sanja/Desktop/Faks/DataScience/MasterThesis/Code/repo/Methods-/Results/Time_const", parameter,"_rep_", repeats,"_outtype_", out.type, ".csv"), row.names = FALSE)
  return(list(AUC = AUC, Time = Time))
}


###################################################################################################################################################
## Experiments
###################################################################################################################################################

simulation("mu",  out.type = "both", repeats = 100)
simulation("mu",  out.type = "cat", repeats = 100)
simulation("mu",  out.type = "cont", repeats = 100)

simulation("rho",  out.type = "both", repeats = 100)
simulation("rho",  out.type = "cat", repeats = 100)
simulation("rho", out.type = "cont", repeats = 100)

simulation("nout",  out.type = "both", repeats = 100)
simulation("nout",  out.type = "cat", repeats = 100)
simulation("nout", out.type = "cont", repeats = 100)

simulation("cat",  out.type = "both", repeats = 100)
simulation("cat",  out.type = "cat", repeats = 100)
simulation("cat",  out.type = "cont", repeats = 100)

simulation("n",  out.type = "both", repeats = 100)
simulation("n",  out.type = "cat", repeats = 100)
simulation("n",  out.type = "cont", repeats = 100)

simulation("p",  out.type = "both", repeats = 100)
simulation("p",  out.type = "cat", repeats = 100)
simulation("p",  out.type = "cont", repeats = 100)


## REAL DATA 
## SANJA PRISERLAC
## 11/05/2022

library(caret)
library(pROC)
library(matlabr)
library(tictoc)

###################################################################################################################################################

real_data <- function(data, outliers, knn = round(0.01*nrow(data)), ntrees = 500,  k = round(log2(nrow(data)+1)), min_sup = 0.5, mc = 2) {
  methods = c("POD",  "ABOD", "SECODA", "ZDisc", "KMeans", "FAMDAD", "PCAmix_T2", "MIX")
  AUC <- rep(NA, length(methods))
  names(AUC) <- methods
  Time = AUC

  tic()
  res <- POD(data, n = knn)
  t <- toc()
  mydata <- data.frame(scores = res$scores, outliers = outliers)
  AUC["POD"] <- auc(outliers~scores, data=mydata, plot = FALSE)
  Time["POD"] <- unname(t$toc - t$tic)
  
  tic()
  res <- FAMDAD(data, ntrees = ntrees)
  t <- toc()
  mydata <- data.frame(scores = res$scores, outliers = outliers)
  AUC["FAMDAD"] <- auc(outliers~scores, data=mydata, plot = FALSE)
  Time["FAMDAD"] <- unname(t$toc - t$tic)
 
  tic()
  res <- SECODA_out(data)
  t <- toc()
  mydata <- data.frame(scores = res$scores, outliers = outliers)
  AUC["SECODA"] <- auc(outliers~scores, data=mydata, plot = FALSE)
  Time["SECODA"] <- unname(t$toc - t$tic)
  
  tic()
  res <- ZDisc_KMeans(data, k = k)
  t <- toc()
  mydata <- data.frame(scores = res$scores, outliers = outliers)
  AUC["ZDisc"] <- auc(outliers~scores, data=mydata, plot = FALSE)
  Time["ZDisc"] <- unname(t$toc - t$tic)
  
  tic()
  res <- ZDisc_KMeans(data, method = "KMeans", k = k)
  t <- toc()
  mydata <- data.frame(scores = res$scores, outliers = outliers)
  AUC["KMeans"] <- auc(outliers~scores, data=mydata, plot = FALSE)
  Time["KMeans"] <- unname(t$toc - t$tic)
  
  tic()
  res <- PCAmix_T2(data)
  t <- toc()
  mydata <- data.frame(scores = res$scores, outliers = outliers)
  AUC["PCAmix_T2"] <- auc(outliers~scores, data=mydata, plot = FALSE)
  Time["PCAmix_T2"] <- unname(t$toc - t$tic)
  
  tic()
  res <- ABOD(data, min_sup = min_sup, mc = mc, k = k)
  t <- toc()
  mydata <- data.frame(scores = res$scores, outliers = outliers)
  AUC["ABOD"] <- auc(outliers~scores, data=mydata, plot = FALSE)
  Time["ABOD"] <- unname(t$toc - t$tic)
  
  data.mix <- cbind(data[, sapply(data, is.factor)], data[, sapply(data, is.numeric)])
  names(data.mix) <- c(paste0("A", 1:sum(sapply(data.mix, is.factor))), paste0("B", 1:sum(sapply(data.mix, is.numeric))))
  data.mix <- cbind(data.mix, class = outliers)
  unlink("C:/Users/Sanja/Desktop/Faks/DataScience/MasterThesis/Code/repo/Methods-/Methods/MIX-master/data/*")
  unlink('C:/Users/Sanja/Desktop/Faks/DataScience/MasterThesis/Code/repo/Methods-/Methods/MIX-master/out.txt')
  write.csv(sapply(data.mix, as.numeric), paste0("C:/Users/Sanja/Desktop/Faks/DataScience/MasterThesis/Code/repo/Methods-/Methods/MIX-master/data/data.csv"), row.names = FALSE)
  setwd("C:/Users/Sanja/anaconda3/envs/MT_MIXnew/")
  system("python C:/Users/Sanja/Desktop/Faks/DataScience/MasterThesis/Code/repo/Methods-/Methods/MIX-master/main.py")
  
  mix.auc <- read.table('C:/Users/Sanja/Desktop/Faks/DataScience/MasterThesis/Code/repo/Methods-/Methods/MIX-master/out.txt', sep = ",", col.names = c("num", "AUC", "rep", "time"))
  AUC["MIX"] <- mix.auc[order(mix.auc$num), "AUC"]
  as.time <- function(i, times) {
    return(as.numeric(sub("s", "", times[i])))
  }
  Time["MIX"] <-  sapply(1:nrow(mix.auc), as.time, mix.auc[order(mix.auc$num), "time"])
  return(list(AUC = AUC, Time = Time))
}


##################################################################################################################################################################
## German Credit
##################################################################################################################################################################
### !!! ADD THE PATH TO DATA 
data.path <- ""

germancredit <- read.table(paste0(data.path,"RealData/GermanCredit/german.data"), sep = " ")
str(germancredit)
sum(is.na(germancredit)) ## no missing values

for (i in 1:ncol(germancredit)) {  if (is.character(germancredit[, i])) { germancredit[, i] <- as.factor(germancredit[, i])}}

sum(sapply(germancredit, is.factor))
## 7 factor variables and 13 factor
sapply(germancredit[, sapply(germancredit, is.numeric)], unique)
## V18, V16 also 

str(germancredit[, sapply(germancredit, is.factor)])
summary(as.factor(germancredit[, 21] - 1))/nrow(germancredit)
#   0   1 
# 700 300 

res.pod.def <- real_data(germancredit[, -21], outliers = as.factor(germancredit[, 21] - 1), k = 2, mc = 10, min_sup = 0.5)

tic()
pod.abod <- ABOD(germancredit[, -21], min_sup = 0.5, mc = 10)
toc()
mydata <- data.frame(scores = pod.abod$scores, outliers = as.factor(germancredit[, 21] - 1))
auc(outliers~scores, data=mydata, plot = FALSE)


res.pod.opt <- real_data(germancredit[, -21], outliers = as.factor(germancredit[, 21] - 1), k = 2, mc = 10)

##################################################################################################################################################################
## Adult
##################################################################################################################################################################

adult <- read.table(paste0(data.path, "RealData/Adult/adult.data"), sep = ",")
names(adult) <- c("age", "workclass", "fnlwgt", "education", "education-num", "marital-status", "occupation", "relationship", "race", "sex", "capital-gain", "capital-loss", "hours-per-week", "native-country", "salary")
str(adult[, sapply(adult[, -15], is.factor)])
for (i in which(sapply(adult, is.character)==TRUE)) {  adult[, i] <- as.factor(adult[, i])}
levels(adult[, "marital-status"])

sum(sapply(adult[, -15], is.factor))  # 8 factors
str(adult[, sapply(adult[, -15], is.factor)])
sum(sapply(adult[, -15], is.numeric)) # 6 numeric


sum(is.na(adult))
summary(as.factor(adult$salary))/nrow(adult)
# <=50K      >50K 
# 0.7591904 0.2408096 
summary(as.factor(as.numeric(adult$`hours-per-week`>80)))/nrow(adult)

res.adult.def <- real_data(adult[, -15], outliers = (as.numeric(adult[, 15])- 1), knn = 100, k = 8, mc = 10)


##################################################################################################################################################################
## Heart
##################################################################################################################################################################

heart <- read.csv(paste0(data.path, "RealData/Cleveland/processed.cleveland.data"), header = FALSE)
## 6 missing values
## some missing values excluded!
sum(heart[, 14]==0)/nrow(heart)  # 0.5412541
sum(heart[, 14]==0)/sum(heart[, 14]%in% c(0, 1))  # 0.7488584
nrow(heart) #303
heart <- heart[which(heart[, 14] %in% c(0, 1)), ]
heart <- heart[-c(which(heart$V13 == "?"), which(heart$V12 == "?")), ]
nrow(heart) #214

for (i in c(2, 3, 6, 7, 9, 11, 12, 13)) {
  heart[, i] <- as.factor(heart[, i])
}

str(heart[, sapply(heart, is.factor)])

summary(as.factor(heart[, 14]))/nrow(heart)

data <- heart[, -14]
outliers <- ifelse(heart[, 14]== 0, 0, 1)

res.heart <- real_data(heart[, -14], outliers = ifelse(heart[, 14]== 0, 0, 1), min_sup = 0.1, k = 7, knn = 3)


##################################################################################################################################################################
## CMC
##################################################################################################################################################################

cmc <- read.csv(paste0(data.path,"RealData/Contraceptive/cmc.data"), header = FALSE)
str(cmc)
sum(is.na(cmc))
cmc[, c(1, 4)] <- sapply(cmc[, c(1, 4)], as.numeric)
for (i in which(sapply(cmc, is.integer)==TRUE)) {  cmc[, i] <- as.factor(cmc[, i])}
str(cmc)

data = cmc[, -10]
outliers <- ifelse(cmc[, 10] == 2, 1, 0)
sum(outliers)/length(outliers)

res.cmc <- real_data(cmc[, -10], outliers = ifelse(cmc[, 10] == 2, 1, 0), k = 8, mc = 10)

##################################################################################################################################################################
## Thoracic Surgery
##################################################################################################################################################################

thoracic <- read.csv(paste0(data.path,"RealData/ThoraricSurgery.arff"), header=FALSE, comment.char = "@")
str(thoracic)

for (i in which(sapply(thoracic, is.logical)==TRUE)) {  thoracic[, i] <- as.factor(thoracic[, i])}
for (i in which(sapply(thoracic, is.character)==TRUE)) {  thoracic[, i] <- as.factor(thoracic[, i])}
str(thoracic)
sum(complete.cases(thoracic)) # no missing values


summary(thoracic[, 17])/nrow(thoracic)
data = thoracic[, -17]
outliers = as.numeric(thoracic[, 17]) -1

res.surgery <- real_data(thoracic[, -17], outlier = (as.numeric(thoracic[, 17]) -1), k = 8)

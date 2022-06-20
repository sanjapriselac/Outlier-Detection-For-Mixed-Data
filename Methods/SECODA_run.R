## Sanja Priselac
## last edit on 31/03/2022
## Method SECODA for run 

## added by Sanja Priselac on 30/03/2022
SECODA_out <- function(datset, cutoff = 0.5, BinningMethod="EW", HighDimMode="IN", MinimumNumberOfIterations=3, MaximumNumberOfIterations=99999, StartHeuristicsAfterIteration=10, FractionOfCasesToRetain=0.2, TestMode=FALSE) {
  source("C:/Users/Sanja/Desktop/Faks/DataScience/MasterThesis/Code/repo/Methods-/Methods/SECODA.R")
  
  secoda <- SECODA(datset, BinningMethod, HighDimMode, MinimumNumberOfIterations, MaximumNumberOfIterations, StartHeuristicsAfterIteration, FractionOfCasesToRetain, TestMode)
  secoda.scores <- 1 - (secoda$AveAnoScore -min(secoda$AveAnoScore))/(max(secoda$AveAnoScore) - min(secoda$AveAnoScore))
  outliers <- rep(0, nrow(datset))
  outliers[which(secoda.scores > cutoff)] <- 1
  return(list(scores = secoda.scores, outliers = outliers))
}

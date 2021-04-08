#' Turn Distance Matrix into Longform Table
#' @param Things character vector: strings
#' @param Matches character vector: matches, comma separated
#' @param allabsnum numeric: total number of abstracts in analysis
#' @export
# Written by John M. A. Wojahn July 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU Affero General Public License v. 3
# Funded by EPSCoR GEM3 at Boise State University

InternalPairwiseDistanceInferreR <- function(Things, Matches, allabsnum)
{
  Things <- as.vector(Things)
  Matches <- as.vector(Matches)
  termsmatches <- data.frame(Things,Matches)
  PWDM <- as.data.frame(matrix(nrow = length(Things), ncol = length(Things)))
  rownames(PWDM) <- Things
  colnames(PWDM) <- Things
  print("...Calculating Pairwise Distances")
  pb <- txtProgressBar(min = 1, max = nrow(PWDM), style = 3)
  #perform distance calculation
  for(i in 1:nrow(PWDM))
  {
    setTxtProgressBar(pb, i)
    for(j in 1:ncol(PWDM))
    {
      #print(paste0("Row ",i," Column ",j))
      FirstNomen <- rownames(PWDM)[i]
      SecondNomen <- colnames(PWDM)[j]
      smatches1 <- termsmatches[which(termsmatches[,1] == FirstNomen),2]
      smatches2 <- termsmatches[which(termsmatches[,1] == SecondNomen),2]
      smatches1 <- as.vector(unlist(strsplit(as.character(smatches1), split = ",")))
      smatches2 <- as.vector(unlist(strsplit(as.character(smatches2), split = ",")))
      smatches1num <- length(smatches1)
      smatches2num <- length(smatches2)
      sallnum <- length(unique(c(smatches1,smatches2)))
      onetwo <- smatches1[which(smatches1 %in% smatches2)]
      twoone <- smatches2[which(smatches2 %in% smatches1)]
      inz <- unique(c(onetwo,twoone))
      inznum <- length(inz)
      vallll <- as.numeric(inznum/sallnum) * sallnum/allabsnum
      Threshold <- 0
      if(vallll >= Threshold)
      {
        PWDM[i,j] <- vallll
      }else{
        PWDM[i,j] <- 0
      }
    }
  }
  #Longed <- LongformTableMakeR(PWDM, decreasing = F)
  #Longed <- as.data.frame(Longed[-which(Longed$Distance == 0),])
  #LongedSummary <- as.data.frame(as.matrix(summary(Longed$Distance)))
  #LongMean <- as.numeric(LongedSummary[4,])
  print("...Normalizing Distances")
  #pb <- txtProgressBar(min = 1, max = nrow(PWDM), style = 3)
  #perform distance calculation
 # for(i in 1:nrow(PWDM))
  #{
   # setTxtProgressBar(pb, i)
   # for(j in 1:ncol(PWDM))
   # {
    #  if(as.numeric(PWDM[i,j]) < .50)
    #  {
      #  PWDM[i,j] <- 0
    #  }
   # }
#  }
  print("...Zeroing Selfs")
  for(i in 1:nrow(PWDM))
  {
    PWDM[i,i] <- 0
  }
  print("...Removing empty nodes")
  RowSumz <- c(1:nrow(PWDM))
  ColSumz <- c(1:ncol(PWDM))
  for(i in 1:nrow(PWDM))
  {
    RowSumz[i] <- sum(PWDM[i,])
  }
  for(i in 1:ncol(PWDM))
  {
    ColSumz[i] <- sum(PWDM[,i])
  }
  BadRows <- which(as.numeric(RowSumz) == 0)
  BadCols <- which(as.numeric(ColSumz) == 0)
  if(length(BadRows) > 0 && length(BadCols) > 0)
  {
    SmallOut <- as.data.frame(PWDM[-BadRows,-BadCols])
  }else if(length(BadRows) > 0 && length(BadCols) == 0){
    SmallOut <- as.data.frame(PWDM[-BadRows,])
  }else if(length(BadRows) == 0 && length(BadCols) > 0){
    SmallOut <- as.data.frame(PWDM[,-BadCols])
  }else{
    SmallOut <- PWDM
  }
  print(sprintf("Returning Matrix with %s Rows and %s Columns.",nrow(SmallOut),ncol(SmallOut)))
  return(SmallOut)
}

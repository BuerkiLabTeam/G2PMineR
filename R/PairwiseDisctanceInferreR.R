#' Finds weighted pairwise distances between two levels using shared abstracts
#' @param terms1 a vector of unique strings
#' @param matches1 a vector of unique identifiers for abstracts collated by commas for each term
#' @param terms2 a vector of unique strings
#' @param matches2 a vector of unique identifiers for abstracts collated by commas for each term
#' @param allabsnum total number of abstracts in analysis
#' @export
# Written by John M. A. Wojahn July 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU Affero General Public License v. 3
# Funded by EPSCoR GEM3 at Boise State University

PairwiseDistanceInferreR <- function(terms1, matches1, terms2, matches2, allabsnum)
{
  #Process Input No 1
  terms1 <- as.vector(terms1)
  matches1 <- as.vector(matches1)
  termsmatches1 <- data.frame(terms1, matches1)
  termsmatches1 <- unique(termsmatches1)
  #Process Input No 2
  terms2 <- as.vector(terms2)
  matches2 <- as.vector(matches2)
  termsmatches2 <- data.frame(terms2, matches2)
  termsmatches2 <- unique(termsmatches2)
  nummatch1 <- (nchar(as.character(matches1)) - nchar(gsub(",","",as.character(matches1)))) + 1
  nummatch2 <- (nchar(as.character(matches2)) - nchar(gsub(",","",as.character(matches2)))) + 1
  allabsnum <- as.numeric(allabsnum)
  #Find lengths of both inputs
  NumRows1 <- nrow(termsmatches1)
  NumRows2 <- nrow(termsmatches2)

  #make distance matrix
  DistMatDF <- as.data.frame(matrix(nrow=NumRows1,ncol=NumRows2))
  rownames(DistMatDF) <- termsmatches1[,1]
  colnames(DistMatDF) <- termsmatches2[,1]
  pb <- txtProgressBar(min = 1, max = nrow(DistMatDF), style = 3)
  print("...Calculating Pairwise Distance")
  #perform distance calculation
  for(i in 1:nrow(DistMatDF))
  {
    setTxtProgressBar(pb, i)
    for(j in 1:ncol(DistMatDF))
    {
      FirstNomen <- rownames(DistMatDF)[i]
      SecondNomen <- colnames(DistMatDF)[j]
      smatches1 <- termsmatches1[which(termsmatches1[,1] == FirstNomen),2]
      smatches2 <- termsmatches2[which(termsmatches2[,1] == SecondNomen),2]
      smatches1 <- as.vector(unlist(strsplit(as.character(smatches1), split = ",")))
      smatches2 <- as.vector(unlist(strsplit(as.character(smatches2), split = ",")))
      smatches1num <- length(smatches1)
      smatches2num <- length(smatches2)
      sallnum <- length(unique(c(smatches1,smatches2)))
      onetwo <- smatches1[which(smatches1 %in% smatches2)]
      twoone <- smatches2[which(smatches2 %in% smatches1)]
      inz <- unique(c(onetwo,twoone))
      inznum <- length(inz)
      inzprop <- inznum/sallnum
      weight <- sallnum/allabsnum
      weightedvalue <- as.numeric(inzprop*weight)
      DistMatDF[i,j] <- weightedvalue
    }
  }
  PWDM <- as.data.frame(DistMatDF)
  print("...Making Longform Table")
  DistMat <- as.data.frame(DistMatDF)
  numcells <- as.numeric(nrow(DistMat) * ncol(DistMat))
  LongOut <- as.data.frame(matrix(nrow = numcells, ncol = 3))
  colnames(LongOut) <- c("Term1", "Term2", "Distance")
  cellcount <- 1
  pb <- txtProgressBar(min = 1, max = nrow(DistMat), style = 3)
  start <- 0
  for(i in 1:nrow(DistMat))
  {
    setTxtProgressBar(pb, i)
    inz <- as.numeric(as.vector(DistMat[i,]))
    inznomencols <- as.character(colnames(DistMat))
    insnomenonerow <- rep_len(as.character(rownames(DistMat)[i]),length(inznomencols))
    #tobind <- data.frame(insnomenonerow,inznomencols,inz)
    if(i == 1)
    {
      start <- 1
      end <- length(inznomencols)
    }else{
      start <- start  + length(inznomencols)
      end <- end  + length(inznomencols)
    }
    LongOut[start:end,1] <- insnomenonerow
    LongOut[start:end,2] <- inznomencols
    LongOut[start:end,3] <- inz
  }
  decreasing <- T
  if(decreasing == T)
  {
    LongOut <- as.data.frame(LongOut[order(LongOut$Distance, decreasing = T),])
  }else{
    LongOut <- as.data.frame(LongOut)
  }
  Longed <- LongOut
  Longed <- as.data.frame(Longed[-which(Longed$Distance == 0),])
  LongedSummary <- as.data.frame(as.matrix(summary(Longed$Distance)))
  LongMean <- as.numeric(LongedSummary[4,])

  #print("...Normalizing Distances")
  #pb <- txtProgressBar(min = 1, max = nrow(PWDM), style = 3)
  #perform distance calculation
  #for(i in 1:nrow(PWDM))
  #{
  #  setTxtProgressBar(pb, i)
  #  for(j in 1:ncol(PWDM))
  #  {
  #    if(as.numeric(PWDM[i,j]) < LongMean)
  #    {
  #      PWDM[i,j] <- 0
  #    }
  #  }
  #}
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
  if(length(BadRows) != 0)
  {
    if(length(BadCols) != 0)
    {
      SmallOut <- as.data.frame(PWDM[-BadRows,-BadCols])
    }else{
      SmallOut <- as.data.frame(PWDM[-BadRows,])
    }
  }else{
    if(length(BadCols) != 0)
    {
      SmallOut <- as.data.frame(PWDM[,-BadCols])
    }else{
      SmallOut <- as.data.frame(PWDM)
    }
  }

  print(sprintf("Returning Matrix with %s Rows and %s Columns.",nrow(SmallOut),ncol(SmallOut)))
  return(SmallOut)
}

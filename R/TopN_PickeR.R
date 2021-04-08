#' Picks n Most Similar Pairs
#' @param matrix data.frame or matrix: a distance matrix dataframe
#' @param n numeric: how many tops t take a numeric
#' @param decreasing boolean: order of output
#' @export
# Written by John M. A. Wojahn July 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU Affero General Public License v. 3
# Funded by EPSCoR GEM3 at Boise State University

TopN_PickeR <- function(matrix, n, decreasing)
{
  DistMat <- as.data.frame(matrix)
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
  if(decreasing == T)
  {
    LongOut <- as.data.frame(LongOut[order(LongOut$Distance, decreasing = T),])
  }else{
    LongOut <- as.data.frame(LongOut)
  }
  matrixLF <- LongOut
  TopMat <- matrixLF[1:n,]
  inz <- c(1:n)
  colnomen <- unique(TopMat[,1])
  rownomen <- unique(TopMat[,2])
  colnomen <- colnomen[!is.na(colnomen)]
  rownomen <- rownomen[!is.na(rownomen)]
  OutMat <- as.data.frame(matrix(nrow=length(rownomen),ncol=length(colnomen)))
  colnames(OutMat) <- colnomen
  rownames(OutMat) <- rownomen
  OutMat[,] <- 0
  pb <- txtProgressBar(min = 1, max = nrow(TopMat), style = 3)
  for(i in 1:nrow(TopMat))
  {
    #print(paste0(i," of ",nrow(TopMat)))
    setTxtProgressBar(pb, i)
    OutMat[which(rownames(OutMat) == as.character(TopMat[i,2])),
           which(colnames(OutMat) == as.character(TopMat[i,1]))] <-
                                                        as.numeric(TopMat[i,3])

  }
  return(OutMat)
}

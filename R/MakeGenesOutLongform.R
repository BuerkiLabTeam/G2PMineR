#' Makes GenesOut Longform
#' @param GenesOut data.frame: output of SynonymReplaceR
#' @export
# Written by John M. A. Wojahn August 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University

MakeGenesOutLongform <- function(GenesOut)
{
  GenesOut <- as.data.frame(GenesOut)
  LFGO <- as.data.frame(matrix(nrow=0, ncol = ncol(GenesOut)))
  colnames(LFGO) <- colnames(GenesOut)
  print("Processing Abstract Matches")
  pb <- txtProgressBar(min = 1, max = nrow(GenesOut), style = 3)
  for(i in 1:nrow(GenesOut))
  {
    setTxtProgressBar(pb, i)
    if(grepl(",",as.character(GenesOut[i,3])))
    {
      splitted <- as.vector(unlist(strsplit(as.character(GenesOut[i,3]), split = ",")))
      tobind <- as.data.frame(matrix(nrow=length(splitted), ncol = ncol(GenesOut)))
      tobind[1:length(splitted),1] <- as.character(GenesOut[i,1])
      tobind[1:length(splitted),2] <- as.character(GenesOut[i,2])
      tobind[1:length(splitted),3] <- splitted
      tobind[1:length(splitted),4] <- as.character(GenesOut[i,4])
      tobind[1:length(splitted),5] <- as.character(GenesOut[i,5])
      tobind[1:length(splitted),6] <- as.character(GenesOut[i,6])
      colnames(tobind) <- colnames(GenesOut)
      LFGO <- as.data.frame(rbind(LFGO, tobind))
    }else{
      LFGO <- as.data.frame(rbind(LFGO, GenesOut[i,]))
    }
  }
  print("Processing Ontology Matches")
  LFGOO <- as.data.frame(matrix(nrow=0, ncol = ncol(LFGO)))
  colnames(LFGOO) <- colnames(LFGO)
  pb <- txtProgressBar(min = 1, max = nrow(LFGO), style = 3)
  for(i in 1:nrow(LFGO))
  {
    setTxtProgressBar(pb, i)
    if(grepl("; ",as.character(LFGO[i,6])))
    {
      splitted <- as.vector(unlist(strsplit(as.character(LFGO[i,6]), split = "; ")))
      tobind <- as.data.frame(matrix(nrow=length(splitted), ncol = ncol(LFGO)))
      tobind[1:length(splitted),1] <- as.character(LFGO[i,1])
      tobind[1:length(splitted),2] <- as.character(LFGO[i,2])
      tobind[1:length(splitted),3] <- as.character(LFGO[i,3])
      tobind[1:length(splitted),4] <- as.character(LFGO[i,4])
      tobind[1:length(splitted),5] <- as.character(LFGO[i,5])
      tobind[1:length(splitted),6] <- splitted
      colnames(tobind) <- colnames(LFGO)
      LFGOO <- as.data.frame(rbind(LFGOO, tobind))
    }else{
      LFGOO <- as.data.frame(rbind(LFGOO, LFGO[i,]))
    }
  }
  return(LFGOO)
}

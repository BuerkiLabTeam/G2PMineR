#' Makes AbstractsSpp Longform
#' @param AbsPhen data.frame: output of SpeciesLookeR
#' @export
# Written by John M. A. Wojahn August 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University

MakeAbsPhenLongform <- function(AbsPhen)
{
  AbsPhen <- as.data.frame(AbsPhen)
  LFAS <- as.data.frame(matrix(nrow=0, ncol = ncol(AbsPhen)))
  colnames(LFAS) <- colnames(AbsPhen)
  print("Processing Abstract Matches")
  pb <- txtProgressBar(min = 1, max = nrow(AbsPhen), style = 3)
  for(i in 1:nrow(AbsPhen))
  {
    setTxtProgressBar(pb, i)
    if(grepl(",",as.character(AbsPhen[i,6])))
    {
      splitted <- as.vector(unlist(strsplit(as.character(AbsPhen[i,6]), split = ",")))
      tobind <- as.data.frame(matrix(nrow=length(splitted), ncol = ncol(AbsPhen)))
      tobind[1:length(splitted),1] <- as.character(AbsPhen[i,1])
      tobind[1:length(splitted),2] <- as.character(AbsPhen[i,2])
      tobind[1:length(splitted),3] <- as.character(AbsPhen[i,3])
      tobind[1:length(splitted),4] <- as.character(AbsPhen[i,4])
      tobind[1:length(splitted),5] <- as.character(AbsPhen[i,5])
      tobind[1:length(splitted),6] <- splitted
      colnames(tobind) <- colnames(AbsPhen)
      LFAS <- as.data.frame(rbind(LFAS, tobind))
    }else{
      LFAS <- as.data.frame(rbind(LFAS, AbsPhen[i,]))
    }
  }
  return(LFAS)
}

#' Makes AbstractsSpp Longform
#' @param AbstractsSpp data.frame: output of SpeciesLookeR
#' @export
# Written by John M. A. Wojahn August 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University

MakeAbstractsSppLongform <- function(AbstractsSpp)
{
  AbstractsSpp <- as.data.frame(AbstractsSpp)
  LFAS <- as.data.frame(matrix(nrow=0, ncol = ncol(AbstractsSpp)))
  colnames(LFAS) <- colnames(AbstractsSpp)
  print("Processing Abstract Matches")
  pb <- txtProgressBar(min = 1, max = nrow(AbstractsSpp), style = 3)
  for(i in 1:nrow(AbstractsSpp))
  {
    setTxtProgressBar(pb, i)
    if(grepl(",",as.character(AbstractsSpp[i,3])))
    {
      splitted <- as.vector(unlist(strsplit(as.character(AbstractsSpp[i,3]), split = ",")))
      tobind <- as.data.frame(matrix(nrow=length(splitted), ncol = ncol(AbstractsSpp)))
      tobind[1:length(splitted),1] <- as.character(AbstractsSpp[i,1])
      tobind[1:length(splitted),2] <- as.character(AbstractsSpp[i,2])
      tobind[1:length(splitted),3] <- splitted
      colnames(tobind) <- colnames(AbstractsSpp)
      LFAS <- as.data.frame(rbind(LFAS, tobind))
    }else{
      LFAS <- as.data.frame(rbind(LFAS, AbstractsSpp[i,]))
    }
  }
  return(LFAS)
}

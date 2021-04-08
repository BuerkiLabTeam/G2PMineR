#' Calculates Proportion of Abstracts from which Results Were Obtained
#' @param Result data.frame: either genes, phenotypes, or species results
#' @param IDs character or numeric vector: unique abstract IDs
#' @export
# Written by John M. A. Wojahn June 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU Affero General Public License v. 3
# Funded by EPSCoR GEM3 at Boise State University

AbstractsProportionCalculator <- function(Result, IDs)
{
  options(warn=-1)
  if(exists("prop"))
  {
    rm(prop)
  }
  Result <- as.data.frame(Result)
  if("Matches" %in% colnames(Result))
  {
    sprintf("using matches")
    IDs <- as.character(IDs)
    IDsInSpp <- as.vector(matrix(nrow=length(IDs),ncol=1))
    pb <- txtProgressBar(min = 1, max = length(IDsInSpp), style = 3)
    for(i in 1:length(IDsInSpp))
    {
      setTxtProgressBar(pb, i)
      tmp <- grepl(IDs[i],Result$Matches)
      if(T %in% tmp)
      {
        IDsInSpp[i] <- T
      }else{
        IDsInSpp[i] <- F
      }
    }
    outab <- as.data.frame(table(IDsInSpp))
  }else{
    sprintf("using absmatches")
    IDs <- as.character(IDs)
    IDsInSpp <- as.vector(matrix(nrow=length(IDs),ncol=1))
    pb <- txtProgressBar(min = 1, max = length(IDsInSpp), style = 3)
    for(i in 1:length(IDsInSpp))
    {
      setTxtProgressBar(pb, i)
      tmp <- grepl(IDs[i],Result$AbsMatches)
      if(T %in% tmp)
      {
        IDsInSpp[i] <- T
      }else{
        IDsInSpp[i] <- F
      }
    }
    outab <- as.data.frame(table(IDsInSpp))
  }
  prop <- outab[which(outab[,1] == T),2]/sum(outab[,2])
  options(warn=0)
  return(prop)
}

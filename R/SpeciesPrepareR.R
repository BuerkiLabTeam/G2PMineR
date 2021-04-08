#' Checks If Species is in String
#' @param AbstractSpecies a vector of strings
#' @export
# Written by John M. A. Wojahn June 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University

SpeciesPrepareR <- function(AbstractSpecies)
{
  UniqueSpp <- unique(AbstractSpecies$Species)
  UniqueSpp <- UniqueSpp[-which(UniqueSpp == "Arabidopsis spp.")]
  UniqueSpp <- UniqueSpp[order(UniqueSpp)]
  AbstractsSpeciesCollapsed <- as.data.frame(matrix(nrow=length(UniqueSpp)))
  AbstractsSpeciesCollapsed[,1] <- UniqueSpp
  pb <- txtProgressBar(min = 1, max = nrow(AbstractsSpeciesCollapsed), style = 3)
  for(i in 1:nrow(AbstractsSpeciesCollapsed))
  {
    setTxtProgressBar(pb, i)
    if(AbstractsSpeciesCollapsed[i,1] == "Arabidopsis thaliana")
    {
      tmp <- unique(as.vector(AbstractSpecies[which(AbstractSpecies$Species %in% c("Arabidopsis thaliana","Arabidopsis spp.")),3]))
      tmp <- paste(tmp, collapse = ",")
      AbstractsSpeciesCollapsed[i,2] <- tmp
    }else{
      tmp <- unique(as.vector(AbstractSpecies[which(AbstractSpecies$Species == as.character(AbstractsSpeciesCollapsed[i,1])),3]))
      tmp <- paste(tmp, collapse = ",")
      AbstractsSpeciesCollapsed[i,2] <- tmp
    }
  }
  colnames(AbstractsSpeciesCollapsed) <- c("Species","Matches")
  return(AbstractsSpeciesCollapsed)
}


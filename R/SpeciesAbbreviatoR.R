#' Gets Abstracts IDs From Results Data Using Search
#' @param AbstractsSpp data.frame: output of SpeciesLookeR
#' @export
# Written by John M. A. Wojahn August 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University

SpeciesAbbreviatoR <- function(AbstractsSpp)
{
  AbstractsSpp <- as.data.frame(AbstractsSpp) #rightclass
  Spp <- as.vector(unique(AbstractsSpp$Species)) #unique and vectorize
  Spp <- Spp[grepl(" ",Spp)] #remove just genera
  splitted <- strsplit(as.character(Spp), split=" ")
  Abbreviations <- c(1:length(Spp))
  Abbreviations[] <- NA
  pb <- txtProgressBar(min = 1, max = length(Abbreviations), style = 3)
  for(i in 1:length(splitted))
  {
    setTxtProgressBar(pb, i)
    this <- splitted[[i]]
    this <- as.character(this)
    Abbreviations[i] <- paste(substring(this,1,1), collapse = "")
  }
  mandatories <- c("At","Gm","Zm","Bd","Nb","Pt","Nt")
  Abbreviations <- as.vector(c(mandatories,Abbreviations))
  Abbreviations <- unique(Abbreviations)
  return(Abbreviations)
}

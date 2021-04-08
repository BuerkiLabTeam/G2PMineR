#' Cleans Phenotypes Library
#' @param MoBotTerms character vector: all phenotypes as strings
#' @export
# Written by John M. A. Wojahn June 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU Affero General Public License v. 3
# Funded by EPSCoR GEM3 at Boise State University

CleanPhenotypesLibrary <- function(MoBotTerms)
{
  print("CLEANING PHENOTYPES LIST")
  if(T %in% grepl("\t",MoBotTerms))
  {
    MoBotTerms <- gsub("\t","",MoBotTerms)
  }
  MoBotTerms <- as.data.frame(MoBotTerms)
  MoBotTerms <- as.vector(MoBotTerms[,1])
  MoBotTerms <- MoBotTerms[-grep("(prefix)", MoBotTerms)]
  MoBotTerms <- gsub(" =.*","",MoBotTerms)
  MoBotTerms <- gsub(" \\(.*","",MoBotTerms)
  MoBotTerms <- as.vector(unlist(strsplit(as.character(MoBotTerms),split = ",")))
  MoBotTerms <- as.vector(unlist(strsplit(as.character(MoBotTerms),split = ";")))
  MoBotTerms <- tolower(MoBotTerms)
  MoBotTerms <- MoBotTerms[!(MoBotTerms == ".")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == ",")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == ":")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == ";")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == " ")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "")]
  MoBotTerms <- MoBotTerms[!is.na(MoBotTerms)]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "a")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "b")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "c")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "d")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "e")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "f")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "g")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "h")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "i")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "j")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "k")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "l")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "m")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "n")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "o")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "p")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "q")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "r")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "s")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "t")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "u")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "v")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "w")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "x")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "y")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "z")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "1")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == " ii")]
  MoBotTerms <- MoBotTerms[!(MoBotTerms == "a | b | c | d | e | f | g | h | i | j | k | l | m | n | o | p | q | r | s | t | u | v | w | x |z")]
  return(MoBotTerms)
}

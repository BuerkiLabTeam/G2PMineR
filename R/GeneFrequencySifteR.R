#' Siftes Genes by Gene Frequency
#' @param GenesOut matrix: swissprot genes abstracts search results
#' @param IDs vector: unique identifiers
#' @export
# Written by John M. A. Wojahn September 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by  EPSCoR GEM3 at Boise State University


GeneFrequencySifteR <- function(GenesOut, IDs)
{
  if(class(GenesOut) != "data.frame")
  {
    stop("GenesOut must be a data frame")
  }
  if(class(IDs) != "character" && class(IDs) != "integer" && class(IDs) != "numeric")
  {
    stop("IDs must be a character vector")
  }
  print("...Preprocessing GenesOut")
  GenesOutShorted <- GenesOut
  nummatches <- nchar(gsub("[[:digit:]]","",GenesOutShorted$Matches)) + 1
  nummatches[which(GenesOutShorted$Matches == "No")] <- 0
  totalnum <- length(IDs)
  print("Here are the summary statistics for your matches:")
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print(summary(nummatches))
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  print ("NOTE: USER MUST NOW ENTER DECISION IN PROMPT BELOW")
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  threshold <- readline(prompt = "PLEASE TYPE A THRESHOLD MINIMUM NUMBER OF MATCHES HERE: ")
  GenesOutShorted <- data.frame(GenesOutShorted,nummatches)
  GenesOutShorted <- as.data.frame(GenesOutShorted[order(GenesOutShorted$nummatches, decreasing = T),])
  GenesOutShorted <- as.data.frame(GenesOutShorted[which(as.numeric(GenesOutShorted$nummatches) > as.numeric(threshold)),])
  return(GenesOutShorted)
}

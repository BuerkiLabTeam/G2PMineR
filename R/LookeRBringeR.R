#' Gets Abstracts IDs From Results Data Using Search
#' @param Terms character or numeric vector: unique abstract IDs
#' @param SearchingData character or numeric vector: unique abstract IDs
#' @param IDs character or numeric vector: unique abstract IDs
#' @param Type charcater: either S, G, or P (Species, Genes, or Species Analysis)
#' @export
# Written by John M. A. Wojahn August 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University

LookeRBringeR <- function(Terms, SearchingData, IDs, Type)
{
  #check classes
  if(class(Terms) != "character")
  {
    stop("ERROR: Terms must be a character vector!")
  }
  if(class(IDs) != "character" && class(IDs) != "numeric")
  {
    stop("ERROR: IDs must be a character or numeric vector!")
  }
  if(class(SearchingData) != "data.frame")
  {
    stop("ERROR: SearchingData must be a data.frame!")
  }
  if(Type == "S" || Type == "s")
  {
    foundlings <- as.vector(unlist(strsplit(as.character(SearchingData[which(SearchingData$Species %in% Terms),3]), split = ",")))
  }else if(Type == "G" || Type == "g"){
    foundlings <- as.vector(unlist(strsplit(as.character(SearchingData[which(SearchingData$Gene %in% Terms),3]), split = ",")))
  }else if(Type == "P" || Type == "p"){
    foundlings <- as.vector(unlist(strsplit(as.character(SearchingData[which(SearchingData$PhenoWord %in% Terms),6]), split = ",")))
  }else{
    stop("ERROR: 'Type' must be either 'S', 'G', or 'P'!")
  }
  return(foundlings)
}

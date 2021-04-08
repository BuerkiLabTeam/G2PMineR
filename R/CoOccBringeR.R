#' Gets Abstracts IDs From Results Data Using Search Term
#' @param Terms1 character or numeric vector: unique abstract IDs
#' @param Terms2 character or numeric vector: unique abstract IDs
#' @param SearchingData1 character or numeric vector: unique abstract IDs
#' @param SearchingData2 character or numeric vector: unique abstract IDs
#' @param IDs character or numeric vector: unique abstract IDs
#' @param Type charcater: either S, G, or P (Species, Genes, or Species Analysis)
#' @export
# Written by John M. A. Wojahn August 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University

CoOcBringeR <- function(Terms1, Terms2, SearchingData1, SearchingData2, IDs, Type)
{
 #check classes
  if(class(Terms1) != "character")
  {
    stop("ERROR: Terms must be a character vector!")
  }
  if(class(Terms2) != "character")
  {
    stop("ERROR: Terms must be a character vector!")
  }
  if(class(IDs) != "character" && class(IDs) != "numeric")
  {
    stop("ERROR: IDs must be a character or numeric vector!")
  }
  if(class(SearchingData1) != "data.frame")
  {
    stop("ERROR: SearchingData must be a data.frame!")
  }
  if(class(SearchingData2) != "data.frame")
  {
    stop("ERROR: SearchingData must be a data.frame!")
  }
  if(Type == "SG" || Type == "sg")
  {
    foundlings1 <- as.vector(unlist(strsplit(as.character(SearchingData1[which(SearchingData1$Species %in% Terms1),3]), split = ",")))
    foundlings2 <- as.vector(unlist(strsplit(as.character(SearchingData2[which(SearchingData2$Gene %in% Terms2),3]), split = ",")))
    foundlings <- foundlings1[foundlings1 %in% foundlings2]
  }else if( Type == "GS" || Type == "gs"){
    foundlings1 <- as.vector(unlist(strsplit(as.character(SearchingData1[which(SearchingData1$Species %in% Terms1),3]), split = ",")))
    foundlings2 <- as.vector(unlist(strsplit(as.character(SearchingData2[which(SearchingData2$Gene %in% Terms2),3]), split = ",")))
    foundlings <-  foundlings1[foundlings1 %in% foundlings2]
  }else if(Type == "GP" || Type == "gp"){
    foundlings1 <- as.vector(unlist(strsplit(as.character(SearchingData1[which(SearchingData1$Gene %in% Terms1),3]), split = ",")))
    foundlings2 <- as.vector(unlist(strsplit(as.character(SearchingData2[which(SearchingData2$PhenoWord %in% Terms2),6]), split = ",")))
    foundlings <-  foundlings1[foundlings1 %in% foundlings2]
  }else if(Type == "PG" || Type == "pg"){
    foundlings1 <- as.vector(unlist(strsplit(as.character(SearchingData1[which(SearchingData1$PhenoWord %in% Terms1),6]), split = ",")))
    foundlings2 <- as.vector(unlist(strsplit(as.character(SearchingData2[which(SearchingData2$Gene %in% Terms2),3]), split = ",")))
    foundlings <-  foundlings1[foundlings1 %in% foundlings2]
  }else if(Type == "SP" || Type == "sp"){
    foundlings1 <- as.vector(unlist(strsplit(as.character(SearchingData1[which(SearchingData1$Species %in% Terms1),3]), split = ",")))
    foundlings2 <- as.vector(unlist(strsplit(as.character(SearchingData2[which(SearchingData2$PhenoWord %in% Terms2),6]), split = ",")))
    foundlings <-  foundlings1[foundlings1 %in% foundlings2]
  }else if(Type == "PS" || Type == "ps" ){
    foundlings1 <- as.vector(unlist(strsplit(as.character(SearchingData1[which(SearchingData1$PhenoWord %in% Terms1),6]), split = ",")))
    foundlings2 <- as.vector(unlist(strsplit(as.character(SearchingData2[which(SearchingData2$Gene %in% Terms2),3]), split = ",")))
    foundlings <-  foundlings1[foundlings1 %in% foundlings2]
  }else{
    stop("ERROR: 'Type' must be either 'SG', 'GS', 'GP', 'PG', 'SP' or 'PS'!")
  }
  if(length(foundlings) == 0)
  {
    stop("No Matches Found")
  }
  return(foundlings)
}

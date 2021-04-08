#' Creates Artificial Gene Broups Based On Non-Numerical Names
#' @param input character vector or data.frame: input
#' @export
# Written by John M. A. Wojahn August 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University

GeneNamesGroupeR <- function(input)
{
  if(class(input) == "character")
  {
    print("Processing vector, GeneGroups will be returned as vector")
    input <- as.vector(input)
    GeneGroups <- c(1:length(input))
    GeneGroups[] <- NA
    for(i in 1:length(input))
    {
      thisgene <- as.character(input[i])
      if(grepl("[[:digit:]]",thisgene))
      {
        counter <- 0
        while(substring(thisgene,nchar(thisgene)-counter,nchar(thisgene)-counter) %in% as.character(0:9))
        {
          counter <- counter + 1
        }
        GeneGroups[i] <- substring(thisgene, 1, nchar(thisgene)-counter)
      }else{
        GeneGroups[i] <- thisgene
      }
    }
    output <- GeneGroups
  }else if(class(input) == "data.frame"){
    print("Processing data.frame, GeneGroups will be appended to right of data.frame")
    input <- as.vector(input$Gene)
    GeneGroups <- c(1:length(input))
    GeneGroups[] <- NA
    for(i in 1:length(input))
    {
      thisgene <- as.character(input[i])
      if(grepl("[[:digit:]]",thisgene))
      {
        counter <- 0
        while(substring(thisgene,nchar(thisgene)-counter,nchar(thisgene)-counter) %in% as.character(0:9))
        {
          counter <- counter + 1
        }
        GeneGroups[i] <- substring(thisgene, 1, nchar(thisgene)-counter)
      }else{
        GeneGroups[i] <- thisgene
      }
    }
    input[,(ncol(input) + 1)] <- GeneGroups
    colnames(input)[ncol(input) + 1] <- "GeneNamesGroup"
    output <- input
  }else{
     stop("ERROR: input must be either a character vector or a data.frame!")
  }
  return(output)
}

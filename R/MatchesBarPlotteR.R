#' Makes Matches Barplot
#' @param strings character vector: abstracts text
#' @param matches character vector: comma separated matched IDs
#' @param n numeric: how many columns, either numeric or character of "all"
#' @export
# Written by John M. A. Wojahn July 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University

MatchesBarPlotteR <- function(strings,matches,n)
{
  strings <- as.vector(strings)
  matches <- as.vector(matches)
  uniquestrings <- unique(strings)
  ForBarplot <- as.data.frame(matrix(nrow = length(uniquestrings), ncol=2))
  ForBarplot[,1] <- uniquestrings
  pb <- txtProgressBar(min = 1, max = nrow(ForBarplot), style = 3)
  for(i in 1:nrow(ForBarplot))
  {
    setTxtProgressBar(pb, i)
    thismatches <- paste(unique(as.character(matches[which(as.character(strings) == as.character(ForBarplot[i,1]))])), collapse = ",")
    ForBarplot[i,2] <- nchar(thismatches) - nchar(gsub(",","",thismatches)) + 1
  }
  ForBarplot <- as.data.frame(ForBarplot[order(ForBarplot[,2], decreasing = T),])
  if(n != "all")
  {
    ForBarplot <- as.data.frame(ForBarplot[1:n,])
  }
  return(ForBarplot)
}

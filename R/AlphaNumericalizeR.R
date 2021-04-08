#' Removes Non-Alphanumerical Characters From Abstracts
#' @param AbstractStrings character vector: of abstracts as strings
#' @export
# Written by John M. A. Wojahn August 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University

AlphaNumericalizeR <- function(AbstractStrings)
{
  #check classes
  if(class(AbstractStrings) != "character")
  {
    stop("ERROR: AbstractStrings must be a character vector!")
  }
  print(sprintf("Removing all non-alphanumerical characters from %s abstracts",length(AbstractStrings)))
  AbstractStrings <- as.vector(AbstractStrings)
  AbstractsStringsFixed <- gsub("[^[:alnum:] ]", "", AbstractStrings)
  return(AbstractsStringsFixed)
}

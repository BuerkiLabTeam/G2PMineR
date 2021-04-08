#' RemovesHTML elements
#' @param AbstractsStrings character vector: character vector: of abstracts as strings
#' @export
# Written by John M. A. Wojahn August 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University

HTMLElementRemoveR <- function(AbstractsStrings)
{
  if(class(AbstractsStrings) != "character")
  {
    stop("ERROR: Input must be a character vector!")
  }
  if(length(AbstractsStrings) == 0)
  {
    stop("ERROR: Input must not be empty!")
  }
  tmp <- gsub("<b>","",as.character(AbstractsStrings))
  tmp <- gsub("</b>","",as.character(tmp))
  tmp <- gsub("<em>","",as.character(tmp))
  tmp <- gsub("</em>","",as.character(tmp))
  tmp <- gsub("<mark>","",as.character(tmp))
  tmp <- gsub("</mark>","",as.character(tmp))
  tmp <- gsub("<ins>","",as.character(tmp))
  tmp <- gsub("</ins>","",as.character(tmp))
  tmp <- gsub("<del>","",as.character(tmp))
  tmp <- gsub("</del>","",as.character(tmp))
  tmp <- gsub("<i>","",as.character(AbstractsStrings))
  tmp <- gsub("</i>","",as.character(tmp))
  tmp <- gsub("<sup>","",as.character(tmp))
  tmp <- gsub("</sup>","",as.character(tmp))
  tmp <- gsub("<sub>","",as.character(tmp))
  tmp <- gsub("</sub>","",as.character(tmp))
  tmp <- gsub("<strong>","",as.character(tmp))
  tmp <- gsub("</strong>","",as.character(tmp))
  tmp <- gsub("<small>","",as.character(tmp))
  tmp <- gsub("</small>","",as.character(tmp))
  return(tmp)
}

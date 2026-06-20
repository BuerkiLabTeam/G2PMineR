#' Gets Abstracts Frm PubMed IDs
#' @param IDs character or numeric vector: unique abstract IDs
#' @export
# Written by John M. A. Wojahn June 2020, Updated June 2026 (Awesome AWD)
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University

AbstractsGetteR <- function(IDs)
{
  if (!is.character(IDs) && !is.numeric(IDs)) {
    stop("ERROR: IDs must be a character or numeric vector!")
  }
  IDs <- as.character(as.vector(IDs))
  outvec <- character(length(IDs))
  
  for (l in 1:length(IDs)) {
    # 1. Take a full 1.5-second breath to fly completely under NCBI's firewall radar
    Sys.sleep(1.5) 
    print(sprintf("Processing PubMedID No. %s: %s", l, IDs[l]))
    
    # 2. Build the exact raw E-utilities URL for fetching an abstract in plain text
    api_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
                      "db=pubmed&id=", IDs[l], "&retmode=text&rettype=abstract")
    
    # 3. Read the web text directly bypassing heavy R package wrappers
    raw_text <- tryCatch({
      suppressWarnings(readLines(api_url, warn = FALSE))
    }, error = function(e) { NULL })
    
    # 4. Process the raw text dump
    if (!is.null(raw_text) && length(raw_text) > 0) {
      # Collapse the lines into a single block of text
      full_abstract <- paste(raw_text, collapse = " ")
      
      # Clean up any trailing/leading white space
      full_abstract <- trimws(full_abstract)
      
      if (nchar(full_abstract) > 50) {
        print("Isolating Abstract... IN PUBMED")
        outvec[l] <- full_abstract
      } else {
        print("NOT IN PUBMED VERIFIED (Empty Response)")
        outvec[l] <- "NOT_IN_PUBMED"
      }
    } else {
      print("NOT IN PUBMED VERIFIED (NCBI Connection Dropped)")
      outvec[l] <- "NOT_IN_PUBMED"
    }
  }
  return(outvec)
}

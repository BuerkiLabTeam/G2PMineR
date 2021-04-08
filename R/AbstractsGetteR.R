#' Gets Abstracts Frm PubMed IDs
#' @param IDs character or numeric vector: unique abstract IDs
#' @export
# Written by John M. A. Wojahn June 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University

AbstractsGetteR <- function(IDs)
{
  if(class(IDs) != "character" && class(IDs) != "numeric")
  {
    stop("ERROR: IDs must be a character or numeric vector!")
  }
  IDs <- as.numeric(as.vector(IDs))
  # Query pubmed and fetch many results
  outvec <- c(1:length(IDs))
  for(l in 1:length(IDs))
  {
    Sys.sleep(1)
    print(sprintf("Processing PubMedID No. %s",l))
    print(IDs[l])
    PubMed_Query <- easyPubMed::get_pubmed_ids(IDs[l])
    if(length(PubMed_Query[["IdList"]][["Id"]]) >= 1)
    {
      # Fetch data
      print("Mining PubMed Data...")
      my_abstracts_xml <- easyPubMed::fetch_pubmed_data(PubMed_Query,
                                                        retstart = 0,
                                                        retmax = 5000)
      # store Pubmed Records as elements of a list
      print("Coercing Articles to list...")
      all_xml <- easyPubMed::articles_to_list(my_abstracts_xml)
      # Perform operation (use lapply here, no further parameters)
      print("Converting...")
      final_df <- do.call(rbind, lapply(all_xml, easyPubMed::article_to_df,
                                  max_chars = -1, getAuthors = FALSE))
      print("Converting...")
      # Perform operation (use lapply here, no further parameters)
      final_authors_df <- do.call(rbind, lapply(all_xml, easyPubMed::article_to_df,
                                  max_chars = -1, getAuthors = TRUE))
      print("Formatting...")
      if(length(final_authors_df) != 0)
      {
        id_abst_journal_year_names_df <- final_authors_df[,c("pmid",
                                                             "abstract",
                                                             "year", "journal",
                                                             "lastname",
                                                             "firstname")]
        print("Isolating Abstract...")
        outvec[l] <- as.character(id_abst_journal_year_names_df[1,2])
        print("IN PUBMED")
      }else{
        print("NOT IN PUBMED VERIFIED")
        outvec[l] <- "NOT_IN_PUBMED"
      }
    }else{
        print("NOT IN PUBMED VERIFIED")
        outvec[l] <- "NOT_IN_PUBMED"
    }
  }
  return(outvec)
}

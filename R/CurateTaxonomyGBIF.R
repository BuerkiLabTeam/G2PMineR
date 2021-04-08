#' Curates Taxonomy By Checking For Synonymy And Publication Validity
#' @param mainDirect Main directory
#' @param spp a data.frame of species "Genus spp"
#' @param nomen chosen file name for output
#' @export
# Written by John M. A. Wojahn August 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University


CurateTaxonomy <- function(list_of_spp)
{
  GBIF_data_set_spp <- as.vector(list_of_spp)
  unique_GBIF_data_set_spp <- as.vector(unique(GBIF_data_set_spp))
  for(i in 1:length(unique_GBIF_data_set_spp))
  {
    print(paste0("~~~Curating GBIF_data_set Taxon No. ",i," of ",length(unique_GBIF_data_set_spp),"~~~"))
    if(unique_GBIF_data_set_spp[i] == "")
    {
      GBIF_data_set_spp[which(GBIF_data_set_spp == unique_GBIF_data_set_spp[i])] <- "BAD"
    }else{
      CorrectedSpp <- Taxonstand::TPL(as.character(unique_GBIF_data_set_spp[i]),
        infrasp = NULL,infra = TRUE, corr = TRUE, diffchar = 2, max.distance = 1,
        version = "1.1", encoding = "UTF-8", author = TRUE,
        drop.lower.level = FALSE, file = "", silent = TRUE, repeats = 6)
      CorrectedSpp <- as.data.frame(CorrectedSpp)
      CorrectedSpp <- as.character(CorrectedSpp$Taxon)
      SplitSpp <- as.data.frame(Biostrings::strsplit(CorrectedSpp," "))
      taxize_results <- taxize::ipni_search(family = NULL, infrafamily = NULL,
        genus = as.character(SplitSpp[1,1]),
        infragenus = NULL, species = as.character(SplitSpp[2,1]), infraspecies = NULL,
        publicationtitle = NULL, authorabbrev = NULL,
        includepublicationauthors = NULL, includebasionymauthors = NULL,
        geounit = NULL, addedsince = NULL, modifiedsince = NULL,
        isapnirecord = NULL, isgcirecord = NULL, isikrecord = NULL,
        ranktoreturn = NULL, output = "minimal")
      taxize_results <- as.data.frame(taxize_results)
      if(nrow(taxize_results) == 0)
      {
        GBIF_data_set_spp[which(GBIF_data_set_spp == unique_GBIF_data_set_spp[i])] <- "BAD"

      }else{
        GBIF_data_set_spp[which(GBIF_data_set_spp == unique_GBIF_data_set_spp[i])] <- as.character(CorrectedSpp)
      }
      gc()
    }
  }
  curated_list_of_spp <- as.vector(GBIF_data_set_spp)
  return(curated_list_of_spp)
}

#2019 John Michael Adrian Wojahn
#All rights reserved/Tous droits réservés
#CeCILL-B Free Software License Agreement/
#Contrat de Licence de Logiciel Libre CeCILL-B

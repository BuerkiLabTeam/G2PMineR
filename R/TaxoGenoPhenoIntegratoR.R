#' Integrate Matches For Taxonomy, Genes, and Phenotypes Analyses
#' @param AbstractsSpp data.frame: output of SpeciesLookeR
#' @param GenesOut data.frame: output of GenesLookeR
#' @param AbsPhen data.frame: output of PhenotypeLookeR
#' @param MakeVenn boolean: make venn diagram of intersection?
#' @export
# Written by John M. A. Wojahn September 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University

TaxoGenoPhenoIntegratoR <- function(AbstractsSpp,GenesOut,AbsPhen,MakeVenn)
{
  #check classes
  if(class(AbstractsSpp) != "data.frame")
  {
    stop("ERROR: AbstractsSpp must be a data.frame!")
  }
  if(class(GenesOut) != "data.frame")
  {
    stop("ERROR: GenesOut must be a data.frame!")
  }
  if(class(AbsPhen) != "data.frame")
  {
    stop("ERROR: ERROR: AbsPhen must be a data.frame!")
  }
  print("Making Typewise Unique IDs Vector")
  AbstractsSppMatches <- unique(as.vector(unlist(strsplit(as.character(AbstractsSpp$Matches),split=","))))
  GenesOutMatches <- unique(as.vector(unlist(strsplit(as.character(GenesOut$Matches),split=","))))
  AbsPhenMatches <- unique(as.vector(unlist(strsplit(as.character(AbsPhen$AbsMatches),split=","))))
  print("Making Typewise Non-Unique IDs Lists")
  ListedAbstractsSppMatches <- strsplit(as.character(AbstractsSpp$Matches),split=",")
  ListedGenesOutMatches <- strsplit(as.character(GenesOut$Matches),split=",")
  ListedAbsPhenMatches <- strsplit(as.character(AbsPhen$AbsMatches),split=",")
  print("Making Unique Combination IDs Vector")
  AllUniqueAbs <- unique(c(AbstractsSppMatches,GenesOutMatches,AbsPhenMatches))
  print("Performing Three-Dimensional Integration")
  IntegratorOut <- as.data.frame(matrix(nrow=length(AllUniqueAbs),ncol=8))
  colnames(IntegratorOut) <- c("ID","TaxonomyCoocs","GenesCoocs","PhenotypesCoocs","TaxonomyCoocsWords","GenesCoocsWords","PhenotypesCoocsWords","Tripartite")
  IntegratorOut[,1] <- AllUniqueAbs
  pb <- txtProgressBar(min = 1, max = nrow(IntegratorOut), style = 3)
  for(i in 1:nrow(IntegratorOut))
  {
    setTxtProgressBar(pb, i)
    #print(paste0("Processing No. ",i," of ",nrow(IntegratorOut)))
    #taxonomy
    if(length(as.character(AbstractsSpp[which(as.character(ListedAbstractsSppMatches) %in% as.character(IntegratorOut[i,1])),3])) != 0)
    {
      IntegratorOut[i,2] <- paste0(as.character(AbstractsSpp[which(as.character(ListedAbstractsSppMatches) %in% as.character(IntegratorOut[i,1])),3]), collapse = "@")
      IntegratorOut[i,5] <- paste0(as.character(AbstractsSpp[which(as.character(ListedAbstractsSppMatches) %in% as.character(IntegratorOut[i,1])),2]), collapse = "@")
    }else{
      IntegratorOut[i,2] <- "NONE"
      IntegratorOut[i,5] <- "NONE"
    }
    if(length(as.character(GenesOut[which(as.character(ListedGenesOutMatches) %in% as.character(IntegratorOut[i,1])),3])) != 0)
    {
      #genes
      IntegratorOut[i,3] <- paste0(as.character(GenesOut[which(as.character(ListedGenesOutMatches) %in% as.character(IntegratorOut[i,1])),3]), collapse = "@")
      IntegratorOut[i,6] <- paste0(as.character(GenesOut[which(as.character(ListedGenesOutMatches) %in% as.character(IntegratorOut[i,1])),1]), collapse = "@")
    }else{
      IntegratorOut[i,3] <- "NONE"
      IntegratorOut[i,6] <- "NONE"
    }
    if(length(as.character(AbsPhen[which(as.character(ListedAbsPhenMatches) %in% as.character(IntegratorOut[i,1])),6])) != 0)
    {
      #phenotypes
      IntegratorOut[i,4] <- paste0(as.character(AbsPhen[which(as.character(ListedAbsPhenMatches) %in% as.character(IntegratorOut[i,1])),6]), collapse = "@")
      IntegratorOut[i,7] <- paste0(as.character(AbsPhen[which(as.character(ListedAbsPhenMatches) %in% as.character(IntegratorOut[i,1])),1]), collapse = "@")
    }else{
      IntegratorOut[i,4] <- "NONE"
      IntegratorOut[i,7] <- "NONE"
    }
    if("NONE" %in% IntegratorOut[i,])
    {
      IntegratorOut[i,8] <- "No"
    }else{
      IntegratorOut[i,8] <- "Yes"
    }
  }
  if(MakeVenn == T)
  {
    vennplot <- plot(eulerr::euler(list(WithTaxonomy = AbstractsSppMatches, WithGenes = GenesOutMatches, WithPhenotypes = AbsPhenMatches)))
    print("Returning Outmatrix as first object and Venn as second in list")
    outlist <- list(IntegratorOut,vennplot)
    return(outlist)
  }else{
    print("Returning Outmatrix")
    return(IntegratorOut)
  }
}

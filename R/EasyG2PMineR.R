#' Easy G2PMineR Analysis Function
#' @param Query character: Search terms
#' @param Kingdom character: which kingdom to use
#' @param Seed numeric: a number
#' @param HowMany numeric: max number of abstracts to sample
#' @export
# Written by John M. A. Wojahn August 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU Affero General Public License v. 3
# Funded by EPSCoR GEM3 at Boise State University

EasyG2PMineR <- function(Query, Kingdom, Seed, HowMany)
{
  {
    message("*********************************************************")
    message("** Welcome to Version 3 (Happy Hyena) of EasyG2PMineR! **")
    message("********* By John. M. A. Wojahn and Sven Buerki *********")
    message("********         ©2023 GNU AGPL v.3            **********")
    message("*********************************************************")
  }
  Sys.sleep(5)
  oldo <- getOption("warn", default = NULL)
  options(warn=1)
  if(length(Query) > 1)
  {
    stop("ERROR: Query must be of length 1!")
  }
  if(class(Query) != "character")
  {
    stop("ERROR: Query must be of class character!")
  }
  if(!(Kingdom %in% c("P","A","F")))
  {
    stop("ERROR: Kingdom must be either 'P','A',or 'F' [i.e. Plantae, Animalia, Fungi]!")
  }
  if(class(Seed) != "numeric")
  {
    stop("ERROR: Seed must be a number [no quotes also]!")
  }
  if(class(HowMany) != "numeric")
  {
    stop("ERROR: HowMany must be a number [no quotes also]!")
  }
  if(HowMany < 10)
  {
    stop("ERROR: HowMany must be 10 or greater!")
  }
  if(Kingdom == "P")
  {
    MyKing <- "P"
  }else if(Kingdom == "A"){
    MyKing <- "A"
  }else if(Kingdom == "F"){
    MyKing <- "F"
  }
  # Search
  message("Searching PubMed using Query")
  res <- suppressWarnings(RISmed::EUtilsSummary(Query,
                             type="esearch", db="pubmed",
                             retmax=10000))
  message("Getting PubMed IDs")
  IDs <- attr(res,"PMID")
  if(length(IDs) > 0)
  {
    #select a different seed number for each group
    set.seed(Seed)
    message("Pseudorandomly Sampling PubMed IDs")
    if(length(IDs) > HowMany)
    {
      IDs <- IDs[sample(1:length(IDs), size = HowMany, replace = F)]
    }
    IDs <- as.data.frame(IDs)
    IDs <- as.numeric(as.character(IDs[,1]))
    message("Extracting abstracts text")
    AbstractsStrings <- G2PMineR::AbstractsGetteR(IDs)
    IDs <- IDs[!is.na(AbstractsStrings)]
    AbstractsStrings <- AbstractsStrings[!is.na(AbstractsStrings)]
    message("Cleaning abstracts")
    # Cleaning Abstracts
    AbstractsStrings <- G2PMineR::HTMLElementRemoveR(AbstractsStrings)
    AbstractsStrings <- G2PMineR::AlphaNumericalizeR(AbstractsStrings)
    # Mining Abstracts for Taxonomy
    message("Searching abstracts for species")
    AbstractsSpp <- G2PMineR::SpeciesLookeR(AbstractsStrings,IDs, Kingdom = MyKing, Add = NULL)

    if(nrow(AbstractsSpp) > 0)
    {
      Sys.sleep(2)
      # Mining Abstracts for Genes
      SpeciesAbbrvs <- G2PMineR::SpeciesAbbreviatoR(AbstractsSpp)
      Sys.sleep(2)
      message("Searching abstracts for genes")
      GenesOut <- G2PMineR::GenesLookeR(AbstractsStrings, IDs, Kingdom = MyKing, Add = NULL, SppAbbr = SpeciesAbbrvs)
      GenesOut <- as.data.frame(GenesOut[!GenesOut$InOrNot == "No",])
      if(nrow(GenesOut) == 1 && is.na(GenesOut[1,1]))
      {
        stop("NO GENES FOUND, BROADEN SEARCH TERMS")
      }else{
        Sys.sleep(2)
        GenesOut <- G2PMineR::SynonymReplaceR(GenesOut, Kingdom = MyKing)
        #to create artificial gene groups
        Sys.sleep(2)
        GeneGroups <- G2PMineR::GeneNamesGroupeR(as.vector(GenesOut[,1]))
        #to grade the usefulness of matches
        Sys.sleep(2)
        GeneGrades <- G2PMineR::UtilityGradeR(GenesOut, Kingdom = MyKing, Add = NULL, Groups=as.data.frame(GeneGroups))
        Sys.sleep(2)
        # Mining Abstracts for Phenotypes
        message("Searching abstracts for phenotypes")
        AbsPhen <- G2PMineR::PhenotypeLookeR(AbstractsStrings, IDs, Kingdom = MyKing, Add = NULL)
        if(nrow(AbsPhen) == 0)
        {
          stop("NO PHENOTYPES FOUND, BROADEN SEARCH TERMS")
        }else{
          # Analyzing genes, taxonomy, and phenotypes data
          #for species matches
          Sys.sleep(2)
          #message("Producing bar plots")
          #SppBarPlotDF <- G2PMineR::MatchesBarPlotteR(AbstractsSpp$Species, AbstractsSpp$Matches,n = 25)
          #for genes matches
          #Sys.sleep(2)
          Genez <- as.data.frame(GenesOut[which(GenesOut$InOrNot == "Yes"),])
          Genez <- data.frame(Genez$Gene,Genez$Matches)
          Genez <- unique(Genez)
          #Sys.sleep(2)
          #GenesBarPlotDF <- G2PMineR::MatchesBarPlotteR(Genez[,1], Genez[,2], n = 25)
          #for phenology matches
          #Sys.sleep(2)
          #PhenoBarPlotDF <- G2PMineR::MatchesBarPlotteR(AbsPhen$PhenoWord,
          #                                              AbsPhen$AbsMatches,n = 25)
          #
          # Internal interrelatins network analysis
          #for species matches
          #message("Conducting internal network analyses")
          #Sys.sleep(2)
          #SppInt <- G2PMineR::InternalPairwiseDistanceInferreR(AbstractsSpp$Species,
          #                                                     AbstractsSpp$Matches,
          #                                                     allabsnum = length(AbstractsStrings))
          #Sys.sleep(2)
          #if(nrow(SppInt) > 1)
          #{
          #  SppIntSmall <- G2PMineR::TopN_PickeR_Internal(SppInt, n = 50, decreasing = T)
          #}else{
          #  SppIntSmall <- SppInt
          #}
          #Sys.sleep(2)
          #rwmsS <- suppressWarnings(as.dist(SppIntSmall, diag = F, upper = FALSE))
          #for genes matches
          #Sys.sleep(2)
          #GenInt <- G2PMineR::InternalPairwiseDistanceInferreR(Genez[,1], Genez[,2],
          #                                                     allabsnum = length(AbstractsStrings))
          #Sys.sleep(2)
          #if(nrow(GenInt) > 1)
          #{
          #  GenIntSmall <- G2PMineR::TopN_PickeR_Internal(GenInt, n = 50, decreasing = T)
          #}else{
          #  GenIntSmall <- GenInt
          #}
          #Sys.sleep(2)
          #rwmsG <- suppressWarnings(as.dist(GenIntSmall, diag = F, upper = FALSE))
          #Sys.sleep(2)
          #for phenology matches
          #PhenInt <- G2PMineR::InternalPairwiseDistanceInferreR(AbsPhen$PhenoWord,
          #                                                      AbsPhen$AbsMatches,
          #                                                      allabsnum = length(AbstractsStrings))
          #Sys.sleep(2)
          #if(nrow(PhenInt) > 1)
          #{
          #  PhenIntSmall <- G2PMineR::TopN_PickeR_Internal(PhenInt, n = 100, decreasing = T)
          #}else{
          #  PhenIntSmall <- PhenInt
          #}
          #Sys.sleep(2)
          #rwmsP <- suppressWarnings(as.dist(PhenIntSmall, diag = FALSE, upper = FALSE))
          #Sys.sleep(2)

          # Linking genome to phenome (within a taxonomic framework)
          #Phenotypes vs Species
          message("Performing bipartite analyses")
          PhenoSpecies <- G2PMineR::PairwiseDistanceInferreR(AbstractsSpp$Species,
                                                             AbstractsSpp$Matches,
                                                             AbsPhen$PhenoWord,
                                                             AbsPhen$AbsMatches,
                                                             allabsnum = length(AbstractsStrings))
          Sys.sleep(2)
          if(nrow(PhenoSpecies) > 1)
          {
            PhenoSpeciesSmall <- G2PMineR::TopN_PickeR(PhenoSpecies, n = 50, decreasing = T)
          }else if(nrow(PhenoSpecies) == 1){
            PhenoSpeciesSmall <- PhenoSpecies
          }else if(nrow(PhenoSpecies) == 0){
            stop("ERROR: No P2S Associations Found, broaden search terms, increase HowMany, or change Seed")
          }
          Sys.sleep(2)
          #Genes vs Species
          GeneSpecies <- G2PMineR::PairwiseDistanceInferreR(AbstractsSpp$Species,
                                                            AbstractsSpp$Matches,
                                                            Genez[,1], Genez[,2],
                                                            allabsnum = length(AbstractsStrings))
          Sys.sleep(2)
          if(nrow(GeneSpecies) > 1)
          {
            GeneSpeciesSmall <- G2PMineR::TopN_PickeR(GeneSpecies, n = 50, decreasing = T)
          }else if(nrow(GeneSpecies) == 1){
            GeneSpeciesSmall <- GeneSpecies
          }else if(nrow(GeneSpecies) == 0){
            stop("ERROR: No G2S Associations Found, broaden search terms, increase HowMany, or change Seed")
          }
          Sys.sleep(2)
          #Phenotypes vs Genes
          PhenoGenes <- G2PMineR::PairwiseDistanceInferreR(AbsPhen$PhenoWord,
                                                           AbsPhen$AbsMatches,
                                                           Genez[,1], Genez[,2],
                                                           allabsnum = length(AbstractsStrings))

          Sys.sleep(2)
          if(nrow(PhenoGenes) > 1)
          {
            PhenoGenesSmall <- G2PMineR::TopN_PickeR(PhenoGenes, n = 50, decreasing = T)
          }else if(nrow(PhenoGenes) == 1){
            PhenoGenesSmall <- PhenoGenes
          }else if(nrow(PhenoGenes) == 0){
            stop("ERROR: No G2P Associations Found, broaden search terms, increase HowMany, or change Seed")
          }
          Sys.sleep(2)
          spmats <- ceiling(G2PMineR::AbstractsProportionCalculator(AbstractsSpp, IDs)*100)
          gnmats <- ceiling(G2PMineR::AbstractsProportionCalculator(GenesOut, IDs)*100)
          phmats <- ceiling(G2PMineR::AbstractsProportionCalculator(AbsPhen, IDs)*100)

          message("Producing output")
          pdf("InterpretableResults.pdf")
          plot.new()
          text(0.3,1,"EasyG2PMineR Output")
          text(0.5,0.75,Query)
          text(0.2,0.9,Sys.Date())
          text(0.5,0.6, sprintf("%s percent abstracts had species matches",spmats))
          text(0.5,0.5, sprintf("%s percent abstracts had gene matches",gnmats))
          text(0.5,0.4, sprintf("%s percent abstracts had phenotype matches",phmats))
          plot.new()
          text(0.3,1,"Bipartite Outputs")
          library(bipartite)
          bipartite::plotweb(PhenoGenesSmall, text.rot=90, col.interaction = "gray",
                             labsize = 0.75, method='normal')
          text(0.05,1.78,"G2P")

          bipartite::plotweb(PhenoSpeciesSmall, text.rot=90, col.interaction = "gray",
                             labsize = 0.75, method='normal')
          text(0.05,1.78,"P2S")

          bipartite::plotweb(GeneSpeciesSmall, text.rot=90, col.interaction = "gray",
                             labsize = 0.75, method='normal')
          text(0.05,1.78,"G2S")
          dev.off()
          write.csv(PhenoGenesSmall, "G2P_data.csv",row.names=F)
        }
      }
    }else{
      message("NO SPECIES FOUND, BROADEN SEARCH TERMS")
    }
  }else{
    message("NO ABSTRACTS FOUND, BROADEN SEARCH TERMS")
}


  options(warn=oldo)
}

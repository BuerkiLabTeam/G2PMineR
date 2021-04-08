#' Explore Inter-Dataset Intersectionality
#' @param Ta data.frame: longform taxonomical results OR NULL
#' @param G data.frame: longform genes results OR NULL
#' @param P data.frame: longform phenotypes results OR NULL
#' @param AbstractsSpp data.frame: shortform taxonomical results OR NULL
#' @param GenesOut data.frame: shortform genes results OR NULL
#' @param AbsPhen data.frame: shortform phenotypes results OR NULL
#' @export
# Written by John M. A. Wojahn January 2021
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University

ConsensusInferreR <- function(Ta, G, P, AbstractsSpp, GenesOut, AbsPhen)
{
  if(is.null(Ta) == F && is.null(G) && is.null(P))
  {
    stop("You must give me more than one dataset to compare!")
  }else if(is.null(Ta) && is.null(G) == F && is.null(P)){
    stop("You must give me more than one dataset to compare!")
  }else if(is.null(Ta) && is.null(G) && is.null(P) == F){
    stop("You must give me more than one dataset to compare!")
  }else if(is.null(Ta) == F && is.null(G) == F && is.null(P)){
    if(is.null(AbstractsSpp) == F && is.null(GenesOut) == F && is.null(AbsPhen))
    {
      TaxoMatchesUnique <- unique(Ta$Matches)
      GenoMatchesUnique <- unique(G$Matches)
      SharedMatches <- TaxoMatchesUnique[TaxoMatchesUnique %in% GenoMatchesUnique]
      OUT <- as.data.frame(matrix(nrow=length(SharedMatches),ncol=3))
      colnames(OUT) <- c("Matches","Genes","Species")
      OUT[,1] <- SharedMatches
      for(i in 1:length(SharedMatches))
      {
        Spp <- unique(Ta[Ta$Matches == SharedMatches[i],2])
        Genes <- unique(G[G$Matches == SharedMatches[i],1])
        OUT[i,2] <- paste(Genes, collapse = ",")
        OUT[i,3] <- paste(Spp, collapse = ",")
      }

      print("Making consensus-only files")
      ConsensusIDs <- as.vector(OUT[,1])
      Ta_Consensus <- as.data.frame(matrix(nrow=0, ncol= ncol(AbstractsSpp)))
      colnames(Ta_Consensus) <- colnames(AbstractsSpp)
      pb <- txtProgressBar(min = 1, max = nrow(AbstractsSpp), style = 3)
      for(i in 1:nrow(AbstractsSpp))
      {
        setTxtProgressBar(pb, i)
        thisrow <- as.vector(unlist(strsplit(AbstractsSpp[i,3], split = ",")))
        goodthisrow <- thisrow[thisrow %in% ConsensusIDs]
        if(length(goodthisrow) != 0)
        {
          fullrow <- as.data.frame(AbstractsSpp[i,])
          fullrow[1,3] <- paste(goodthisrow, collapse = ",")
          Ta_Consensus <- as.data.frame(rbind(Ta_Consensus,fullrow))
        }
      }

      G_Consensus <- as.data.frame(matrix(nrow=0, ncol= ncol(GenesOut)))
      colnames(G_Consensus) <- colnames(GenesOut)
      pb <- txtProgressBar(min = 1, max = nrow(GenesOut), style = 3)
      for(i in 1:nrow(GenesOut))
      {
        setTxtProgressBar(pb, i)
        thisrow <- as.vector(unlist(strsplit(GenesOut[i,3], split = ",")))
        goodthisrow <- thisrow[thisrow %in% ConsensusIDs]
        if(length(goodthisrow) != 0)
        {
          fullrow <- as.data.frame(GenesOut[i,])
          fullrow[1,3] <- paste(goodthisrow, collapse = ",")
          G_Consensus <- as.data.frame(rbind(G_Consensus,fullrow))
        }
      }

      TaxoLength <- length(TaxoMatchesUnique)
      GenoLength <- length(GenoMatchesUnique)
      both <- length(SharedMatches)
      OUTvenn <- VennDiagram::draw.pairwise.venn(TaxoLength, GenoLength, both, c("Abstracts with Species Matches","Abstracts with Genes Matches"))
      OUTlist <- list(OUT,OUTvenn,Ta_Consensus,G_Consensus,ConsensusIDs)
      names(OUTlist) <- c("IntersectionMatrix", "Venn", "Taxonomy Consensus-Only","Genes Consensus-Only","ConsensusIDs")
      print("Output is a list with matrix as first element, venn as second element, and consensus-only files as third (Taxonomy) and fourth (Genes) elements, and ConsensusIDs as the last element")
      print("To draw venn use library(VennDiagram) grid.draw()")
      print("To make bipartite graphs feed the consensus-only data.frames into the pairwise distance functions and bipartite::plotweb function. See vignette for instructions")
      return(OUTlist)
    }else{
      stop("You must give me the shortform datasets too!")
    }
  }else if(is.null(Ta) && is.null(G) == F && is.null(P) == F){
    if(is.null(AbstractsSpp) && is.null(GenesOut) == F && is.null(AbsPhen) == F)
    {
      PhenoMatchesUnique <- unique(P$AbsMatches)
      GenoMatchesUnique <- unique(G$Matches)
      SharedMatches <- PhenoMatchesUnique[PhenoMatchesUnique %in% GenoMatchesUnique]
      OUT <- as.data.frame(matrix(nrow=length(SharedMatches),ncol=3))
      colnames(OUT) <- c("Matches","Genes","Species")
      OUT[,1] <- SharedMatches
      for(i in 1:length(SharedMatches))
      {
        Phens <- unique(P[P$AbsMatches == SharedMatches[i],1])
        Genes <- unique(G[G$Matches == SharedMatches[i],1])
        OUT[i,2] <- paste(Genes, collapse = ",")
        OUT[i,3] <- paste(Phens, collapse = ",")
      }

      print("Making consensus-only files")
      ConsensusIDs <- as.vector(OUT[,1])
      G_Consensus <- as.data.frame(matrix(nrow=0, ncol= ncol(GenesOut)))
      colnames(G_Consensus) <- colnames(GenesOut)
      pb <- txtProgressBar(min = 1, max = nrow(GenesOut), style = 3)
      for(i in 1:nrow(GenesOut))
      {
        setTxtProgressBar(pb, i)
        thisrow <- as.vector(unlist(strsplit(GenesOut[i,3], split = ",")))
        goodthisrow <- thisrow[thisrow %in% ConsensusIDs]
        if(length(goodthisrow) != 0)
        {
          fullrow <- as.data.frame(GenesOut[i,])
          fullrow[1,3] <- paste(goodthisrow, collapse = ",")
          G_Consensus <- as.data.frame(rbind(G_Consensus,fullrow))
        }
      }

      P_Consensus <- as.data.frame(matrix(nrow=0, ncol= ncol(AbsPhen)))
      colnames(P_Consensus) <- colnames(AbsPhen)
      pb <- txtProgressBar(min = 1, max = nrow(AbsPhen), style = 3)
      for(i in 1:nrow(AbsPhen))
      {
        setTxtProgressBar(pb, i)
        thisrow <- as.vector(unlist(strsplit(AbsPhen[i,6], split = ",")))
        goodthisrow <- thisrow[thisrow %in% ConsensusIDs]
        if(length(goodthisrow) != 0)
        {
          fullrow <- as.data.frame(AbsPhen[i,])
          fullrow[1,6] <- paste(goodthisrow, collapse = ",")
          fullrow[1,2] <- length(goodthisrow)
          P_Consensus <- as.data.frame(rbind(P_Consensus,fullrow))
        }
      }

      PhensLength <- length(PhenoMatchesUnique)
      GenoLength <- length(GenoMatchesUnique)
      both <- length(SharedMatches)
      OUTvenn <- VennDiagram::draw.pairwise.venn(PhensLength, GenoLength, both, c("Abstracts with Phenotypes Matches","Abstracts with Genes Matches"))
      OUTlist <- list(OUT,OUTvenn,G_Consensus,P_Consensus,ConsensusIDs)
      names(OUTlist) <- c("IntersectionMatrix", "Venn", "Genes Consensus-Only","Phenotypes Consensus-Only","ConsensusIDs")
      print("Output is a list with matrix as first element, venn as second element, and consensus-only files as third (Genes) and fourth (Phenotypes) elements, and ConsensusIDs as the last element")
      print("To draw venn use library(VennDiagram) grid.draw()")
      print("To make bipartite graphs feed the consensus-only data.frames into the pairwise distance functions and bipartite::plotweb function. See vignette for instructions")
      return(OUTlist)
    }else{
      stop("You must give me the shortform datasets too!")
    }
  }else if(is.null(Ta) == F && is.null(G) && is.null(P) == F){
    if(is.null(AbstractsSpp) == F && is.null(GenesOut) && is.null(AbsPhen) == F)
    {
      PhenoMatchesUnique <- unique(P$AbsMatches)
      TaxoMatchesUnique <- unique(Ta$Matches)
      SharedMatches <- PhenoMatchesUnique[PhenoMatchesUnique %in% TaxoMatchesUnique]
      OUT <- as.data.frame(matrix(nrow=length(SharedMatches),ncol=3))
      colnames(OUT) <- c("Matches","Genes","Species")
      OUT[,1] <- SharedMatches
      for(i in 1:length(SharedMatches))
      {
        Phens <- unique(P[P$AbsMatches == SharedMatches[i],1])
        Spp <- unique(Ta[Ta$Matches == SharedMatches[i],2])
        OUT[i,2] <- paste(Spp, collapse = ",")
        OUT[i,3] <- paste(Phens, collapse = ",")
      }
          print("Making consensus-only files")
      ConsensusIDs <- as.vector(OUT[,1])
      Ta_Consensus <- as.data.frame(matrix(nrow=0, ncol= ncol(AbstractsSpp)))
      colnames(Ta_Consensus) <- colnames(AbstractsSpp)
      pb <- txtProgressBar(min = 1, max = nrow(AbstractsSpp), style = 3)
      for(i in 1:nrow(AbstractsSpp))
      {
        setTxtProgressBar(pb, i)
        thisrow <- as.vector(unlist(strsplit(AbstractsSpp[i,3], split = ",")))
        goodthisrow <- thisrow[thisrow %in% ConsensusIDs]
        if(length(goodthisrow) != 0)
        {
          fullrow <- as.data.frame(AbstractsSpp[i,])
          fullrow[1,3] <- paste(goodthisrow, collapse = ",")
          Ta_Consensus <- as.data.frame(rbind(Ta_Consensus,fullrow))
        }
      }

      P_Consensus <- as.data.frame(matrix(nrow=0, ncol= ncol(AbsPhen)))
      colnames(P_Consensus) <- colnames(AbsPhen)
      pb <- txtProgressBar(min = 1, max = nrow(AbsPhen), style = 3)
      for(i in 1:nrow(AbsPhen))
      {
        setTxtProgressBar(pb, i)
        thisrow <- as.vector(unlist(strsplit(AbsPhen[i,6], split = ",")))
        goodthisrow <- thisrow[thisrow %in% ConsensusIDs]
        if(length(goodthisrow) != 0)
        {
          fullrow <- as.data.frame(AbsPhen[i,])
          fullrow[1,6] <- paste(goodthisrow, collapse = ",")
          fullrow[1,2] <- length(goodthisrow)
          P_Consensus <- as.data.frame(rbind(P_Consensus,fullrow))
        }
      }
      PhensLength <- length(PhenoMatchesUnique)
      TaxoLength <- length(TaxoMatchesUnique)
      both <- length(SharedMatches)
      OUTvenn <- VennDiagram::draw.pairwise.venn(PhensLength, TaxoLength, both, c("Abstracts with Phenotypes Matches","Abstracts with Species Matches"))
      OUTlist <- list(OUT,OUTvenn,Ta_Consensus,P_Consensus,ConsensusIDs)
      names(OUTlist) <- c("IntersectionMatrix", "Venn", "Taxonomy Consensus-Only","Phenotypes Consensus-Only","ConsensusIDs")
      print("Output is a list with matrix as first element, venn as second element, and consensus-only files as third (Taxonomy) amd fourth (Phenotypes) elements, and ConsensusIDs as the last element")
      print("To draw venn use library(VennDiagram) grid.draw()")
      print("To make bipartite graphs feed the consensus-only data.frames into the pairwise distance functions and bipartite::plotweb function. See vignette for instructions")
      return(OUTlist)
    }else{
      stop("You must give me the shortform datasets too!")
    }
  }else if(is.null(Ta) == F && is.null(G) == F && is.null(P) == F){
    if(is.null(AbstractsSpp) == F && is.null(GenesOut) ==F && is.null(AbsPhen) == F)
    {
      PhenoMatchesUnique <- unique(P$AbsMatches)
      TaxoMatchesUnique <- unique(Ta$Matches)
      GenoMatchesUnique <- unique(G$Matches)
      SharedMatches <- PhenoMatchesUnique[PhenoMatchesUnique %in% TaxoMatchesUnique]
      SharedMatches <- SharedMatches[SharedMatches %in% GenoMatchesUnique]
      OUT <- as.data.frame(matrix(nrow=length(SharedMatches),ncol=4))
      colnames(OUT) <- c("Matches","Genes","Species","Phenotypes")
      OUT[,1] <- SharedMatches
      for(i in 1:length(SharedMatches))
      {
        Genes <- unique(G[G$Matches == SharedMatches[i],1])
        Phens <- unique(P[P$AbsMatches == SharedMatches[i],1])
        Spp <- unique(Ta[Ta$Matches == SharedMatches[i],2])
        OUT[i,2] <- paste(Genes, collapse = ",")
        OUT[i,3] <- paste(Spp, collapse = ",")
        OUT[i,4] <- paste(Phens, collapse = ",")
      }

      print("Making consensus-only files")
      ConsensusIDs <- as.vector(OUT[,1])
      Ta_Consensus <- as.data.frame(matrix(nrow=0, ncol= ncol(AbstractsSpp)))
      colnames(Ta_Consensus) <- colnames(AbstractsSpp)
      pb <- txtProgressBar(min = 1, max = nrow(AbstractsSpp), style = 3)
      for(i in 1:nrow(AbstractsSpp))
      {
        setTxtProgressBar(pb, i)
        thisrow <- as.vector(unlist(strsplit(AbstractsSpp[i,3], split = ",")))
        goodthisrow <- thisrow[thisrow %in% ConsensusIDs]
        if(length(goodthisrow) != 0)
        {
          fullrow <- as.data.frame(AbstractsSpp[i,])
          fullrow[1,3] <- paste(goodthisrow, collapse = ",")
          Ta_Consensus <- as.data.frame(rbind(Ta_Consensus,fullrow))
        }
      }

      G_Consensus <- as.data.frame(matrix(nrow=0, ncol= ncol(GenesOut)))
      colnames(G_Consensus) <- colnames(GenesOut)
      pb <- txtProgressBar(min = 1, max = nrow(GenesOut), style = 3)
      for(i in 1:nrow(GenesOut))
      {
        setTxtProgressBar(pb, i)
        thisrow <- as.vector(unlist(strsplit(GenesOut[i,3], split = ",")))
        goodthisrow <- thisrow[thisrow %in% ConsensusIDs]
        if(length(goodthisrow) != 0)
        {
          fullrow <- as.data.frame(GenesOut[i,])
          fullrow[1,3] <- paste(goodthisrow, collapse = ",")
          G_Consensus <- as.data.frame(rbind(G_Consensus,fullrow))
        }
      }

      P_Consensus <- as.data.frame(matrix(nrow=0, ncol= ncol(AbsPhen)))
      colnames(P_Consensus) <- colnames(AbsPhen)
      pb <- txtProgressBar(min = 1, max = nrow(AbsPhen), style = 3)
      for(i in 1:nrow(AbsPhen))
      {
        setTxtProgressBar(pb, i)
        thisrow <- as.vector(unlist(strsplit(AbsPhen[i,6], split = ",")))
        goodthisrow <- thisrow[thisrow %in% ConsensusIDs]
        if(length(goodthisrow) != 0)
        {
          fullrow <- as.data.frame(AbsPhen[i,])
          fullrow[1,6] <- paste(goodthisrow, collapse = ",")
          fullrow[1,2] <- length(goodthisrow)
          P_Consensus <- as.data.frame(rbind(P_Consensus,fullrow))
        }
      }

      vennlist <- list(Genes = as.vector(GenoMatchesUnique),Species = as.vector(TaxoMatchesUnique),Phenotypes = as.vector(PhenoMatchesUnique))
      OUTvenn <- VennDiagram::venn.diagram(vennlist, filename = NULL)
      OUTlist <- list(OUT,OUTvenn,Ta_Consensus,G_Consensus,P_Consensus,ConsensusIDs)
      names(OUTlist) <- c("IntersectionMatrix", "Venn", "Taxonomy Consensus-Only","Genes Consensus-Only", "Phenotypes Consensus-Only","ConsensusIDs")
      print("Output is a list with matrix as first element, venn as second element, and consensus-only files as third (Taxonomy), fourth (Genes), and fifth (Phenotypes) elements, and ConsensusIDs as the last element")
      print("To draw venn use library(VennDiagram) grid.draw()")
      print("To make bipartite graphs feed the consensus-only data.frames into the pairwise distance functions and bipartite::plotweb function. See vignette for instructions")
      return(OUTlist)
    }else{
      stop("You must give me the shortform datasets too!")
    }
  }else if(is.null(Ta) && is.null(G) && is.null(P)){
    stop("You gave me only nulls.  Please feed me something different!")
  }
}

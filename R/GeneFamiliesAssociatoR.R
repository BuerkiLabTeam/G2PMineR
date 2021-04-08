#' Preprocesses genes (OPTIONAL)
#' @param SwissGenesRaw matrix: swissprot genes
#' @export
# Written by John M. A. Wojahn June 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University

GeneFamiliesAssociatoR <- function(SwissGenesRaw)
{
  print("...cleaning family nomenclature")
  SwissGenesRaw <- as.data.frame(SwissGenesRaw)
  fams <- c(1:length(SwissGenesRaw[,6]))
  fams[] <- NA
  pb <- txtProgressBar(min = 1, max = nrow(SwissGenesRaw), style = 3)
  for(i in 1:nrow(SwissGenesRaw))
  {
    setTxtProgressBar(pb, i)
    if(as.character(SwissGenesRaw[i,6]) == "")
    {
      fams[i] <- paste0(SwissGenesRaw[i,5]," family (one gene)")
    }else{
      fams[i] <-as.character(SwissGenesRaw[i,6])
    }
  }
  SwissGenesRaw[,6] <- as.character(fams)
  SwissGenesRaw <- SwissGenesRaw[!(SwissGenesRaw$Gene.names...primary.. == ""),]
  SwissGenesRaw[SwissGenesRaw$Gene.names...synonym.. == ";",4] <- ""
  Synonyms <- c()
  Accepteds <- c()
  Fams <- c()
  Ontz <- c()
  SoleAccepteds <- c()
  SoleFams <- c()
  SoleOntz <- c()
  print("...inferring family names")
  pb <- txtProgressBar(min = 1, max = nrow(SwissGenesRaw), style = 3)
  for(i in 1:nrow(SwissGenesRaw))
  {
    setTxtProgressBar(pb, i)
    #print(paste0(i," of ",nrow(SwissGenesRaw)))
    if(as.character(SwissGenesRaw[i,4]) != "")
    {
      splitz <- as.vector(unlist(strsplit(as.character(SwissGenesRaw[i,4]), split = " ")))
      acceptedrepped <- rep_len(as.character(SwissGenesRaw[i,5]), length(splitz))
      famsrepped <- rep_len(as.character(SwissGenesRaw[i,6]), length(splitz))
      ontzrepped <- rep_len(as.character(SwissGenesRaw[i,3]), length(splitz))
      Synonyms <- c(Synonyms,splitz)
      Accepteds <- c(Accepteds,acceptedrepped)
      Fams <- c(Fams,famsrepped)
      Ontz <- c(Ontz,ontzrepped)
    }else{
      SoleAccepteds <- c(SoleAccepteds,as.character(SwissGenesRaw[i,5]))
      SoleFams <- c(SoleFams,as.character(SwissGenesRaw[i,6]))
      SoleOntz <- c(SoleOntz,as.character(SwissGenesRaw[i,3]))
    }
  }
  AccSyns <- data.frame(Synonyms,Accepteds,Fams,Ontz)
  AccSyns <- AccSyns[!(AccSyns$Synonyms == ";"),]
  AccSyns <- unique(AccSyns)
  syntab <- as.data.frame(table(AccSyns$Synonyms))
  bads <- unique(syntab[which(syntab$Freq > 1),1])
  AccSynsGoods <- as.data.frame(AccSyns[!(AccSyns$Synonyms %in% bads),])
  SoleAccFams <- data.frame(SoleAccepteds,SoleFams,SoleOntz)
  AccPartFams <- as.data.frame(AccSynsGoods[,-1])
  SynPartFams <- as.data.frame(AccSynsGoods[,-2])
  colnames(SoleAccFams) <- c("SwissGenes", "SwissGenesFamilies", "SwissOnts")
  colnames(AccPartFams) <- c("SwissGenes", "SwissGenesFamilies", "SwissOnts")
  colnames(SynPartFams) <- c("SwissGenes", "SwissGenesFamilies", "SwissOnts")
  SwissGenesCombo <- rbind(SoleAccFams,AccPartFams,SynPartFams)
  SwissGenesCombo <- unique(SwissGenesCombo)
  print("...removing ambiguous gene names")
  gt <- as.data.frame(table(SwissGenesCombo$SwissGenes))
  ambiguousgenes <- as.character(gt[which(gt$Freq > 1),1])
  SwissGenesCombo <- as.data.frame(SwissGenesCombo[-which(SwissGenesCombo$SwissGenes %in% ambiguousgenes),])
  GeneSynonymyKey <- as.data.frame(AccSynsGoods)
  print("...making output list")
  OutList <- list(SwissGenesCombo,GeneSynonymyKey)
  return(OutList)
}

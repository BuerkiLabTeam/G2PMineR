#' Replaces Gene Synonyms With Accepted Nomenclature
#' @param GenesOut matrix: swissprot genes abstracts search results
#' @param Kingdom boolean: which kingdom to use
#' @export
# Written by John M. A. Wojahn July-August 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University


SynonymReplaceR <- function(GenesOut, Kingdom)
{

  require(G2PMineR)

  if(Kingdom == "P")
  {
    print("Using Plantae Internal Data")
  }else if(Kingdom == "A"){
    print("Using Animalia Internal Data")
    GeneSynonymyKey <- AnimalGeneSynonymyKey
  }else if(Kingdom == "F"){
    print("Using Fungi Internal Data")
    GeneSynonymyKey <- FungiGeneSynonymyKey
  }else if(Kingdom == "H"){
    message("Using Human Internal Data")
    GeneSynonymyKey <- HumanGeneSynonymyKey
  }else{
    stop("ERROR: You need to choose a kingdom!")
  }
  GenesOut <- as.data.frame(GenesOut)
  GenesOut <- as.data.frame(GenesOut[!(GenesOut$InOrNot == "No"),])
  GeneSynonymyKey <- as.data.frame(GeneSynonymyKey)
  NewGenesOut <- as.data.frame(matrix(nrow=0, ncol=ncol(GenesOut)))
  colnames(NewGenesOut) <- colnames(GenesOut)
  print("...replacing synonyms")
  pb <- txtProgressBar(min = 1, max = nrow(GenesOut), style = 3)
  for(i in 1:nrow(GenesOut))
  {
    setTxtProgressBar(pb, i)
    NewGenesOutSmall <- as.data.frame(matrix(nrow=0, ncol=ncol(GenesOut)))
    colnames(NewGenesOutSmall) <- colnames(GenesOut)
    if(as.character(GenesOut[i,1]) %in% GeneSynonymyKey$Synonyms)
    {
      AcceptedNomen <- as.character(GeneSynonymyKey[(GeneSynonymyKey$Synonyms == as.character(GenesOut[i,1])),2])
      if(length(AcceptedNomen) > 1)
      {
        AcceptedNomen <- sprintf("AMBIGUOUS SYNONYM: %s",paste(AcceptedNomen, collapse = " or "))
      }
      NewGenesOutSmall[1,1] <- AcceptedNomen
      NewGenesOutSmall[1,2] <- as.character(GenesOut[i,2])
      NewGenesOutSmall[1,3] <- as.character(GenesOut[i,3])
      NewGenesOutSmall[1,4] <- as.character(GenesOut[i,4])
      NewGenesOutSmall[1,5] <- as.character(GenesOut[i,5])
      NewGenesOutSmall[1,6] <- as.character(GenesOut[i,6])
    }else{
      NewGenesOutSmall[1,1] <- as.character(GenesOut[i,1])
      NewGenesOutSmall[1,2] <- as.character(GenesOut[i,2])
      NewGenesOutSmall[1,3] <- as.character(GenesOut[i,3])
      NewGenesOutSmall[1,4] <- as.character(GenesOut[i,4])
      NewGenesOutSmall[1,5] <- as.character(GenesOut[i,5])
      NewGenesOutSmall[1,6] <- as.character(GenesOut[i,6])
    }
    NewGenesOut <- as.data.frame(rbind(NewGenesOut,NewGenesOutSmall))
  }
  print("...Combining synames")
  genenums <- as.data.frame(table(NewGenesOut[,1]))
  if(length(which(genenums$Freq >1)) != 0)
  {
    tocombines <- genenums[which(genenums$Freq >1),1]
    NewOutComboSmall <- as.data.frame(NewGenesOut[-which(NewGenesOut$Gene %in% tocombines),])
    if(length(tocombines) > 0)
    {
      #pb <- txtProgressBar(min = 1, max = length(tocombines), style = 3)
      for(i in 1:length(tocombines))
      {
      #setTxtProgressBar(pb, i)
      SmallOut <- as.data.frame(matrix(nrow=1,ncol=ncol(NewOutComboSmall)))
      colnames(SmallOut) <- colnames(NewOutComboSmall)
      tmp <- NewGenesOut[which(NewGenesOut$Gene == as.character(tocombines[i])),]
      SmallOut [1,1] <- paste(unique(as.vector(unlist(strsplit(as.character(tmp$Gene),split=",")))), collapse = ",")
      SmallOut [1,2] <- paste(unique(as.vector(unlist(strsplit(as.character(tmp$InOrNot),split=",")))), collapse = ",")
      SmallOut [1,3] <- paste(unique(as.vector(unlist(strsplit(as.character(tmp$Matches),split=",")))), collapse = ",")
      SmallOut [1,4] <- paste(unique(as.vector(unlist(strsplit(as.character(tmp$InSitus),split=",")))), collapse = ",")
      SmallOut [1,5] <- paste(unique(as.vector(unlist(strsplit(as.character(tmp$Family),split=",")))), collapse = ",")
      SmallOut [1,6] <- paste(unique(as.vector(unlist(strsplit(as.character(tmp$Ontology),split=",")))), collapse = ",")
      NewOutComboSmall <- as.data.frame(rbind(NewOutComboSmall,SmallOut))
      }
    }
  }else{
    NewOutComboSmall <- GenesOut
  }
  return(NewOutComboSmall)
  }


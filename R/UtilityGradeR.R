#' Assigns A Grade To Gene Matches Based On Utility
#' @param GenesOut character or numeric vector: unique abstract IDs
#' @param Kingdom character: which kingdom to use
#' @param Add data.frame: user gene names, families, and ontologies to add to library
#' @param Groups data.frame: gene groups, one column output of GeneNamesGroupeR
#' @export
# Written by John M. A. Wojahn August 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University

UtilityGradeR <- function(GenesOut, Kingdom, Add, Groups)
{
  if(class(GenesOut) != "data.frame")
  {
     stop("ERROR: GenesOut must be a data.frame!")
  }
  if(class(Groups) != "data.frame")
  {
     stop("ERROR: Groups must be a data.frame!")
  }
  require(G2PMineR)
  if(Kingdom == "P")
  {
    print("Using Plantae Internal Data")
  }else if(Kingdom == "A"){
    print("Using Animalia Internal Data")
    SwissGenesCombo <- AnimalSwissGenesCombo
  }else if(Kingdom == "F"){
    print("Using Fungi Internal Data")
    SwissGenesCombo <- FungiSwissGenesCombo
  }else{
    stop("ERROR: You need to choose a kingdom!")
  }
  if(!is.null(Add))
  {
    if(ncol(Add) != 3)
    {
      stop("ERROR: If not NULL, Add must be a data.frame of three columns!")
    }
    print("Adding user terms to library")
    Add <- as.data.frame(Add)
    colnames(Add) <- colnames(SwissGenesCombo)
    SwissGenesCombo <- as.data.frame(rbind(SwissGenesCombo, Add))
  }
  Grades <- as.data.frame(matrix(nrow = nrow(GenesOut), ncol = 2))
  colnames(Grades) <- c("Gene","Grade")
  Grades[,1] <- GenesOut$Gene
  pb <- txtProgressBar(min = 1, max = nrow(Grades), style = 3)
  for(i in 1:nrow(GenesOut))
  {
    setTxtProgressBar(pb, i)
    if(GenesOut[i,2] == "No")
    {
      Grades[i,2] <- "F"
    }else if(as.character(GenesOut[i,6]) != "" && !grepl("(one gene)",as.character(GenesOut[i,5]))){
      Grades[i,2] <- "A"
    }else if(as.character(GenesOut[i,6]) != "" && grepl("(one gene)",as.character(GenesOut[i,5]))){
      Grades[i,2] <- "B"
    }else if(as.character(GenesOut[i,6]) == "" && !grepl("(one gene)",as.character(GenesOut[i,5]))){
      Grades[i,2] <- "C"
    }else if(as.character(GenesOut[i,6]) == "" && grepl("(one gene)",as.character(GenesOut[i,5]))){
      Grades[i,2] <- "D"
    }
  }
  print("Remember:")
  print("A = In text AND has ontologies AND has SwissProt family")
  print("B = In text AND has ontologies BUT NOT SwissProt family")
  print("C = In text AND has SwissProt family BUT NOT ontologies")
  print("D = In text AND DOES NOT HAVE SwissProt family AND DOES NOT HAVE ontologies")
  print("F = NOT in text")
  Grades <- as.data.frame(cbind(Grades,GenesOut[,3:6],Groups))
  colnames(Grades)[7] <- "Group"
  return(Grades)
}


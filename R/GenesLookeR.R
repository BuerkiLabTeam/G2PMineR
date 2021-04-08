#' Mines Abstract for Gene Names
#' @param AbstractStrings character vector: of abstracts as strings
#' @param IDs character or numeric vector: unique abstract IDs
#' @param Kingdom character: which kingdom to use
#' @param Add data.frame: user gene names, families, and ontologies to add to library
#' @param SppAbbr character vector: list of species abbreviations to use to make isoforms
#' @export
# Written by John M. A. Wojahn June 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University

GenesLookeR <- function(AbstractStrings,IDs,Kingdom,Add,SppAbbr)
{
  #check classes
  if(class(AbstractStrings) != "character")
  {
    stop("ERROR: AbstractStrings must be a character vector!")
  }
  if(class(IDs) != "character" && class(IDs) != "numeric")
  {
    stop("ERROR: IDs must be a character or numeric vector!")
  }
  if(length(AbstractStrings) != length(IDs))
  {
    stop("ERROR: Your AbstractIDs and IDs must be the same length!")
  }
  options(warn=-1)
  if(exists("SwissGenesCombo"))
  {
    rm(SwissGenesCombo)
  }
  options(warn=0)
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
  if(class(SppAbbr) != "character")
  {
    stop("ERROR: SppAbbr must be character!")
  }
  if(length(SppAbbr) > 0)
  {
    TaxaCombos <- unique(SppAbbr)
  }else{
    "NOITCE: User provided no SppAbbr, not performing isoform analysis!"
  }
  #coerce to right class
  AbstractStrings <- as.vector(AbstractStrings) #coerce to vector
  AbstractIDs <- as.vector(IDs) #coerce to vector
  SwissGenesCombo <- as.data.frame(SwissGenesCombo) #coerce to dataframe
  if("; " %in% SwissGenesCombo$SwissGenes)
  {
    SwissGenesCombo <- SwissGenesCombo[!(SwissGenesCombo$SwissGenes == "; "),]
  }
  SwissGenesComboX <- unique(SwissGenesCombo)
  rm(SwissGenesCombo)
  SwissGenesCombo <- as.data.frame(SwissGenesComboX)
  UniqueGeneNames <- SwissGenesCombo$SwissGenes
  tmp <- gsub("\\,","",as.character(AbstractStrings)) #NOT tolowered!!!
  tmp <- gsub("\\.","",tolower(tmp))
  tmp <- gsub(":","",tolower(tmp))
  tmp <- gsub("<i>","",tolower(tmp))
  tmp <- gsub("</i>","",tolower(tmp))
  tmp <- gsub("<sup>","",tolower(tmp))
  tmp <- gsub("</sup>","",tolower(tmp))
  tmp <- gsub(";","",tolower(tmp))
  bigsstrings <- unique(as.vector(unlist(strsplit((tmp), split = " "))))
  onestring <- paste(bigsstrings, collapse = " ")
  print("...Inferring Hypotheses of Genes In Abstracts")
  UniqueGeneNamesDets <- c(1:length(UniqueGeneNames))
  pb <- txtProgressBar(min = 1, max = length(UniqueGeneNames), style = 3)
  for(i in 1:length(UniqueGeneNames))
  {
    setTxtProgressBar(pb, i)
    UniqueGeneNamesDets[i] <- grepl(tolower(as.character(UniqueGeneNames[i])), tolower(as.character(onestring)))
  }
  SwissGenesComboX <- SwissGenesCombo[which(UniqueGeneNamesDets == 1),]
  rm(SwissGenesCombo)
  SwissGenesCombo <- as.data.frame(SwissGenesComboX)
  SwissGenesIn <- SwissGenesCombo$SwissGenes #make genes only vector
  OutCombo <- as.data.frame(matrix(nrow=length(SwissGenesIn),ncol=6)) #make output dataframe
  colnames(OutCombo) <- c("Gene","InOrNot","Matches","InSitus","Family","Ontology") #columns
  OutCombo$Gene <- SwissGenesIn
  for(i in 1:length(SwissGenesIn)) #for each gene...
  {
    print(sprintf("Processing Gene No. %s of %s (%s)",i,length(SwissGenesIn),as.character(SwissGenesIn[i])))
    tmp <- gsub(",","",as.character((AbstractStrings))) #NOT tolowered!!!
    tmp <- gsub("\\.","",as.character(tmp))
    tmp <- gsub(":","",as.character(tmp))
    tmp <- gsub(";","",as.character(tmp))
    tmp <- gsub("<i>","",as.character(tmp))
    tmp <- gsub("</i>","",as.character(tmp))
    tmp <- gsub("<sup>","",as.character(tmp))
    tmp <- gsub("</sup>","",as.character(tmp))
    splitlist <- strsplit((tmp), split = " ")
    ingenes <- c()
    insitus <- c()
    matches <- c()
    pb <- txtProgressBar(min = 1, max = length(splitlist), style = 3)
    for(k in 1:length(splitlist))
    {
      setTxtProgressBar(pb, k)
      strict <- as.character(splitlist[[k]]) %in% as.character(SwissGenesIn[i])
      if(T %in% strict)
      {
        ingenes <- c(ingenes,as.character(SwissGenesIn[i]))
        insitus <- c(insitus,as.character(SwissGenesIn[i]))
        matches <- c(matches,AbstractIDs[k])
      }
      #now do for putative isoforms
      Isoforms <- paste0(TaxaCombos,as.character(SwissGenesIn[i]))
      for(m in 1:length(Isoforms))
      {
        strict <- as.character(splitlist[[k]]) %in% as.character(Isoforms[m])
        if(T %in% strict)
        {
          ingenes <- c(ingenes,as.character(SwissGenesIn[i]))
          insitus <- c(insitus,as.character(Isoforms[m]))
          matches <- c(matches,AbstractIDs[k])
        }
      }
    }
    ingenes <- unique(ingenes)
    if(length(ingenes) > 0)
    {
      OutCombo[i,2] <- "Yes"
      OutCombo[i,4] <- paste(unique(insitus), collapse = ",")
      OutCombo[i,3] <- paste(unique(matches), collapse = ",")
      OutCombo[i,5] <- as.character(SwissGenesCombo[i,2])
      OutCombo[i,6] <- as.character(SwissGenesCombo[i,3])
    }else{
      OutCombo[i,2] <- "No"
      OutCombo[i,4] <- "No"
      OutCombo[i,3] <- "No"
      OutCombo[i,5] <- "No"
      OutCombo[i,6] <- "No"
    }
  }
  return(OutCombo) #return output to main
}

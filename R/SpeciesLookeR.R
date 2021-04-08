#' Checks If Species is in String
#' @param AbstractStrings character vector: of abstracts as strings
#' @param IDs character or numeric vector: unique abstract IDs
#' @param Kingdom character: which kingdom to use
#' @param Add character vector: user terms to add to phenology library
#' @export
# Written by John M. A. Wojahn June, July 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPLv3 License
# Funded by EPSCoR GEM3 at Boise State University

SpeciesLookeR <- function(AbstractStrings,IDs,Kingdom,Add)
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
  options(warn=-1)
  if(exists("TracheophytaSpp"))
  {
    rm(TracheophytaSpp)
  }
  options(warn=0)
  if(length(AbstractStrings) != length(IDs))
  {
    stop("ERROR: Your AbstractIDs and IDs must be the same length!")
  }
  require(G2PMineR)
  if(Kingdom == "P")
  {
    print("Using Plantae Internal Data")
  }else if(Kingdom == "A"){
    print("Using Animalia Internal Data")
    TracheophytaSpp <- AnimalSpp
  }else if(Kingdom == "F"){
    print("Using Fungi Internal Data")
    TracheophytaSpp <- FungiSpp
  }else{
    stop("ERROR: You need to choose a kingdom!")
  }
  if(!is.null(Add))
  {
    if(class(Add) != "character")
    {
      stop("ERROR: If not NULL, Add must be a character vector!")
    }
    print("Adding user terms to library")
    Add <- as.data.frame(Add)
    colnames(Add) <- colnames(TracheophytaSpp)
    TracheophytaSpp <- as.data.frame(rbind(TracheophytaSpp, Add))
  }
  bigsstrings <- as.vector(AbstractStrings)
  onestring <- paste(bigsstrings, collapse = " ")
  #CropSpp <- as.data.frame(CropSpp)
  AbstractIDs <- as.vector(IDs)
  SppNames <- as.vector(TracheophytaSpp[,1]) #internal object
  GeneraNamesFull <- gsub(" .*","",SppNames)
  UniqueGeneraNames <- unique(gsub(" .*","",SppNames))
  UniqueGeneraNamesDets <- c(1:length(UniqueGeneraNames))
  print("...Inferring Hypotheses of Genera In Abstracts")
  pb <- txtProgressBar(min = 1, max = length(UniqueGeneraNames), style = 3)
  for(i in 1:length(UniqueGeneraNames))
  {
    setTxtProgressBar(pb, i)
    UniqueGeneraNamesDets[i] <- grepl(as.character(UniqueGeneraNames[i]), as.character(onestring))
  }
  InGenera <- sort(UniqueGeneraNames[which(UniqueGeneraNamesDets == 1)])
  tmp <- gsub(",","",tolower(AbstractStrings))
  tmp <- gsub("\\.","",tolower(tmp))
  tmp <- gsub(":","",tolower(tmp))
  tmp <- gsub(";","",tolower(tmp))
  tmp <- gsub("\\(","",tolower(tmp))
  tmp <- gsub("\\)","",tolower(tmp))
  splitlist <- strsplit((tmp), split = " ")
  cleanbigstrings <- tmp
  grandgenlooks <- as.data.frame(matrix(nrow=0,ncol = 4)) #make outmatrix
  colnames(grandgenlooks) <- c("Genus","Species","Matches","Crops")
  print("...Searching Abstracts for Species")
  pb <- txtProgressBar(min = 1, max = length(InGenera), style = 3)
  for(i in 1:length(InGenera))
  {
    setTxtProgressBar(pb, i)
    ThisGenusSpp <- SppNames[which(as.character(GeneraNamesFull) %in% as.character(InGenera[i]))]
    inspp <- c()
    inmatch <- c()
    perfectmatches <- c()
    for(j in 1:length(ThisGenusSpp))
    {
      #First do exact matching for all species
      dets <- c(1:length(splitlist))
      for(k in 1:length(splitlist))
      {
        #print(paste0("...",k))
        dets[k]  <- grepl(tolower(as.character(ThisGenusSpp[j])), tolower(cleanbigstrings[k]))
      }
      if(1 %in% dets)
      {
        perfectmatches <- unique(c(perfectmatches, AbstractIDs[which(dets == T)]))
        grandgenlooksSmall <- as.data.frame(matrix(nrow=1,ncol = 4)) #make outmatrix
        colnames(grandgenlooksSmall) <- c("Genus","Species","Matches","Crops")
        grandgenlooksSmall[1,1] <- gsub(" .*","",ThisGenusSpp[j])
        grandgenlooksSmall[1,2] <- ThisGenusSpp[j]
        grandgenlooksSmall[1,3] <- paste(AbstractIDs[which(dets == T)], collapse = ",")
        grandgenlooksSmall[1,4] <- "No"
        grandgenlooks <- as.data.frame(rbind(grandgenlooks,grandgenlooksSmall))
      }
    }
    generaldets <- c(1:length(splitlist))
    for(k in 1:length(splitlist))
    {
      #print(paste0("...",k))
      generaldets[k]  <- grepl(tolower(as.character(InGenera[i])), tolower(cleanbigstrings[[k]]))
    }
    GeneralIDs <- AbstractIDs[which(generaldets == 1)]
    if(is.null(perfectmatches))
    {
      #see if part of separate word
      realfake <- c(1:length(GeneralIDs))
      for(u in 1:length(GeneralIDs))
      {
        fooby <- which(AbstractIDs %in% GeneralIDs[u])
        realfake[u] <- as.character(InGenera[i]) %in% as.character(splitlist[[fooby]])
      }
      if(1 %in% realfake)
      {
        grandgenlooksSmall <- as.data.frame(matrix(nrow=1,ncol = 4)) #make outmatrix
        colnames(grandgenlooksSmall) <- c("Genus","Species","Matches","Crops")
        grandgenlooksSmall[1,1] <- gsub(" .*","",ThisGenusSpp[j])
        grandgenlooksSmall[1,2] <- "No Species"
        grandgenlooksSmall[1,3] <- paste(GeneralIDs[which(realfake == 1)], collapse = ",")
        grandgenlooksSmall[1,4] <- "No"
        grandgenlooks <- as.data.frame(rbind(grandgenlooks,grandgenlooksSmall))
      }
    }else{
      missingmatches <- GeneralIDs[-which(GeneralIDs %in% perfectmatches)]
      if(length(missingmatches) != 0)
      {
        realfake <- c(1:length(missingmatches))
        for(u in 1:length(missingmatches))
        {
          fooby <- which(AbstractIDs %in% missingmatches[u])
          realfake[u] <- tolower(as.character(InGenera[i])) %in% tolower(as.character(splitlist[[fooby]]))
        }
        if(1 %in% realfake)
        {
          grandgenlooksSmall <- as.data.frame(matrix(nrow=1,ncol = 4)) #make outmatrix
          colnames(grandgenlooksSmall) <- c("Genus","Species","Matches","Crops")
          grandgenlooksSmall[1,1] <- gsub(" .*","",ThisGenusSpp[j])
          grandgenlooksSmall[1,2] <- "No Species"
          grandgenlooksSmall[1,3] <- paste(missingmatches[which(realfake == 1)], collapse = ",")
          grandgenlooksSmall[1,4] <- "No"
          grandgenlooks <- as.data.frame(rbind(grandgenlooks,grandgenlooksSmall))
        }
      }
    }
  }
  grandgenlooks <- as.data.frame(grandgenlooks[,1:3])
  grandgenlooks[which(grandgenlooks$Species == "No Species"),2] <- grandgenlooks[which(grandgenlooks$Species == "No Species"),1]
  return(grandgenlooks)
}

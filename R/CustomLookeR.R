#' Extracts Custom Words From Abstracts
#' @param AbstractStrings character vector: of abstracts as strings
#' @param IDs character or numeric vector: unique abstract IDs
#' @param CustomWords character vector: user terms to use for search library
#' @param CustomCategories character vector: user categories to use for search library
# Written by John M. A. Wojahn December 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU Affero General Public License v. 3
# Funded by EPSCoR GEM3 at Boise State University

CustomLookeR <- function(AbstractStrings, IDs, CustomWords, CustomCategories)
{
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
  if(class(CustomWords) != "character")
  {
    stop("ERROR: CustomWords must be a character vector!")
  }
  if(is.null(CustomCategories) == F)
  {
    if(class(CustomCategories) != "character")
    {
      stop("ERROR: CustomCategories must be a character vector!")
    }
    if(length(CustomCategories) != length(CustomWords))
    {
      stop("ERROR: Your CustomWords and CustomCategories must be the same length!")
    }
  }else{
    print("USER HAS CHOSEN 'NO CATEGORIES' OPTION")
  }
  require(G2PMineR)
  # reduce possible matches
  MoBotTerms <- as.vector(CustomWords)
  AbstractsStrings <- as.vector(AbstractStrings)
  PhenOUT <- as.data.frame(matrix(nrow=length(AbstractsStrings),ncol=5))
  colnames(PhenOUT) <- c("IDs","Strings","PhenoWords","ForwardPhenoArounds","BackwardsPhenoArounds")
  PhenOUT[,1] <- IDs
  PhenOUT[,2] <- AbstractsStrings
  AllAbs <- tolower(paste(AbstractsStrings, collapse = " "))
  GoodMoBotTerms <- c()
  #sieve out the phenotypes that occur in no abstracts
  print("PRE-PROCESSING CUSTOM TERMS")
  pb <- txtProgressBar(min = 1, max = length(MoBotTerms), style = 3)
  for(i in 1:length(MoBotTerms))
  {
    setTxtProgressBar(pb, i)
    if(grepl(MoBotTerms[i], AllAbs))
    {
      GoodMoBotTerms <- c(GoodMoBotTerms,as.character(MoBotTerms[i]))
    }
  }
  #extract phenotypes from abstracts
  print("PROCESSING ABSTRACTS")
  pb <- txtProgressBar(min = 1, max = nrow(PhenOUT), style = 3)
  for(i in 1:nrow(PhenOUT))
  {
    setTxtProgressBar(pb, i)
    inabstermz <- c()
    for(j in 1:length(GoodMoBotTerms))
    {
      if(grepl(" ",GoodMoBotTerms[j]))
      {
        if(grepl(GoodMoBotTerms[j],tolower(as.character(PhenOUT[i,2]))))
        {
          inabstermz <- c(inabstermz, as.character(GoodMoBotTerms[j]))
        }
      }else{
        bloog <- gsub("\\.","",as.character(PhenOUT[i,2]))
        bloog <- gsub(",","",bloog)
        bloog <- gsub(":","",bloog)
        bloog <- gsub(";","",bloog)
        tmp <- as.vector(unlist(strsplit(bloog, split = " ")))
        if(GoodMoBotTerms[j] %in% tmp)
        {
          inabstermz <- c(inabstermz, as.character(GoodMoBotTerms[j]))
        }
      }
    }
    #get around words for each phenotype match in each abstract
    if(is.null(inabstermz))
    {
      PhenOUT[i,3] <- "NO MATCHES"
      PhenOUT[i,4] <- "NO MATCHES"
      PhenOUT[i,5] <- "NO MATCHES"
    }else{
      PhenOUT[i,3] <- paste(inabstermz, collapse = "@")
      for(d in 1:length(inabstermz))
      {
        tmp <- as.vector(unlist(strsplit(as.character(inabstermz[d]), split = "")))
        if(tmp[1] == " ")
        {
          tmp <- tmp[-1]
        }
        if(tmp[length(tmp)] == " ")
        {
          tmp <- tmp[-length(tmp)]
        }
        inabstermz[d] <- as.character(paste(tmp, collapse = ""))
      }
      invisible(capture.output(AWds <- G2PMineR::AroundWordsSearcheR(as.character(PhenOUT[i,2]),as.character(PhenOUT[i,1]),Anchors = inabstermz,n=1,RemoveCommons = T)))
      PhenOUT[i,4] <- AWds$WordsAroundForward <- gsub(",;,",";",AWds$WordsAroundForward)
      PhenOUT[i,5] <- AWds$WordsAroundBackwards <- gsub(",;,",";",AWds$WordsAroundBackwards)
    }
  }
  #infer stats for each around word bigram and choose 3 commonest
  PhenOUT <- PhenOUT[!(PhenOUT$PhenoWords == "NO MATCHES"),]
  matchers <- as.vector(unlist(strsplit(PhenOUT[,3],spli = "@")))
  mtab <- as.data.frame(table(matchers))
  mtab <- as.data.frame(mtab[order(mtab[,2], decreasing = T),])
  coocs <- as.data.frame(matrix(nrow=nrow(mtab),ncol=6))
  colnames(coocs) <- c("PhenoWord","NumberAbs","1stWordPair","2ndWordPair","3rdWordPair","AbsMatches")
  coocs[,1] <- as.character(mtab[,1])
  coocs[,2] <- as.character(mtab[,2])
  for(g in 1:nrow(coocs))
  {
    coocs[g,6] <- paste(PhenOUT[grep(coocs[g,1], PhenOUT[,3]),1], collapse = ",")
  }
  pb <- txtProgressBar(min = 1, max = nrow(coocs), style = 3)
  for(i in 1:nrow(coocs))
  {
    setTxtProgressBar(pb, i)
    PW <- as.vector(PhenOUT[grep(coocs[i,1],PhenOUT[,3]),3])
    whchWrd <- c(1:length(PW))
    whchWrd[1:length(PW)] <- NA
    for(j in 1:length(PW))
    {
      goo <- as.vector(unlist(strsplit(PW[j], split = "@")))
      if((T %in% (goo %in% coocs[i,1])) == F)
      {
        whchWrd[j] <- "BAD"
      }else{
        whchWrd[j] <- paste(which(goo %in% coocs[i,1]), collapse = ",")
      }
    }
    Badz <- grep("BAD",whchWrd)
    if(T %in% grepl("BAD",whchWrd))
    {
      whchWrd <- whchWrd[-Badz]
    }
    FD <- PhenOUT[grep(coocs[i,1],PhenOUT[,3]),4]
    if(T %in% grepl("BAD",whchWrd))
    {
      FD <- FD[-Badz]
    }
    FDSplits <- c()
    for(j in 1:length(FD))
    {
      goo <- as.vector(unlist(strsplit(FD[j], split = "@")))
      floop <- unique(as.vector(unlist(strsplit(goo[as.numeric(as.vector(unlist(strsplit(as.character(whchWrd[j]),split = ","))))], split = ";"))))
      floop <- floop[!(floop == "END_OF_SENTENCE")]
      floop <- floop[!(floop == "START_OF_SENTENCE")]
      FDSplits <- c(FDSplits, floop)
    }
    BD <- PhenOUT[grep(coocs[i,1],PhenOUT[,3]),5]
    if(T %in% grepl("BAD",whchWrd))
    {
      BD <- BD[-Badz]
    }
    BDSplits <- c()
    for(j in 1:length(BD))
    {
      goo <- as.vector(unlist(strsplit(BD[j], split = "@")))
      floop <- as.vector(unlist(strsplit(goo[as.numeric(as.vector(unlist(strsplit(as.character(whchWrd[j]),split = ","))))], split = ";")))
      floop <- floop[!(floop == "END_OF_SENTENCE")]
      floop <- floop[!(floop == "START_OF_SENTENCE")]
      BDSplits <- c(BDSplits, floop)
    }
    FDSplits <- paste0(FDSplits, " ", as.character(coocs[i,1]))
    BDSplits <- paste0(as.character(coocs[i,1]), " ",BDSplits)
    AllSplits <- c(FDSplits,BDSplits)
    Alltab <- as.data.frame(table(AllSplits))
    Alltab <- Alltab[order(Alltab[,2], decreasing = T),]
    coocs[i,3] <- as.character(Alltab[1,1])
    coocs[i,4] <- as.character(Alltab[2,1])
    coocs[i,5] <- as.character(Alltab[3,1])
  }
  if(is.null(CustomCategories) == F)
  {
    key <- data.frame(CustomWords, CustomCategories)
    categories <- as.data.frame(matrix(nrow=nrow(coocs),ncol=1))
    print("ASSIGNING CATEGORIES")
    pb <- txtProgressBar(min = 1, max = nrow(coocs), style = 3)
    for(i in 1:nrow(coocs))
    {
      setTxtProgressBar(pb, i)
      category <- paste(unique(key[which(key$CustomWords == coocs[i,1]),2]),collapse = ";")
      categories[i,1] <- category
    }
    OUT <- as.data.frame(cbind(categories,coocs))
    colnames(OUT) <- c("Category","CustomWord","NumberAbs","1stWordPair","2ndWordPair","3rdWordPair","AbsMatches")
    return(OUT)
  }else{
    return(coocs)
  }
}

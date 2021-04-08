#' Checks If Gene Matches Are English Words
#' @param membership data.frame or matrix: output of AbstractsClusterMakeR
#' @param threshold numerical: between 0 and 1, proportion of abstracts that must contain match to be considered widespread
#' @param singularize boolean: if T then all words are singularized
#' @export
# Written by John M. A. Wojahn June 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPL v. 3 License
# Funded by EPSCoR GEM3 at Boise State University

MembershipInvestigatoR <- function(membership, threshold, singularize)
{
  print("...preparing data")
  membership <- as.data.frame(membership)
  GrOUT <- as.data.frame(matrix(nrow=length(unique(membership[,3])),ncol=6))
  GrOUT[,1] <- c(1:length(unique(membership[,3])))
  colnames(GrOUT) <- c("Group", "NumberNonStopWords", "NumberNonStopWordsOverThreshold", "WordsOverThreshold","WordsOverThresholdAbstractCounts",
                       "NumberWordsUnderThreshold")
  print("...Starting groupwise loops")
  for(i in 1:nrow(GrOUT))
  {
    print(sprintf("Processinf Group No. %s of %s",i,nrow(GrOUT)))
    print("...vectorizing")
    abs <- as.vector(membership[which(membership[,3] == i),2])

    print("...cleaning strings")
    absNoStop <- tolower(abs)
    absNoStopNoPunct <-  gsub("[[:punct:]]","",absNoStop)
    print("...splitting strings")
    #we unique abstract to analyze ONLY interabstract sharing, not intraabstract commonness!
    absNoStopNoPunctNoRep <- c(1:length(absNoStopNoPunct))
    for(j in 1:length(absNoStopNoPunct))
    {
      print(sprintf("...Uniquing Abstract No. %s of %s",j,length(absNoStopNoPunct)))
      print("...Removing common/stopwords")
      tmp <- G2PMineR::CommonWordsRemoveR(as.character(absNoStopNoPunct[j]))
      if(singularize == T)
      {
        print("...singularizing words")
        vecchar <- as.vector(unlist(strsplit(tmp, split = " ")))
        singlevec <- c(1:length(vecchar))
        require(pluralize)
        for(e in 1:length(vecchar))
        {
          if(tolower(vecchar[e]) == "species")
          {
            singlevec[e] <- vecchar[e]
          }else{
            singlevec[e] <- pluralize::singularize(vecchar[e])
          }
        }
        tmp <- as.character(paste(singlevec, collapse = " "))
      }else{
         print("...not singularizing words")
      }
      absNoStopNoPunctNoRep[j] <- paste(unique(as.vector(unlist(strsplit(tmp, split = " ")))),collapse=" ")
    }
    Unlisted <- as.vector(unlist(strsplit(absNoStopNoPunctNoRep,split = " ")))

    print("...calculating loopwise satistics")
    SharedTable <- as.data.frame(table(Unlisted))
    print("...ordering table")
    SharedTable <- as.data.frame(SharedTable[order(SharedTable$Freq,
                                                   decreasing = T),])
    #get rid of pure numbers
    print("...removing pure numbers")
    SharedTable <- as.data.frame(SharedTable[grepl("[A-Za-z]",SharedTable$Unlisted),])
    minamt <- length(absNoStopNoPunct)*threshold
    SharedTableOver <- SharedTable[which(SharedTable$Freq >= minamt),]
    SharedTableUnder <- SharedTable[which(SharedTable$Freq < minamt),]

    print("...populating output")
    GrOUT[i,2] <- length(unique(Unlisted))
    GrOUT[i,3] <- nrow(SharedTableOver)
    GrOUT[i,4] <- paste(SharedTableOver[,1], collapse = "@")
    GrOUT[i,5] <- paste(SharedTableOver[,2], collapse = "@")
    GrOUT[i,6] <- nrow(SharedTableUnder)
  }
  #now do unshared
  return(GrOUT)
}





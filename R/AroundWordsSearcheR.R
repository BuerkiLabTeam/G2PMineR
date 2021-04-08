#' N-gram abstract word search function
#' @param AbstractStrings character vector: of abstracts as strings
#' @param IDs character or numeric vector: unique abstract IDs
#' @param Anchors haracter vector: a vector of words for which to search
#' @param n numeric: the number of words before and after anchors to search
#' @param RemoveCommons boolean: remove common English words
#' @export
# Written by John M. A. Wojahn June 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU Affero General Public License v. 3
# Funded by EPSCoR GEM3 at Boise State University

AroundWordsSearcheR <- function(AbstractStrings,IDs,Anchors,n,RemoveCommons)
{
  Abstracts <- as.vector(AbstractStrings) #vectorize
  Anchors <- as.vector(Anchors) #vectorize
  Anchors <- tolower(Anchors) #tolower
  AbstractsIDs <- as.vector(IDs) #vectorize
  n <- as.numeric(n) #numericalize
  OUT <- as.data.frame(matrix(nrow=length(Abstracts),ncol=5)) #make output
  colnames(OUT) <- c("AbstractID","AbstractText","MatchingAnchors","WordsAroundForward","WordsAroundBackwards")
  OUT[,1] <- AbstractsIDs #populate IDs
  OUT[,2] <- Abstracts #populate text
  #process dashes so that anchorwordtypedashes are preserved
  if(T %in% (grepl("-",Anchors)))
  {
    splittedz <- as.vector(unlist(strsplit(Abstracts,split = " ")))
    dashanchors <- Anchors[grep("-",Anchors)]
    for(i in 1:length(splittedz))
    {
      deter <- c(1:length(dashanchors))
      for(j in 1:length(dashanchors))
      {
        if(grepl(tolower(dashanchors[j]),tolower(splittedz[i])))
        {
          deter[j] <- "in"
        }else{
          deter[j] <- "out"
        }
      }
      if(("in" %in% deter) == F)
      {
        splittedz[i] <- gsub("-"," ",splittedz[i])
      }
    }
    Abstracts <- paste(splittedz, collapse = " ")
  }else{
    Abstracts <- gsub("-"," ",Abstracts)
  }
  Abstracts <- G2PMineR::CommonWordsRemoveR(Abstracts)
	for(i in 1:length(Abstracts)) #for each abstract
	{
	  print(sprintf("Processing Abstract No. %s of %s",i,length(Abstracts)))
	  AbstractText <- Abstracts[i] #get text
	  AbstractText <- as.character(AbstractText) #characterize
    AbstractTextNoPunct <- gsub("[[:punct:]]","",AbstractText) #remove punct.
    #remove double spaces if exists
    AbstractTextNoPunctOneSpace <- gsub("  "," ",AbstractTextNoPunct)
    ReadyString <- tolower(AbstractTextNoPunctOneSpace) #to lower
    for(j in 1:length(Anchors))
    {
      if(grepl(" ",Anchors[j]))
      {
        ReadyString <- gsub(Anchors[j],gsub(" ","_",Anchors[j]),ReadyString)
        Anchors[j] <- gsub(" ","_",Anchors[j])
      }
    }
	  splits <- as.vector(unlist(strsplit(ReadyString,split=" "))) #split
	  uniquesplits <- unique(splits) #unique
	  AnchorBools <- c(1:length(Anchors))
	  AnchorBools[1:length(Anchors)] <- T
	  #AnchorBools <- c(Anchors %in% uniquesplits) #see if any anchors are in
	  if(T %in% AnchorBools) #if at least one is in
	  {
	    GoodAnchors <- Anchors[which(AnchorBools == T)] #for those that are in
	    GRANDNEGOUT <- c() #bigoutputvec
	    GRANDPOSOUT <- c() #bigoutputvec
	    OUT[i,3] <- paste(GoodAnchors, collapse = "@") #collapse with @
	    for(j in 1:length(GoodAnchors)) #for each anchor in abstract
	    {
	      print(sprintf(".....Processing Good Anchor No. %s of %s",j,length(GoodAnchors)))
	      AnchorMatches <- grep(gsub("_"," ",GoodAnchors[j]),splits) #get matches
	      if(length(AnchorMatches) != 0)
	      {
	        NEGOUT <- c() #littleoutputvec
	        POSOUT <- c() #littleoutputvec
          for(k in 1:length(AnchorMatches)) #for each match for this anchor
          {
            #look after
            if((AnchorMatches[k] - 3) >= 1) #if in abstract
            {
              #get matches around
              smallsplits <- splits[(AnchorMatches[k] - n):(AnchorMatches[k])]
              if(length(smallsplits) > 1) #if not end
              {
                smallsplits <- smallsplits[-length(smallsplits)]
                if(length(NEGOUT) > 0)
                {
                  NEGOUT <- c(NEGOUT,";",as.vector(smallsplits))
                }else{
                  NEGOUT <- c(NEGOUT,as.vector(smallsplits))
                }
              }else{
                if(length(NEGOUT) > 0)
                {
                  NEGOUT <- c(NEGOUT,";","END_OF_SENTENCE")
                }else{
                  NEGOUT <- c(NEGOUT,"END_OF_SENTENCE")
                }
              }
            }else{ #if after abstract
              #get matches around
              smallsplits <- splits[1:(AnchorMatches[k])]
              if(length(smallsplits) > 1) #if not end
              {
                smallsplits <- smallsplits[-length(smallsplits)]
                if(length(NEGOUT) > 0)
                {
                  NEGOUT <- c(NEGOUT,";",as.vector(smallsplits))
                }else{
                  NEGOUT <- c(NEGOUT,as.vector(smallsplits))
                }
              }else{
                if(length(NEGOUT) > 0)
                {
                  NEGOUT <- c(NEGOUT,";","END_OF_SENTENCE")
                }else{
                  NEGOUT <- c(NEGOUT,"END_OF_SENTENCE")
                }
              }
            }
            #look before...
            if((AnchorMatches[k] + 3) <= length(splits)) #if in abstract
            {
              #get matches around
              smallsplits <- splits[(AnchorMatches[k]):(AnchorMatches[k] + n)]
              if(length(smallsplits) > 1) #if not start
              {
                smallsplits <- smallsplits[-1]
                if(length(POSOUT) > 0)
                {
                  POSOUT <- c(POSOUT,";",as.vector(smallsplits))
                }else{
                  POSOUT <- c(POSOUT,as.vector(smallsplits))
                }
              }else{
                if(length(POSOUT) > 0)
                {
                  POSOUT <- c(POSOUT,";","START_OF_SENTENCE")
                }else{
                  POSOUT <- c(POSOUT,"START_OF_SENTENCE")
                }
              }
            }else{ #if before abstract
              #get matches around
              smallsplits <- splits[(AnchorMatches[k]):length(splits)]
              if(length(smallsplits) > 1) #if not start
              {
                smallsplits <- smallsplits[-1]
                if(length(POSOUT) > 0)
                {
                  POSOUT <- c(POSOUT,";",as.vector(smallsplits))
                }else{
                  POSOUT <- c(POSOUT,as.vector(smallsplits))
                }
              }else{
                if(length(POSOUT) > 0)
                {
                  POSOUT <- c(POSOUT,";","START_OF_SENTENCE")
                }else{
                  POSOUT <- c(POSOUT,"START_OF_SENTENCE")
                }
              }
            }
          } #anchormatches
	          NEGOUT <- paste(NEGOUT,collapse = ",") #collapse all same anchortype
            POSOUT <- paste(POSOUT,collapse = ",") #collapse all same anchortype
	          GRANDNEGOUT <- c(GRANDNEGOUT,NEGOUT) #add to multianchor vector
	          GRANDPOSOUT <- c(GRANDPOSOUT,POSOUT) #add to multianchor vector
	      }else{
	        tmp <- as.vector(unlist(strsplit(OUT[i,3], "@")))
	        tmpsml <- tmp[-which(tmp == GoodAnchors[j])]
	        OUT[i,3] <- paste(tmpsml, collapse = "@")
	      }
	      OUT[,4] <- paste(GRANDNEGOUT, collapse = "@") #collapse by anchortype
	      OUT[,5] <- paste(GRANDPOSOUT, collapse = "@") #collapse by anchortype
	    }
	  }else{
      OUT[,3] <- "NO_ANCHOR_MATCHES" #sad!
      OUT[,4] <- "NO_ANCHOR_MATCHES" #sad!
      OUT[,5] <- "NO_ANCHOR_MATCHES" #sad!
	  }
	}
  print("Remember Intra-Anchor Matches Are Sep By , And Inter-Anchor Matches Are Sep By @")
  return(OUT) #return output
}

#' Checks If Gene Matches Are English Words
#' @param String character: a string
#' @export
# Written by John M. A. Wojahn June 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU AGPL v. 3 License
# Funded by EPSCoR GEM3 at Boise State University

CommonWordsRemoveR<- function(String)
{
    #to make sure it is lowercase with no punctuation
    WikiWords <-c("the","be", "to", "of", "and", "a",
                  "in", "that", "have", "I", "it", "for",
                  "not", "on", "with","he", "as", "you",
                  "do", "at", "this", "but", "his", "by", "from",
                  "they", "we", "say", "her", "she", "or",
                  "an", "will", "my", "one","all", "would",
                  "there", "their", "what", "so", "up", "out", "if",
                  "about", "who", "get", "which",
                  "go", "me", "when", "make", "can",
                  "like", "time", "no", "just", "him",
                  "know", "take", "people","into", "year", "your",
                  "good", "some", "could", "them", "see",
                  "other", "than", "then", "now", "look",
                  "only", "come", "its", "over", "think", "also",
                  "back", "after", "use", "two", "how",
                  "our", "work", "first", "well", "way",
                  "even", "new", "want", "because", "any",
                  "these", "give", "day", "most", "us",
                 "am","are","is","be","were","was","being","been","became",
                    "do","have","has", "did","do","doing","done")
    CommonWords <- as.vector(WikiWords)
    SplitString <- as.vector(unlist(strsplit(String,split = " ")))
    UncommonStringVector <- as.vector(matrix(nrow = length(SplitString), ncol = 1))
    for(i in 1:length(SplitString))
    {
      if(tolower(gsub("[[:punct:]]","",as.character(SplitString[i]))) %in% CommonWords)
      {
        if(grepl("[[:punct:]]",as.character(SplitString[i])))
        {
          SplitWord <- as.vector(unlist(strsplit(SplitString[i],split = "")))
          WordPunct <- SplitWord[grep("[[:punct:]]",SplitWord)]
          if(T %in% (WordPunct %in% c(".",",",":",";","!","?")))
          {
            NumBack <- 1
            while((UncommonStringVector[i-NumBack] == "" || is.na(UncommonStringVector[i-NumBack])) && (i-NumBack) != 1)
            {
              NumBack <- NumBack + 1
            }
            UncommonStringVector[i-NumBack] <- as.character(paste0(UncommonStringVector[i-NumBack],paste(WordPunct, collapse = "")))
          }else{
            UncommonStringVector[i] <- as.character(paste(WordPunct,collapse = ""))
          }
        }else{
          UncommonStringVector[i] <- ""
        }
      }else{
        UncommonStringVector[i] <- as.character(SplitString[i])
      }
    }
    UncommonStringVector <- UncommonStringVector[!is.na(UncommonStringVector)]
    UncommonStringVector <- UncommonStringVector[UncommonStringVector != ""]
    UncommonString <- paste(UncommonStringVector, collapse = " ")
    return(UncommonString)
}

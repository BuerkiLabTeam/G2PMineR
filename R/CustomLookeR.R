GenuswiseLookeR <- function(genz,use_partial,key)
{
  if(ncol(genz) != 3)
  {
    stop("ERROR: genz must have three columns, first is family, second is genus, third is number of species per genus")
  }
  if(typeof(use_partial) != "logical")
  {
    stop("ERROR: use_partial must be either T or F")
  }
  if(nrow(genz) == 0)
  {
    stop("ERROR: genz must have at least one row")
  }
  rentrez::set_entrez_key(as.character(key))
  OUT <- as.data.frame(matrix(nrow=0,ncol=9))
  colnames(OUT) <- c("Family","Genus","PubMedID","GenBankIDs","NumTaxa","ListOfTaxa","Dated","Abstract","FullPaperLinks")
  filez <- list.files()
  if("OUT_partial.csv" %in% filez)
  {
    if(use_partial == T)
    {
      message("Using partial results")
      OUT <- read.csv("OUT_partial.csv")
    }
  }
  if(use_partial == T)
  {
    if(length(unique(OUT$Genus))+1 < nrow(genz))
    {
      numzz <- length(unique(OUT$Genus))+1:nrow(genz)
    }else{
      stop("Already done! If you want to redo, please set use_partial = F")
    }
  }else{
    numzz <- 1:nrow(genz)
  }
  for(i in numzz)
  {
    message(sprintf("Processing Genus No. %s of %s",i,nrow(genz)))
    NumGenBankTot <- rentrez::entrez_search("nucleotide",sprintf("%s [Organism]",genz[i,2]),limit=10000,retmax=10000)
    #print number genbank seqs
    if(NumGenBankTot[["count"]] > 0.25*genz[i,3])
    {
      PubMedIDz <- rentrez::entrez_search("PubMed",sprintf("phylogeny AND phylogenetics AND %s",genz[i,2]),limit=10000)
      PubMedIDz <- PubMedIDz$ids
      tobind <- as.data.frame(matrix(nrow=length(PubMedIDz),ncol=9))
      colnames(tobind) <- c("Family","Genus","PubMedID","GenBankIDs","NumTaxa","ListOfTaxa","Dated","Abstract","FullPaperLinks")
      tobind$Family <- genz[i,1]
      tobind$Genus <- genz[i,2]
      tobind$PubMedID <- PubMedIDz
      #AbstractsStrings <- G2PMineR::AbstractsGetteR(PubMedIDz)
      #tobind$Abstract <- AbstractsStrings
      GenBankIDzAll <- c()
      for(j in 1:length(PubMedIDz))
      {
        message(sprintf("Processing PubMed No. %s of %s",j,length(PubMedIDz)))
        AbstractsStrings <- G2PMineR::AbstractsGetteR(PubMedIDz[j])
        tobind[j,8] <- AbstractsStrings
        paper_links <- paste(unlist(rentrez::linkout_urls(rentrez::entrez_link(dbfrom="pubmed", id=PubMedIDz[j], cmd="llinks"))),collapse="  ;  ")
        tobind[j,9] <- paper_links
        them <- RISmed::EUtilsGet(as.numeric(PubMedIDz[j]),type="efetch",db="pubmed")
        title <- them@ArticleTitle
        GenBankIDz <- rentrez::entrez_search("nucleotide",title,limit=10000)
        GenBankIDz <- GenBankIDz[["ids"]]
        tobind[j,4] <- paste(GenBankIDz,collapse=",")
        GenBankIDzAll <- c(GenBankIDzAll,GenBankIDz)
        CustomWords <- c("dated","BEAST","fossil","calibrated","time-calibrated","diversification","trait-diversification")
        CustomCategories <- NULL
        minedz <- G2PMineR::CustomLookeR(AbstractsStrings,1,CustomWords,CustomCategories)
        if(is.null(minedz))
        {
          print("Likely Not Dated")
          tobind[j,7] <- "No"
        }else{
          print("Likely Dated")
          tobind[j,7] <- paste(unique(minedz$PhenoWord),collapse=",")
        }
      }
      if(length(GenBankIDzAll) > 0)
      {
        for(q in 1:length(PubMedIDz))
        {
          taxa <- c()
          message(sprintf("Processing PubMed No. %s of %s",q,length(PubMedIDz)))
          GenBankIDzSplitz <- unlist(strsplit(tobind[q,4],split=","))
          if(length(GenBankIDzSplitz) == 0)
          {
            print("No GenBankIDs Found, moving on...")
            tobind[q,5] <- 0
            tobind[q,6] <- "None"
          }else{
            message(sprintf("Found %s GenBankIDs",length(GenBankIDzSplitz)))
            for(l in 1:length(GenBankIDzSplitz))
            {
              Sys.sleep(0.34)
              message(sprintf("Processing PubMed No. %s of %s; GenBank No. %s of %s",q,length(PubMedIDz),l,length(GenBankIDzSplitz)))
              thisacc <- NULL
              countr <- 1
              while(is.null(thisacc))
              {
                message(sprintf("EUtils atttempt No. %s of 5",countr))
                countr <- countr + 1
                if(countr == 5){break}
                tryCatch({
                  thisacc <- RISmed::EUtilsGet(GenBankIDzSplitz[l],db="nucleotide")
                }, error=function(e){})
              }
              if(!is.null(thisacc))
              {
                message("EUtils suceeded!")
                taxon <- thisacc[["GBSeq"]][["GBSeq_organism"]][[1]]
                taxa <- c(taxa,taxon)
              }else{
                message("EUtilsGet failed 5 times, may be corrupted, moving on...")
              }
            }
            taxa <- taxa[grepl(genz[i,2],taxa)]
            if(length(taxa) == 0)
            {
              tobind[q,5] <- 0
              tobind[q,6] <- "None"
            }else{
              tobind[q,5] <- length(unique(taxa))
              tobind[q,6] <- paste(unique(taxa),collapse=",")
            }
          }
        }
      }else{
        tobind[,5] <- "None"
        tobind[,6] <- "None"
      }
    }else{
      tobind <- as.data.frame(matrix(nrow=length(PubMedIDz),ncol=9))
      colnames(tobind) <- c("Family","Genus","PubMedID","GenBankIDs","NumTaxa","ListOfTaxa","Dated","Abstract","FullPaperLinks")
      tobind[1,] <- "Fewer than 25% Genus"
      tobind[1,1] <- genz[i,1]
      tobind[1,2] <- genz[i,2]
    }
    OUT <- rbind(OUT,tobind)
    write.csv(OUT,"OUT_partial.csv",row.names=F)
    if(i == nrow(genz))
    {
      write.csv(OUT,"OUT_final.csv",row.names=F)
    }
  }
}

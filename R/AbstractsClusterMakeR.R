#' Make Preliminary Network for QC
#' @param AbstractStrings character vector: of abstracts as strings
#' @param IDs character or numeric vector: unique abstract IDs
#' @export
# Written by John M. A. Wojahn July 2020
# This is Free and Open-Source Software (F.O.S.S.)
# © J.M.A. Wojahn, S.J. Galla, A.E. Melton, S. Buerki
# Provided under the GNU Affero General Public License v. 3
# Funded by EPSCoR GEM3 at Boise State University
AbstractsClusterMakeR <- function(AbstractStrings, IDs)
{
  require(text2vec)
  AbstractsStrings <- AbstractStrings
  IDs <- as.vector(IDs)
  # tokenize the abstract strings
  print("...tokenizing the abstract strings")
  tokens <- text2vec::word_tokenizer(AbstractsStrings)
  # Get rid of stopwords
  print("...getting rid of stopwords")
  v = create_vocabulary(itoken(tokens))
  v = prune_vocabulary(v, term_count_min = 5, doc_proportion_max = 0.5)
  it = itoken(tokens)
  vectorizer = vocab_vectorizer(v)
  # make the model
  print("...making model")
  dtm = create_dtm(it, vectorizer)
  tcm = create_tcm(it, vectorizer, skip_grams_window = 5)
  glove_model = GloVe$new(rank = 50, x_max = 10)
  # train the model
  wv = glove_model$fit_transform(tcm, n_iter = 5)
  # get average of main and context vectors as proposed in GloVe paper
  wv = wv + t(glove_model$components)
  # run the model
  print("...running relaxed word movers distance model")
  rwmd_model = RelaxedWordMoversDistance$new(dtm, wv)
  # get distance matrix
  print("...calculating distance matrix")
  rwms = rwmd_model$sim2(dtm)
  # turn into actual distance matrix
  print("...getting distance matrix")
  rwms <- as.dist(rwms, diag = FALSE, upper = FALSE)
  # make the network
  print("...making the network")
  net <- qgraph::qgraph(rwms, layout ="spring", labels = colnames(rwms),threshold = 0,DoNotPlot=T)
  #convert to igraph object
  print("...converting to igraph object")
  inet <- igraph::as.igraph(net)
  # infer group members
  # 4 is based on empirical data
  print("...performing cluster walktrap")
  cw <- igraph::cluster_walktrap(inet,steps = 4)
  print("...inferring membership")
  members <- igraph::membership(cw)
  membership <- data.frame(as.data.frame(IDs),as.data.frame(AbstractsStrings),as.vector(members))
  #network <- igraph::plot.igraph(inet,mark.groups = members)
  colnames(membership) <- c("IDs","Strings","Membership")
  print("...creating outlist")
  OUTLIST <- list(membership,rwms,net)
  names(OUTLIST) <- c("Membership","DistanceMatrix","Network")
  return(OUTLIST)
}

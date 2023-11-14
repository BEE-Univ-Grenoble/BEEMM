#' markov_likelihood
#'
#' @param sequence the sequence tested
#' @param table_kmer_proba table which contain kmer, count, proba of the model
#' @param log
#'
#' @return
#' @export
#'
#' @examples

markov_likelihood<-function(sequence,table_kmer_proba,log=TRUE){

  ws<-kmer_size(table_kmer_proba)
  table<-kmer_table(sequence,ws)



  proba <- left_join(table,table_kmer_proba,by="kmer") %>%
             mutate(p=count.x * lprob) %>%
             summarise(p=sum(p)) %>%
             pull(p)



  Init<-str_to_lower(substr(sequence,1,ws-1))

  Initlprob <- table_kmer_proba$knowing_lprob[table_kmer_proba$knowing==Init]

  proba <- proba + Initlprob[1]


  if(log==FALSE){
    proba <- exp(proba)
  }
  return(proba)
}

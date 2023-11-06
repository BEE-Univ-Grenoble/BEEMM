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

  ordre<-kmer_size(table_kmer_proba)
  table<-kmer_table(sequence,ordre)

  Tableau <- left_join(table,table_kmer_proba,by="kmer")%>%
             mutate(p=count.x * lprob) %>%
             summarise(p=sum(p))


  Init<-substr(sequence,1,ordre-1)

  Initlprob<-table_kmer_proba$knowing_lprob[table_kmer_proba$knowing==Init,]

  p<-p+Initlprob

  if(log==FALSE){
    p<-exp(p)
  }
  return(p)
}

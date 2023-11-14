#' A function to return conditional probabilities from a Markov model of a given order.
#'
#' @param fasta_file A fasta file, containing a list of genetic sequences.
#' @param order_markov Order of the Markov model.
#'
#' @return A data frame with the following columns:
#'   \describe{
#'     \item{kmer}{The kmer sequence.}
#'     \item{knowing}{The prefix of the kmer (substring from the beginning to k-1).}
#'     \item{emission}{The last character of the kmer (substring from k to the end).}
#'     \item{count}{The count of occurrences of the kmer.}
#'     \item{knowing_prob}{The probability of the prefix (knowing) occurring.}
#'     \item{knowing_lprob}{The logarithm of the knowing_prob.}
#'     \item{prob}{The conditional probability of the kmer given the prefix.}
#'     \item{lprob}{The logarithm of the prob.}
#'   }
#'
#' @examples
#'   #To build a Markov model of order 3, on a fasta file named "mus_cpg_app.fa" containing sequences.
#'   proba_markov <- build_markov("mus_cpg_app.fa", 3)
#'
#' @export

build_markov <- function(sequences, order_markov) {
  k = order_markov + 1
  kmer_table_seq(sequences, k) %>%
    proba <- kmer_proba(kmer_tab)
}




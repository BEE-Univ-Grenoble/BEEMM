#' Function to calculate the conditional probability of each emitted nucleotide
#'
#' This function calculates the probability of each kmer sequence based on its prefix and emission.
#' It takes a vector of kmers and an optional pseudocount value as input.
#'
#' @param kmers A vector of kmers.
#' @param pseudocount The pseudocount value to be used (default = 1).
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
#' kmers <- as_kmer_table(data.frame(kmer = c("AAA", "AAC", "AAG"), count = c(10, 20, 30)))
#' kmer_proba(kmers)
#'
#' @export
#' @importFrom stringr str_sub
#' @importFrom dplyr mutate group_by ungroup select
kmer_proba <- function(kmers, pseudocount = 1) {
    k <- kmer_size(kmers)

    probs <-  sum_table(
        kmers,
        every_kmer_table(k,
            pseudocount = pseudocount
        )
    ) %>%
        mutate(
            knowing = str_sub(kmer, 1, k - 1),
            emission = str_sub(kmer, k, -1),
        ) %>%
        group_by(knowing) %>%
        mutate(
            knowing_prob = sum(count),
            prob = count / knowing_prob,
            prob = ifelse(is.na(prob), 0, prob),
            lprob = log(prob)
        ) %>%
        ungroup() %>%
        mutate(
            knowing_prob = knowing_prob / sum(count),
            knowing_prob = ifelse(is.na(knowing_prob), 0, knowing_prob),
            knowing_lprob = log(knowing_prob)
        ) %>%
        select(kmer, knowing, emission, count, knowing_prob, knowing_lprob, prob, lprob)

    class(probs) <- class(kmers)
    class(probs)[1] <- "beemm_proba_table"
    attr(probs, "kmer_size") <- kmer_size(kmers)

    probs
}

#' Check if an object is a kmer probability table
#'
#' @param x An object to check
#' @return Logical value indicating if the object is a kmer probability table
#' 
#' @export 
is_kmer_proba <- function(x) {
  "beemm_proba_table" %in% class(x)
}


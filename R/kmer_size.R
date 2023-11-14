#' Get the kmer size of a kmer table
#'
#' @param x A kmer table object
#'
#' @return The kmer size of the kmer table
#' @export
kmer_size <- function(x) {
  # Check if x is of class beemm_kmer_table
  if (!(is_kmer_table(x) || is_kmer_proba(x))) {
    stop("kmer_size requires an x value of class beemm_kmer_table")
  }

  # Retrieve the kmer size attribute from x
  attr(x, "kmer_size")
}

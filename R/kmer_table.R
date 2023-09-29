#' Convert a data frame into a k-mer table
#'
#' This function takes a data frame and converts it into a k-mer table. The input data frame must have columns named "kmer" and "count".
#'
#' @param x The input data frame.
#'
#' @return A k-mer table.
#'
#' @export
#'
#' @examples
#' data <- data.frame(kmer = c("A", "T", "C"), count = c(2, 5, 1))
#' as_kmer_table(data)
#'
#' @importFrom dplyr select filter
#' @importFrom tibble as_tibble
as_kmer_table <- function(x) {
  if (!all(c("kmer", "count") %in% names(x))) {
    stop("A beemm_kmer_table must have columns named : kmer & count")
  }

  kmers <- x %>%
    select(kmer, count) %>%
    filter(str_detect(kmer, pattern = "^[acgt]+$")) %>%
    as_tibble()

  k <- unique(nchar(kmers$kmer))

  if (length(k) > 1) {
    stop(glue(paste(
      "All the kmers do not have the same size",
      "({paste(k, collapse = ', ')})"
    )))
  }

  class(kmers) <- c("beemm_kmer_table", class(kmers))
  attr(kmers, "kmer_size") <- k
  attr(kmers, "missing_kmers") <- 4^k - nrow(kmers)

  kmers
}

#' Function to create a table with all possible k-mers and their counts
#'
#' @param k The length of the k-mers
#' @param pseudocount The initial count for each k-mer (default = 1)
#'
#' @return A k-mer table with k-mers and their counts
#'
#' @export
#'
#' @examples
#' every_kmer_table(3)
every_kmer_table <- function(k, pseudocount = 1) {
  nucleotides <- c("a", "c", "g", "t")
  kmers <- expand.grid(replicate(k, nucleotides, simplify = FALSE))
  kmers <- apply(kmers, 1, paste, collapse = "")
  as_kmer_table(tibble(
    kmer = kmers,
    count = rep(pseudocount, length(kmers))
  ))
}

#' Check if an object is a kmer table
#'
#' @param x An object to check
#' @return Logical value indicating if the object is a kmer table
is_kmer_table <- function(x) {
  "beemm_kmer_table" %in% class(x)
}


#' Generate k-mer table from a given sequence
#'
#' This function takes a DNA sequence and generates a k-mer table,
#' which counts the occurrence of each k-mer in the sequence.
#'
#' @param sequence The DNA sequence
#' @param k The length of the k-mers
#'
#' @return A k-mer table with k-mer counts
#'
#' @examples
#' sequence <- "ATCGATCGATCG"
#' kmer_table(sequence, 3)
#'
#' @export
#'
#' @importFrom tibble tibble
#' @importFrom dplyr mutate group_by summarize
#' @importFrom stringr str_to_lower str_sub
#'
kmer_table <- function(sequence, k) {
  # Get the length of the sequence
  n <- nchar(sequence)

  # Generate all possible k-mers
  kmers <- tibble(start = 1:(n - k + 1)) %>%
    mutate(end = start + k - 1) %>%
    mutate(kmer = str_to_lower(str_sub(sequence, start, end))) %>%
    group_by(kmer) %>%
    summarize(count = n())

  # Convert the k-mers to a k-mer table
  as_kmer_table(kmers)
}



#' Get the kmer size of a kmer table
#'
#' @param x A kmer table object
#'
#' @return The kmer size of the kmer table
#' @export
kmer_size <- function(x) {
  # Check if x is of class beemm_kmer_table
  if (!is_kmer_table(x)) {
    stop("kmer_size requires an x value of class beemm_kmer_table")
  }

  # Retrieve the kmer size attribute from x
  attr(x, "kmer_size")
}



#' Calculate the sum of counts for each kmer in multiple kmer tables
#'
#' @param ... A variable number of beemm_kmer_table objects
#'
#' @return A beemm_kmer_table object with the summed counts for each kmer
#'
#' @details This function takes a variable number of beemm_kmer_table objects
#' and calculates the sum of counts for each kmer across all tables. It checks
#' if all the input arguments are beemm_kmer_table objects and if they have the
#' same kmer size. If not, an error is thrown.
#'
#' @examples
#' # Create example kmer tables
#' table1 <- as_kmer_table(data.frame(kmer = c("AAA", "AAC", "AAG"), count = c(10, 20, 30)))
#' table2 <- as_kmer_table(data.frame(kmer = c("AAA", "AAT"), count = c(5, 15)))
#'
#' # Calculate the sum of counts for each kmer
#' sum_table(table1, table2)
#'
#' @export
#'
sum_table <- function(...) {
  arguments <- list(...) # Convert the arguments to a list

  if (!all(sapply(arguments, is_kmer_table))) {
    stop("At least one argument is not a beemm_kmer_table")
  }

  sizes <- unique(sapply(arguments, kmer_size)) # Get the unique kmer sizes

  if (length(sizes) > 1) {
    stop(glue(paste(
      "All the kmers do not have the same size",
      "({paste(sizes, collapse = ', ')})"
    )))
  }

  kmers <- do.call(bind_rows, arguments) %>% # Combine all kmer tables into one
    group_by(kmer) %>% # Group by kmer
    summarise(count = sum(count)) # Sum the counts for each kmer

  as_kmer_table(kmers) # Convert the result to a beemm_kmer_table object
}


#' Generate k-mer table for a list of sequences
#'
#' This function takes a list of sequences and a value of k as input and returns a table of k-mer frequencies for each sequence combined into a single table.
#' The function expects a data frame with a column named 'sequence' that contains the sequences.
#'
#' @param sequences A data frame with a column named 'sequence' containing the sequences
#' @param k The value of k for generating k-mer table
#'
#' @return A table of k-mer frequencies for each sequence combined into a single table
#'
#' @examples
#' sequences <- data.frame(sequence = c("ACGT", "CGTA"))
#' kmer_table_seq(sequences, 2)
#'
#' @export
kmer_table_seq <- function(sequences, k) {
  sequences$sequence %>%
    lapply(kmer_table, k = k) %>% # Generate k-mer table for each sequence
    do.call(sum_table, .) # Combine k-mer tables into a single table
}

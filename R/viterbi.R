#' HMM Viterbi segmentation of a DNA sequence among two Markov chain models
#'
#' @param sequence a character corresponding to the DNA sequence to segment
#' @param M0 the first Markov model
#' @param M1 the second Markov model
#' @param transitions the transition matrix between models as a 2x2 matrix
#'
#' @return
#' @export
#'
#' @examples
viterbi <- function(sequence, M0, M1, transitions) {
  tm <- log(transitions)

  # Extract order (o) and kmer size (ws) associated to the models
  o <- nchar(M0$knowing[1])
  ws <- nchar(M0$kmer[1])


  # Associate to each kmer its emission probability according M0 & M1

  kmer_vector(sequence, ws) %>%
    left_join(M0, by = "kmer") %>%
    select(kmer, log_M0 = lprob) %>%
    left_join(M1, by = "kmer") %>%
    select(kmer, log_M0, log_M1 = lprob) -> table_proba

  # Scoring

  ## Estimates initial probability according to M0

  p <- M0 %>% filter(kmer == table_proba$kmer[1]) %>%
    pull(knowing_lprob)

  proba_M0 <- p + tm[1, 1]

  ## Estimates initial probability according to M1

  p <- M1 %>% filter(kmer == table_proba$kmer[1]) %>%
    pull(knowing_lprob)

  proba_M1 <- p + tm[2, 2]

  ## Fill the scoring & path matrix

  M0_path <- character(nrow(table_proba))
  M1_path <- character(nrow(table_proba))

  print(proba_M0)
  print(proba_M1)

  for (i in 1:length(M0_path)) {
    P_M0_M0 <- proba_M0 + table_proba$log_M0[i] + tm[1, 1]
    P_M0_M1 <- proba_M1 + table_proba$log_M0[i] + tm[2, 1]

    if (P_M0_M1 > P_M0_M0) {
      proba_M0 <- P_M0_M1
      M0_path[i] <- "M1"
    } else {
      proba_M0 <- P_M0_M0
      M0_path[i] <- "M0"
    }

    P_M1_M1 <- proba_M1 + table_proba$log_M1[i] + tm[2, 2]
    P_M1_M0 <- proba_M0 + table_proba$log_M1[i] + tm[1, 2]

    if (P_M1_M0 > P_M1_M1) {
      proba_M1 <- P_M1_M0
      M1_path[i] <- "M0"
    } else {
      proba_M1 <- P_M1_M1
      M1_path[i] <- "M1"
    }
  }

  # Backtracking

  path <- character(nrow(table_proba))

  if (proba_M0 > proba_M1) {
    current_state <- "M0"
    p <- proba_M0
  } else {
    current_state <- "M1"
    p <- proba_M1
  }

  for (i in length(path):1) {
    path[i] <- current_state
    if (current_state == "M0") {
      current_state <- M0_path[i]
    } else {
      current_state <- M1_path[i]
    }
  }

  unkown <- max(o - 1, 0)
  path <- c(rep("-", unkown), current_state, path)

  kmer_vector(sequence, 1) %>%
    rename(nuc=kmer) %>%
    mutate(position = 1:length(path),
           hidden = path) %>%
    select(position,nuc,hidden) -> result

  cl <- class(result)
  cl[1] <- "beemm_kmer_viterbi"
  class(result) <- cl

  attr(result, "log_likelihood") <- p

  result
}


#' Extracts every segment of a sequence as an individual sequence.
#'
#' @param segmentation : a tibble returned by the `viterbi` function
#'
#' @return a sequence tibble containing four columns
#'   - id : the identifier of the sequence
#'   - model : the model corresponding to that segment
#'   - from : the start position of the fragment in the segmented sequence
#'   - to : the end position of the fragment in the segmented sequence
#'   - sequence : the nucleotide sequence of the fragment
#'
#' @export
#'
#' @md
#' @examples
HMM_as_sequence <- function(segmentation) {
  i <- 0
  map <- list()
  previous <- "-"
  first <- TRUE
  frg <- 1

  for (h in segmentation$hidden) {
    i <- i + 1
    if (h == "-") next

    if (h != previous) {
      if (!first) {
        seq <- paste(segmentation$nuc[start:(i - 1)], collapse = "")
        map <- c(map, list(list(
          id = sprintf("Segment_%d", frg),
          model = previous,
          from = start,
          to = i - 1,
          sequence = seq
        )))
        frg <- frg + 1
      }
      first <- FALSE
      start <- i
      previous <- h
    }
  }

  seq <- paste(segmentation$nuc[start:i], collapse = "")
  map <- c(map, list(list(
    id = sprintf("Segment_%d", frg),
    model = previous,
    from = start,
    to = i,
    sequence = seq
  )))

  map %>%
    lapply(as_tibble) %>%
    do.call(bind_rows, .)
}

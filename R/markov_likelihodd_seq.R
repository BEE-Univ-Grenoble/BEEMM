#' Function to calculate the log likelihood of every sequences
#'
#' This function calculate the log likelihood of multiples sequences, with a chosen Markov model order.
#'
#' @param seqs any table of multiple colones were the sequences are contain in a colone named "sequence"
#' @param ... table of probability of a markov model created with function build_markov, model can be named : model1 = modelprobtable.
#'
#' @return a table of columns
#'   - `...` : all the colomns contain in the original table
#'   - `sequence` : the DNA sequence (character)
#'   - `lprob` : the (log)likelihood for each sequence, if the model name is precised, colone name = model name.
#' @export
#'
#' @importFrom dplyr mutate
#' @return a vector of likelihood values
#'
#' @examples
#' Applying with one markov model order :
#' loglikelihood_table <- markov_likelihood(cpg_sequences, model1 = M1)
#'
#' Applying with multiple markov model order :
#' loglikelihood_table <- markov_likelihood_seq(cpg_sequences, model1 = M1) %>%
#' markov_likelihood_seq(model2 = M2)
#'
markov_likelihood_seq <- function(seqs,...) {
  models <- list(...)
  model <- models[[1]]

  model_names <- names(models)

  if (is.null(model_names)) {
    lprob_name <- "lprob"
    prob_name <- "prob"
  } else {
    lprob_name <- glue("lprob_{model_names[1]}")
    prob_name <- glue("prob_{model_names[1]}")
  }

output <- seqs %>%
  mutate(loglikelihood = sapply(sequence,markov_likelihood,model,
                                log = TRUE,
                                USE.NAMES = FALSE),
         likelihood = exp(loglikelihood)
         ) %>%
  rename({{lprob_name}} := loglikelihood,
         {{prob_name}} := likelihood)

  output %>% select(id,{{prob_name}},{{lprob_name}},sequence)
}

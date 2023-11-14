#' Title
#'
#' @param M0
#' @param M1
#' @param sequence_testing
#' @param priorM0
#'
#' @return
#' @export
#'
#' @importFrom matrixStats logSumExp
#'
#' @examples
markov_bayes <- function (M0, M1, sequence_testing, priorM0 = 0.5) {

  #compute lihelihood for each sequence for both models
  l_M0 <- markov_likelihood_seq(seqs = sequence_testing, M0 = M0)
  l_M1 <- markov_likelihood_seq(seqs = sequence_testing, M1 = M1)

  l_M1 <- l_M1 %>%
    select(lprob_M1, prob_M1 ) #select useful columns

  l_models <- bind_cols(l_M0, l_M1) %>%
    select(id, prob_M0, lprob_M0, prob_M1, lprob_M1, sequence)

  #compute likelihood for each model

  priorM1 = 1 - priorM0

  lpm0 <- log(priorM0)
  lpm1 <- log(priorM1)

  l_models %>%
    mutate(xx_lp0 = lprob_M0 + lpm0,
           xx_lp1 = lprob_M1 + lpm1,
           xx_pseq = mapply(function(x,y) logSumExp(c(x,y)),xx_lp0,xx_lp1),
           lpmodel_M0 = xx_lp0 - xx_pseq,
           lpmodel_M1 = xx_lp1 - xx_pseq,
           pmodel_M0 = exp(lpmodel_M0),
           pmodel_M1 = exp(lpmodel_M1)
    ) %>%
    select(id, prob_M0, lprob_M0,
           pmodel_M0, lpmodel_M0,
           prob_M1, lprob_M1,
           pmodel_M1, lpmodel_M1,
           sequence)
}

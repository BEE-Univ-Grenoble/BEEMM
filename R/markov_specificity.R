.markov_specificity <- function(M0,M1, sequence_testing, priorM1 =0.5) {

  table <- markov_bayes(M0, M1, sequence_testing, priorM1 = priorM1) %>%
    mutate(M1_win = lpmodel_M1 > lpmodel_M0) %>%
    summarise(win = sum(M1_win),n=n(),freq = win/n) %>%
    pull(freq)
}

#' Title
#'
#' @param M0s
#' @param M1s
#' @param sequence_testing_M1
#' @param priorM1
#'
#' @return
#' @export
#'
#' @examples
markov_specificity <- function(M0s,M1s, sequence_testing_M1, priorM1 =0.5) {

  sp <- mapply(.markov_specificity, M0s,M1s,
               list(sequence_testing_M1),
               priorM1 = priorM1)
  sp
}

#' Title
#'
#' @param M0s
#' @param M1s
#' @param sequence_testing_M0
#' @param priorM1
#'
#' @return
#' @export
#'
#' @examples
markov_sensibility <- function(M0s,M1s, sequence_testing_M0, priorM1 =0.5) {

  sp <- mapply(.markov_specificity, M1s,M0s,
               list(sequence_testing_M0),
               priorM1 = 1 - priorM1)
  sp
}

# lapply
# mutate(VN = ifelse(lpmodel_M1 > lpmodel_M0, 1, 0), # nombre sup/nb lignes
#        FP = ifelse(lpmodel_M1 <= lpmodel_M0, 1, 0))


#faire ka sesibilité en premier
#pourquoi inverser les roles de MO et M1




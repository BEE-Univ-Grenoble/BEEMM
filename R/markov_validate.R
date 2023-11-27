#' makov_validate
#'
#' @param sequence_learning_M0 sequences used to build the M0 model
#' @param sequence_learning_M1 sequences used to build the M1 model
#' @param sequence_testing_M0 sequences used to test the M0 model
#' @param sequence_testing_M1 sequences used to test the M1 model
#' @param order_min the smallest order taken to build each model
#' @param order_max the biggest order taken to build each model
#' @param priorM0 the "a priori" probability of the M0 model
#'
#' @return a table with the specificity and sensibility for each order and each model
#' @export
#'
#' @examples

markov_validate<-function(sequence_learning_M0,sequence_learning_M1, sequence_testing_M0,sequence_testing_M1,order_min,order_max, priorM1=0.5){

  # model construction
  orders <- c(order_min:order_max)

    for (i in orders) {
      M0s[i] <- build_markov(sequence_learning_M0,orders[i])
      M1s[i] <- build_markov(sequence_learning_M1,orders[i])
    }

  # specificity and sensibility computation

  specificity <- markov_specificity(M0s, M1s,sequence_testing_M0,priorM1=priorM1)
  sensibility <- markov_sensibility(M0s, M1s,sequence_testing_M1,priorM1=priorM1)

  # dataframe concatenation

  order_comparison <- merge(specificity, sensibility, by="model")
  order_comparison$order<-orders

return(order_comparison)
}

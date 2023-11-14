
#' Markov_validate_kfold
#'
#' Function to validate sensibility and specificity of 2 models with kfold method
#'
#' This function divides randomly sequence data in 2 sub-samples (learning and testing) multiple times
#' and call the markov_validate function each time.
#'
#' @param sequence_M0
#' @param sequence_M1
#' @param order_min
#' @param order_max
#' @param priorM0
#' @param learning_fraction
#' @param nrand
#'
#' @return a concatenate table of markov_validate outputs
#' @export
#'
#' @examples
#' exemple_output <- markov_validate_kfold(sequence_M0,sequence_M1,order_min=1,order_max=8,priorM0 = 0.7,
#' learning_fraction = 0.5, nrand = 30)
#'
#'
markov_validate_kfold <- function(sequence_M0,sequence_M1,
                                  order_min,order_max,
                                  priorM0 = 0.5,
                                  learning_fraction = 0.8, nrand = 10){

  for (i in 1:nrand) {
    nM0 <- length(sequence_M0)
    nM1 <- length(sequence_M1)
    learningM0 <- sample(x= 1:nM0,size=nM0*learning_fraction)
    learningM1 <- sample(x= 1:nM1,size=nM1*learning_fraction)

    seq_learning_M0<- sequence_M0[,learningM0]
    seq_learning_M1<- sequence_M1[,learningM1]

    seq_test_M0<- sequence_M0[,-learningM0]
    seq_test_M1<- sequence_M1[,-learningM1]

    validate <- markov_validate(seq_learning_M0,seq_learning_M1,seq_test_M0,seq_test_M1,order_min,order_max)

    if(i==1){
      output <- validate
    } else {
      output <- rowbind(output,validate)
    }
  }

}


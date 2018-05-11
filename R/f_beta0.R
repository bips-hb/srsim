#' @export
f_beta0 <- function(beta0, margprob, beta1, margprob_parent) { 
  abs(margprob - (1 - margprob_parent)*(exp(beta0) / (1 + exp(beta0))) -
    margprob_parent*(exp(beta0 + beta1) / (1 + exp(beta0 + beta1))))
}
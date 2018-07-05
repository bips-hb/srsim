#' Determining the Intercept
#' 
#' \code{f_beta0} is used by \code{\link{simulateSRS}} to find the appropriate
#' intercept for a logistic regression model. 
#' 
#' @param beta0 The intercept
#' @param margprob The marginal probability 
#' @param beta1 The regression coefficient of the parent variable 
#' @param margprob_parent The marginal probability of the parent 
#' 
#' @return The absolute difference between the desired marginal probability 
#'         and the current marginal probability
#' @export
f_beta0 <- function(beta0, margprob, beta1, margprob_parent) { 
  abs(margprob - (1 - margprob_parent)*(exp(beta0) / (1 + exp(beta0))) -
    margprob_parent*(exp(beta0 + beta1) / (1 + exp(beta0 + beta1))))
}
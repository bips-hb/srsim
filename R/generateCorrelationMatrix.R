#' Correlation Matrix 
#' 
#' Generates a correlation matrix for the SR data set. It randomly selects 
#' a number of drug-event, drug-drug, and event-event pairs to exhibit a 
#' given correlation. It takes the marginal probabilities of the drugs and 
#' events into account, to make sure that the pairs can indeed show the 
#' wanted correlation (see the function \code{\link{validCorrelation}}).
#' \cr\cr 
#' Since the matrix is generated randomly, it might happen that the resulting matrix
#' is not a correlation matrix (i.e., it is not positive definite). In that case, 
#' the process is repeated. One can set the maximal number of iterations with 
#' the parameter \code{max_iter} (the default is set to 1000). 
#' 
#' @param prob_drugs List with the marginal probabilities of the drugs 
#' @param prob_events List with the marginal probabilities of the events
#' @param n_correlated_drugs The number of drug-drug pairs that will have correlation \code{rho_drugs}
#' @param rho_drugs The correlation for the drug-drug pairs 
#' @param n_correlated_events The number of event-event pairs that will have correlation \code{rho_events}
#' @param rho_events The correlation for the event-event pairs 
#' @param n_correlated_pairs The number of drug-event pairs that will have correlation \code{rho}
#' @param rho The correlation for the drug-event pairs 
#' @param max_iter Maximum number of iterations (Default: 1000) 
#' @param verbose Prints the number of iterations (Default: \code{FALSE})
#' 
#' @return \item{corrmat}{The correlation matrix}
#'         \item{valid}{\code{TRUE} when \code{corrmat} is a correlation matrix, \code{FALSE} otherwise} 
#'           
#' 
#' @seealso \code{\link{validCorrelation}}
#' @export
generateCorrelationMatrix <- function(prob_drugs,
                                      prob_events,
                                      n_correlated_drugs,
                                      rho_drugs,
                                      n_correlated_events,
                                      rho_events,
                                      n_correlated_pairs,
                                      rho,
                                      max_iter = 1000, 
                                      verbose = FALSE) {

  # keep on generating correlation matrices until one is valid or 
  # the maximum number of iterations is reached
  for (iter in 1:max_iter) {
    
    if (verbose) {
      cat("Attempt no.", iter, "to generate a valid correlation matrix\n")
    }
    
    mat = generateCorrelationMatrixRcpp(prob_drugs, 
                                        prob_events, 
                                        n_correlated_drugs, 
                                        rho_drugs, 
                                        n_correlated_events, 
                                        rho_events, 
                                        n_correlated_pairs, 
                                        rho) 
    
    # We use the Cholesky decomposition to test for positive definiteness. 
    # Runs much faster for larger matrices than computing all eigenvalues
    if (!is(chol(mat), "error")) {
      return(
        list(
          corrmat = mat, 
          valid = TRUE
        )
      )
    }
  }
  
  # did not find a valid correlation matrix
  return(
    list(
      corrmat = mat, 
      valid = FALSE
    )
  )
}

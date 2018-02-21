#' Block Correlation Matrix 
#' 
#' Generates a block correlation matrix for the SR data set. It randomly selects 
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
#' @param blocksize_drugs The block size for the drugs
#' @param rho_drugs The correlation for the drug-drug pairs in the blocks
#' @param blocksize_events The block size for the events
#' @param rho_drugs The correlation for the event-event pairs in the blocks
#' @param n_correlated_pairs The number of drug-event pairs that will have correlation \code{rho}
#' @param rho The correlation for the drug-event pairs 
#' @param max_iter Maximum number of iterations (Default: 1000) 
#' @param verbose Prints the number of iterations (Default: \code{FALSE})
#' 
#' @return \item{corrmat}{The correlation matrix}
#'         \item{covmat}{The covariance matrix}
#'         \item{valid}{\code{TRUE} when \code{corrmat} is a correlation matrix, \code{FALSE} otherwise} 
#'           
#' @seealso \code{\link{validCorrelation}}
#' @export
generateBlockCorrelationMatrix <- function(prob_drugs,
                                           prob_events,
                                           blocksize_drugs,
                                           rho_drugs,
                                           blocksize_events,
                                           rho_events,
                                           n_correlated_pairs,
                                           rho,
                                           max_iter = 1000,
                                           verbose = FALSE) {
  
  n_drugs <- length(prob_drugs)
  n_events <- length(prob_events)
  n <- n_drugs + n_events
  
  margprob <- c(prob_drugs, prob_events)

  corrmat_drugs <- generateBlockCorrelationMatrixRcpp(prob_drugs,
                                                      blocksize_drugs,
                                                      rho_drugs)
  corrmat_events <- generateBlockCorrelationMatrixRcpp(prob_events,
                                                       blocksize_events,
                                                       rho_events)

  # keep on generating correlation matrices until one is valid or 
  # the maximum number of iterations is reached
  for (iter in 1:max_iter) {
    
    if (verbose) {
      cat("Attempt no.", iter, "to generate a valid correlation matrix\n")
    }
    
    pairs <- returnRandomDrugEventPairsRcpp (prob_drugs, prob_events, n_correlated_pairs, rho) ; 
    
    mat <- diag(rep(1, n))
    
    mat[1:n_drugs, 1:n_drugs] <- corrmat_drugs
    mat[n_drugs+1:n_events, n_drugs+1:n_events] <- corrmat_events
    
    for (i in 1:nrow(pairs)) { 
      mat[pairs[i,1] + 1, pairs[i,2] + 1 + n_drugs] <- rho
      mat[pairs[i,2] + 1 + n_drugs, pairs[i,1] + 1] <- rho
    }
    
    # determine the covariance matrix
    covmat <- corr2cov(mat, margprob)   
    
    # We use the Cholesky decomposition to test for positive definiteness. 
    # Runs much faster for larger matrices than computing all eigenvalues
    if (corpcor::is.positive.definite(covmat, method = "chol")) {
      return(
        list(
          corrmat = mat, 
          covmat = covmat, 
          valid = TRUE
        )
      )
    }
  }
  
  # did not find a valid correlation matrix
  return(
    list(
      corrmat = mat, 
      covmat = covmat, 
      valid = FALSE
    )
  )
}

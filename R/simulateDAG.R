#' Spontaneous Reporting Data
#'
#' \strong{Valid Reports} Not any binary sequence is a valid report. Each report 
#' should contain at least one drug and at least one event  (otherwise it would 
#' never been sent to the spontaneous reporitng sytem). 
#' While generating reports, we make sure that this is indeed the case. When one does not
#' want to check the validity and wants to allow any binary sequence, one can set 
#' \code{valid_reports} to \code{FALSE}. 
#' 
#' @param n_reports Number of reports (Default: 100)
#' @param n_drugs Number of drugs (Default: 10)
#' @param n_events Number of adverse drug events (Default: 10)
#' @param alpha_drugs Alpha parameter for the drug marginal probabilities (Default: 1.0)
#' @param beta_drugs Beta parameter for the drug marginal probabilities (Default: 5.0)
#' @param alpha_events Alpha parameter for the event marginal probabilities (Default: 1.0)
#' @param beta_events Beta parameter for the event marginal probabilities (Default: 5.0)
#' @param n_correlated_pairs Number of drug-event pairs that will exhibit a nonzero correlation (Default: 2)
#' @param valid_reports If \code{TRUE}, only valid reports (with at least one drug and at least one event) are accepted. (Default: \code{TRUE})
#' @param max_iter Maximum number of iterations (Default: 1000) 
#' @param seed The seed used by the RNG (Default: automatically set)
#' @param verbose Verbosity (Default: \code{TRUE})
#'
#' @return \item{sr}{A data frame with the simulated reports. The columns are named \code{drug1}, \code{drug2} ..., \code{event1}, \code{event2}, ...}
#'         \item{prob_drugs}{A list with marginal probabilities of the drugs}
#'         \item{prob_events}{A list with marginal probabilities of the events}
#'         \item{corrmat}{The correlation matrix}
#'
#' @seealso \code{\link{create2x2Tables}}, 
#'          \code{\link{validReport}}              
#' @export
simulateDAG <- function(n_reports = 100, 
                        n_drugs = 10, 
                        n_events = 10,
                        alpha_drugs = 1.0,
                        beta_drugs = 5.0,
                        alpha_events = 1.0,
                        beta_events = 5.0,
                        method = 'er', 
                        exp_degree = 3, 
                        theta_drugs = 1.5,
                        n_correlated_pairs = 2,
                        theta = 2, 
                        valid_reports = TRUE, 
                        max_iter = 1000, 
                        seed = NULL,
                        verbose = TRUE) { 
  
  # set the seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # initialize a list with the marginal probabilities for 
  # the drugs and the events
  prob_drugs <- rep(NA, n_drugs) 
  prob_events <- rep(NA, n_events)
  
  # create the labels used for the drug and the event nodes
  drug_labels <- sprintf("drug%d", 1:n_drugs)
  event_labels <- sprintf("event%d", 1:n_events)
  
  DAG <- randDAG(n_drugs, exp_degree, method = method, weighted = TRUE, wFUN = function(m){theta_drugs})  
  
  # change the node labels for the drugs 
  graph::nodes(DAG) <- drug_labels
  
  # add the event nodes to the DAG
  DAG <- graph::addNode(event_labels, DAG)
  
  # add random drug and edge connections 
  drug_event_combinations <- expand.grid(drug_labels, event_labels) 
  colnames(drug_event_combinations) <- c("drug_id", "event_id")
  drug_event_pairs <- sample_n(drug_event_combinations, n_correlated_pairs, replace = FALSE)
  
  # all connections between drugs and events get the weight theta
  DAG <- graph::addEdge(sprintf("%s", drug_event_pairs$drug_id), sprintf("%s", drug_event_pairs$event_id), DAG, rep(theta, n_correlated_pairs))
  
  return(
    list(
      DAG = DAG 
    )
  )
}
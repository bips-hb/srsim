#' Spontaneous Reporting Data based on a DAG
#'
#' \code{simulateSRData} simulates a spontaneous reporting (SR) data set. It allows for setting
#' the partial probabilities of the drugs and the events, \strong{and} for setting the 
#' correlations not only between the drugs and the events, but also between the drugs and the 
#' events themselves. 
#' \cr\cr
#' Each report to a SR system contains two lists: 
#' \enumerate{
#'   \item the drugs to which the patient was (thought to be) exposed to, and 
#'   \item the adverse events that the patient experienced. 
#' }
#' We will represent each report as a binary vector. The first items represent whether 
#' the patient was exposed to the drug (\code{1} if he/she was, and \code{0} otherwise). 
#' The second part represents whether the patient experienced the event or not 
#' (\code{1} if he/she did, and \code{0} otherwise). For example, if there are 3 drugs and 
#' 4 events in total, a typical report could be 
#' \deqn{0 1 0 1 1 0 0} 
#' which represents that the patient was exposed to drug 2 (but not to drug 1 and 3), and 
#' experienced event 1 and 2 (but not 3 and 4). The simulation results in a binary matrix
#' where each row is a report. 
#' \cr\cr
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
#' @param prob_drugs Specific array with the marginal probability of the drugs (Default: Beta distribution)
#' @param alpha_events Alpha parameter for the event marginal probabilities (Default: 1.0)
#' @param beta_events Beta parameter for the event marginal probabilities (Default: 5.0)
#' @param prob_events Specific array with the marginal probability of the events (Default: Beta distribution)
#' @param method The method used for generating the random DAG for the drugs (Default: \code{'er'}, Erdos-Renyi graph)
#' @param exp_degree Average degree between the drugs in the random DAG (Default: 3)
#' @param theta_drugs The increase in odds-ratio when there is a connection between two drugs (Default: 1.5)
#' @param n_correlated_pairs Number of drug-event pairs that will exhibit a nonzero correlation (Default: 2)
#' @param theta Increase in odds-ratio when there is an edge going from a drug to an event (Default: 2.0)
#' @param valid_reports If \code{TRUE}, only valid reports (with at least one drug and at least one event) are accepted. (Default: \code{TRUE})
#' @param max_iter Maximum number of iterations (Default: 1000) 
#' @param seed The seed used by the RNG (Default: automatically set)
#' @param DAG The DAG used for the simulation. By default generated with the other parameters given to this function
#' @param verbose Verbosity (Default: \code{TRUE})
#'
#' @return \item{sr}{A data frame with the simulated reports. The columns are named \code{drug1}, \code{drug2} ..., \code{event1}, \code{event2}, ...}
#'         \item{prob_drugs}{A list with marginal probabilities of the drugs}
#'         \item{prob_events}{A list with marginal probabilities of the events}
#'         \item{adjencency_matrix}{The adjecency matrix of the DAG used}
#'
#' @seealso \code{\link{create2x2Tables}}, 
#'          \code{\link{generateCorrelationMatrix}}, 
#'          \code{\link{validReport}},
#'          \code{\link{validCorrelation}}               
#' @export
simulateSRDAG <- function(n_reports = 100, 
                          n_drugs = 10, 
                          n_events = 10,
                          alpha_drugs = 1.0,
                          beta_drugs = 5.0,
                          prob_drugs = rbeta(n_drugs, alpha_drugs, beta_drugs), 
                          alpha_events = 1.0,
                          beta_events = 5.0,
                          prob_events = rbeta(n_events, alpha_events, beta_events), 
                          method = 'er', 
                          exp_degree = 3, 
                          theta_drugs = 1.5,
                          n_correlated_pairs = 2,
                          theta = 2,
                          seed = NULL,
                          valid_reports = TRUE, 
                          verbose = TRUE) { 
  
  n_drugs <- length(prob_drugs)
  n_events <- length(prob_events)
  n <- n_drugs + n_events
  
  if (verbose) {
    cat("Creating DAG...") 
  }
  
  DAG <- SRSim::simulateDAG(n_drugs = length(prob_drugs), 
                           n_events = length(prob_events), 
                           prob_drugs = prob_drugs, 
                           prob_events = prob_events, 
                           method = method, 
                           exp_degree = exp_degree, 
                           theta_drugs = theta_drugs,
                           n_correlated_pairs = n_correlated_pairs,
                           theta = theta,
                           seed = seed)
  
  if (verbose) {
    cat("DONE Creating DAG...") 
  }
  
  # the marginal probabilities
  margprob <- c(prob_drugs, prob_events)
  
  #sr <- dplyr::as_tibble(data.frame(matrix(NA, nrow = n_reports, ncol = n)))
  #colnames(sr) <- c(sprintf("drug%d", 1:n_drugs), sprintf("event%d", 1:n_events))
  sr <- matrix(NA, nrow = n_reports, ncol = n)
  
  # get the nodes and sort them to in degree
  nodes <- DAG$nodes %>% 
    mutate(
      id = 1:n()
    ) %>% 
    arrange(in_degree) 
  
  if (verbose) { 
    pb <- txtProgressBar(min = 0, max = n_reports, initial = 0, char = "=",
                 style = 3, file = "")
  }
  
  # information needed:
  # nodes id 
  # nodes label
  # no. parents
  # parents of each node 
  # betas of each node 
  
  # matrix with the betas
  if (verbose) {
    cat("Obtaining the betas...")
  }
  
  betas <- matrix(log(DAG$adjacency_matrix), n, n)
  betas[betas == -Inf] <- 0
  
  if (verbose) {
    cat("DONE Obtaining the betas...")
  }
  
  if (verbose) {
    cat("Creating the reports...")
  }
  
  n_reports_generated <- 0
  
  while (n_reports_generated <= n_reports) {
    report <-
      simulateReportDAG(
        n_drugs,
        n_events,
        nodes$id,
        nodes$in_degree,
        DAG$max_in_degree,
        nodes$beta0,
        betas,
        verbose
      )
    if (valid_reports) {
      if (SRSim::validReport(t(matrix(report == 1)), n_drugs, n_events)) {
        sr[n_reports_generated,] <- report
        #sr <- dplyr::add_row(sr, report)
        n_reports_generated <- n_reports_generated + 1
      }
    } else {
      sr[n_reports_generated,] <- report
      #sr <- dplyr::add_row(sr, report)
      n_reports_generated <- n_reports_generated + 1
    }
    
    if (verbose) {
      setTxtProgressBar(pb, n_reports_generated)
    }
    
  }
  
  if (verbose) {
    cat("DONE generating SR data")
  }
  
  return(
    list( 
      sr = sr, 
      prob_drugs = prob_drugs,
      prob_events = prob_events, 
      adjecency_matrix = DAG$adjacency_matrix
    )
  )
  
  return(report)
  
  n_reports_generated <- 1
  
  while (n_reports_generated <= n_reports) {
    
    report <- rep(NA, n)
    filled_in <-1:n # keeps track of the binary variables that have not been drawn yet
    
    # continue till report is filled
    while (anyNA(report)) {
      # go through all indices that are not filled yet
      for (i in filled_in) {
        id <- nodes[i,]$id
        logreg <- nodes[i, ]$beta0
        # go over the parents
        for (parent_id in 1:DAG$max_in_degree) {
          parent_label <- nodes[i, ][[sprintf("parent%d", parent_id)]]
          if (is.na(parent_label)) {
            break
          }
          beta <- nodes[i, ][[sprintf("beta%d", parent_id)]]
          parent_id <- nodes[nodes$label == parent_label,]$id
          logreg <- logreg + beta * report[parent_id]
        }
        if (!is.na(logreg)) {
          report[id] <- rbinom(1, 1, exp(logreg) / (1 + exp(logreg)))
          filled_in <- filled_in[filled_in != i]
        }
      }
    }
    
    if (valid_reports) {
      if (SRSim::validReport(t(matrix(report == 1)), n_drugs, n_events)) {
        sr[n_reports_generated, ] <- report 
        n_reports_generated <- n_reports_generated + 1  
      }
    } else {
      sr[n_reports_generated, ] <- report 
      n_reports_generated <- n_reports_generated + 1 
    }
    
    if (verbose) { 
      setTxtProgressBar(pb, n_reports_generated)
    }
  }
  
  return(
    list( 
      sr = sr, 
      prob_drugs = prob_drugs,
      prob_events = prob_events, 
      adjecency_matrix = DAG$adjacency_matrix
    )
  )
}

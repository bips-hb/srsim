#' Simulating a Spontaneous Reporting System
#'
#' \code{simulateSRS} simulates a spontaneous reporting (SR) data set. The relationships between
#' the drugs and the adverse events (AEs) are specified by a directed acyclic graph (DAG),
#' see \code{\link{generateDAG}}.
#' \cr\cr
#' Each report to a SRS contains two lists:
#' \enumerate{
#'   \item the drugs to which the patient was (thought to be) exposed to, and
#'   \item the AEs that the patient experienced.
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
#' @param beta_drugs Beta parameter for the drug marginal probabilities (Default: 20.0)
#' @param alpha_events Alpha parameter for the event marginal probabilities (Default: 1.0)
#' @param beta_events Beta parameter for the event marginal probabilities (Default: 20.0)
#' @param n_innocent_bystanders Number of innocent bystanders (Default: 5)
#' @param bystander_prob The conditional probability of the innocent bystander being one when
#'                       the drug that is actually causing the AE is equal to 1. This parameter
#'                       corresponds to \eqn{\gamma} in the paper (Default: .9)
#' @param n_correlated_pairs Number of drug-AE pairs that are associated (Default: 2)
#' @param theta Increase in odds-ratio when there is an edge going from a drug to an AE (Default: 2.0).
#'              In case theta is a vector of length two, the odds ratio is drawn from a truncated
#'              Normal distribution with mean \code{theta[1]} and variance \code{theta[2]}
#' @param valid_reports If \code{TRUE}, only valid reports (with at least one drug and at least one AE)
#'                      are accepted. (Default: \code{TRUE})
#' @param seed The seed used by the RNG (Default: automatically set)
#' @param verbose Verbosity (Default: \code{TRUE})
#'
#' @return \item{sr}{A binary data frame with the simulated reports. The columns are
#'                   named \code{drug1}, \code{drug2} ..., \code{event1}, \code{event2}, ...}
#'         \item{dag}{The directed acycled graph as an \code{igraph} object}
#'         \item{nodes}{A tibble with all the information on each node/variate:
#'             \itemize{
#'                \item{\code{label}}{ The label for each node/variate}
#'                \item{\code{in_degree}}{ The number of edges pointing to the node}
#'                \item{\code{id}}{ The ID of each node (simple integer)}
#'                \item{\code{parent_id}}{ The ID of the parent node - if any. Otherwise equal to \code{-1}}
#'                \item{\code{margprob}}{ The marginal probability of the node/variate}
#'                \item{\code{beta0}}{ The intercept in the logistic regression model for that node}
#'                \item{\code{beta1}}{ The regression coefficient in the logistic regression model for the parent}
#'             }
#'         }
#'         \item{prob_drugs}{A vector with marginal probabilities of the drugs}
#'         \item{prob_events}{A vector with marginal probabilities of the events}
#'
#'
#' @seealso \code{\link{convert2Tables}},
#'          \code{\link{generateDAG}}
#' @export
simulateSRS <- function(n_reports = 100,
                        n_drugs = 10,
                        n_events = 10,
                        alpha_drugs = 1.0,
                        beta_drugs = 20.0,
                        alpha_events = 1.0,
                        beta_events = 20.0,
                        n_innocent_bystanders = 5,
                        bystander_prob = 0.9,
                        n_correlated_pairs = 2,
                        theta = 2,
                        valid_reports = TRUE,
                        seed = NULL,
                        verbose = TRUE) {

  # set the seed ---
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Silence global variable NOTE (could use globalVariables() if need be)
  margprob <- NULL

  n        <- n_drugs + n_events
  sr       <- matrix(NA, nrow = n_reports, ncol = n) # matrix that will contain the reports

  beta <- log(theta) # coefficient for the logistic model given the desired OR, theta

  # generate the DAG ---

  if (verbose) { cat("Creating DAG...\n") }

  DAG <- generateDAG(
    n_drugs = n_drugs,
    n_events = n_events,
    n_innocent_bystanders = n_innocent_bystanders,
    n_correlated_pairs = n_correlated_pairs
  )

  if (verbose) { cat("DONE Creating DAG...\n") }

  drug_labels <- sprintf("drug%d", 1:n_drugs)
  event_labels <- sprintf("event%d", 1:n_events)

  # a data frame that contains all the data for each node (or variate)
  nodes <- tibble::tibble(
    label = c(drug_labels, event_labels),
    in_degree = igraph::degree(DAG, mode = "in"),
    id = 1:n,
    parent_id = -1,
    margprob = c(
      rbeta(n_drugs, alpha_drugs, beta_drugs),
      rbeta(n_events, alpha_events, beta_events)
    ),
    beta0 = log(margprob / (1 - margprob)),
    beta1 = 0
  )

  # create a tibble with all the edges
  edgelist <- igraph::as_edgelist(DAG)
  edges <- tibble::tibble(
    from = edgelist[,1],
    to = edgelist[,2]
  )

  # walk through the edges and add the information to the nodes tibble
  for (i in 1:nrow(edges)) {
    from <- edges[i,]$from
    to   <- edges[i,]$to

    # add the parent to the nodes tibble
    nodes[nodes$id == to,]$parent_id <- from

    # add the appropriate beta value for this edge
    if (to <= n_drugs) { # a drug
      no_bystander_prob <- rbeta(1, alpha_drugs, beta_drugs)
      nodes[nodes$id == to,]$beta0 <-
        log(no_bystander_prob / (1 - no_bystander_prob))
      nodes[nodes$id == to,]$beta1 <-
        log(bystander_prob / (1 - bystander_prob)) - log(no_bystander_prob / (1 - no_bystander_prob))

      # recompute the marginal probability
      margprob_parent <-
        nodes[nodes$id == from, ]$margprob
      nodes[nodes$id == to, ]$margprob <-
        margprob_parent*bystander_prob + (1 - margprob_parent)*no_bystander_prob
    } else { # an event
      # determine beta1
      if (length(theta) == 1) {
        nodes[nodes$id == to,]$beta1 <- log(theta)
      } else {
        nodes[nodes$id == to,]$beta1 <- log(max(1, rnorm(1, theta[1], theta[2])))
      }
    }
  }

  # determine the beta0's (the intercepts) for the events
  for (i in (n_drugs + 1):n) {
    node <- nodes[i,]
    if (node$in_degree != 0) {
      # optimize the beta0
      margprob_parent <-
        nodes[nodes$id == node$parent_id, ]$margprob
      res <-
        optim(
          node$beta0,
          f_beta0,
          method = "BFGS",
          margprob = node$margprob,
          beta1 = node$beta1,
          margprob_parent = margprob_parent
        )
      nodes[i, ]$beta0 <- res$par
    }
  }

  nodes$in_degree <- as.integer(nodes$in_degree)
  nodes$parent_id <- as.integer(nodes$parent_id)

  # generating the actual reports
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = n_reports, initial = 0, char = "=",
                         style = 3, file = "")
  }

  n_reports_generated <- 0
  while (n_reports_generated <= n_reports) {
    report <-
      simulateReport(
        n_drugs,
        n_events,
        nodes$in_degree,
        nodes$beta0,
        nodes$beta1,
        nodes$parent_id
      )

    if (valid_reports) {
      # check whether the report is valid
      if (validReport(t(matrix(report == 1)), n_drugs, n_events)) {
        sr[n_reports_generated,] <- report
        n_reports_generated <- n_reports_generated + 1
      }
    } else {
      sr[n_reports_generated,] <- report
      n_reports_generated <- n_reports_generated + 1
    }

    if (verbose) { setTxtProgressBar(pb, n_reports_generated) }
  }

  colnames(sr) <- c(sprintf("drug%d", 1:n_drugs), sprintf("event%d", 1:n_events))
  sr <- tibble::as_tibble(sr)

  return(
    list(
      sr = sr,
      dag = DAG,
      nodes = nodes,
      prob_drugs = nodes$margprob[1:n_drugs],
      prob_events = nodes$margprob[(n_drugs+1):n]
    )
  )
}

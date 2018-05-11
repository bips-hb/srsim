' Random Directed Acyclic Graph for SR Data
#'
#' \code{simulateDAG} creates a random directed acyclic graph that can be used
#' for simulating SR data
#' 
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
#' @param seed The seed used by the RNG (Default: automatically set)
#' @param verbose Verbosity (Default: \code{TRUE})
#'
#' @return \item{nodes}{A tibble containing all the information on each node/variate}
#'         \item{DAG}{The DAG as an \code{igraph} object} 
#'         \item{prob_drugs}{A list with marginal probabilities of the drugs}
#'         \item{prob_events}{A list with marginal probabilities of the events}
#'         \item{max_in_degree}{The maximal in-degree found in the graph \code{DAG}}
#'         \item{adjacency_matrix}{The adjencency matrix for the graph \code{DAG}} 
#'  
#' @seealso \code{\link{create2x2Tables}}, 
#'          \code{\link{validReport}}              
#' @export
createDAG <- function(n_drugs = 10, 
                        n_events = 10,
                        n_innocent_bystanders = 5, 
                        n_correlated_pairs = 2,
                        seed = NULL,
                        verbose = TRUE) { 
  
  if (n_drugs %% 2 == 1) { 
    stop("Number of drugs should be even") 
  }
  
  if (n_innocent_bystanders > n_drugs / 2) {
    stop("Number of innocent bystanders cannot be larger than half of the drugs") 
  }
  
  # set the seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n <- n_drugs + n_events # total number of variates in the graph 
  
  DAG <- igraph::graph.empty(n = n, directed = TRUE)
  
  # add the innocent bystanders 
  if (n_innocent_bystanders == 0) { 
    edges <- c() 
  } else { 
    edges <- c(matrix(c(1:n_innocent_bystanders, (n_drugs/2 + 1):(n_drugs/2+n_innocent_bystanders)), 2, byrow = TRUE)) 
  }
  
  
  DAG <- igraph::graph(edges, n = n, directed = TRUE)
  
  # add the random connections between the drugs and the events
  drug_event_pairs <- expand.grid(1:(n_drugs/2), (n_drugs+1):n) 
  drug_event_pairs <- sample_n(drug_event_pairs, n_correlated_pairs, replace = FALSE)
  edges <- c(matrix(c(drug_event_pairs$Var1, drug_event_pairs$Var2), 2, byrow = TRUE))
  
  DAG <- igraph::add.edges(DAG, edges) 
  
  # change the labels of the vertices 
  # drug_labels <- sprintf("drug%d", 1:n_drugs)
  # event_labels <- sprintf("event%d", 1:n_events)
  # DAG <- igraph::set.vertex.attribute(DAG, "name", value = c(drug_labels, event_labels))
  # 
  return(DAG)
  
  g<-graph(edges, n=max(edges), directed=TRUE)
  
  
  
  # add random drug and edge connections 
  drug_event_pairs <- expand.grid(drug_labels, event_labels) 
  colnames(drug_event_pairs) <- c("drug_id", "event_id")
  drug_event_pairs <- sample_n(drug_event_pairs, n_correlated_pairs, replace = FALSE)
  
  # all connections between drugs and events get the weight theta
  graph::addEdge(sprintf("%s", drug_event_pairs$drug_id), sprintf("%s", drug_event_pairs$event_id), DAG, rep(theta, n_correlated_pairs))
  
  
  
  
  
  
  
  DAG <- pcalg::randDAG(n_drugs, exp_degree, method = method, weighted = TRUE, wFUN = function(m){theta_drugs})  
  graph::nodes(DAG) <- drug_labels            # change the node labels for the drugs 
  DAG <- graph::addNode(event_labels, DAG)    # add the event nodes to the DAG
  
  # add random drug and edge connections 
  drug_event_pairs <- expand.grid(drug_labels, event_labels) 
  colnames(drug_event_pairs) <- c("drug_id", "event_id")
  drug_event_pairs <- sample_n(drug_event_pairs, n_correlated_pairs, replace = FALSE)
  
  # all connections between drugs and events get the weight theta
  DAG <- graph::addEdge(sprintf("%s", drug_event_pairs$drug_id), sprintf("%s", drug_event_pairs$event_id), DAG, rep(theta, n_correlated_pairs))
  
  adjacency_matrix <- as(DAG, "matrix")
  
  DAG <- igraph::igraph.from.graphNEL(DAG)    # transform to igraph
  
  max_in_degree <- max(igraph::degree(DAG, mode = "in")) # maximal in-degree in the random DAG
}
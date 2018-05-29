#' Directed Acyclic Graph for SR Data
#'
#' \code{createDAG} creates a directed acyclic graph that can be used
#' for simulating SR data
#' 
#' @param n_drugs Number of drugs (Default: 10)
#' @param n_events Number of adverse drug events (Default: 10)
#' @param n_innocent_bystanders Number of innocent bystanders (Default: 5)
#' @param n_correlated_pairs Number of drug-event pairs that will be associated (Default: 2)
#'
#' @return The DAG as an \code{igraph} object
#' @export
createDAG <- function(n_drugs = 10,
                      n_events = 10,
                      n_innocent_bystanders = 5,
                      n_correlated_pairs = 2) {
  
  if (n_drugs %% 2 == 1) { 
    stop("Number of drugs should be even") 
  }
  
  if (n_correlated_pairs > n_events) { 
    stop("Number of events should be larger than the number of correlated pairs") 
  }
  
  if (n_innocent_bystanders > n_drugs / 2) {
    stop("Number of innocent bystanders cannot be larger than half of the drugs") 
  }
  
  n <- n_drugs + n_events # total number of variates in the graph 
  
  DAG <- igraph::graph.empty(n = n, directed = TRUE)
  
  # add the innocent bystanders 
  if (n_innocent_bystanders != 0) { 
    edges <- c(matrix(c(1:n_innocent_bystanders, (n_drugs/2 + 1):(n_drugs/2+n_innocent_bystanders)), 2, byrow = TRUE)) 
    DAG <- igraph::graph(edges, n = n, directed = TRUE)
  }
  
  # add the random connections between the drugs and the events
  edges <- c(matrix(c(1:n_correlated_pairs, (n_drugs+1):(n_drugs+n_correlated_pairs)), 2, byrow = TRUE))
  DAG <- igraph::add.edges(DAG, edges) 
  
  return(DAG)
}
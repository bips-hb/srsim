#' Random Directed Acyclic Graph for SR Data
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
simulateDAG <- function(n_drugs = 10, 
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
                        verbose = TRUE) { 
  
  # set the seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n <- n_drugs + n_events # total number of variates in the graph 
  
  # create the labels used for the drug and the event nodes
  drug_labels <- sprintf("drug%d", 1:n_drugs)
  event_labels <- sprintf("event%d", 1:n_events)
  
  beta_drugs <- log(theta_drugs) # the regression coefficient for the drug-drug pairs
  beta <- log(theta) # the regression coefficient for the drug-event pairs 
  
  ### create a random DAG for the drugs --------------
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
  ### DONE with creating the random DAG --------------
  
  # a data frame that contains all the data for each node (or variate)
  nodes <- dplyr::tibble(
     label = c(drug_labels, event_labels), 
     in_degree = igraph::degree(DAG, mode = "in"), # as.factor? 
     n_done = 0, # number of incoming edges processed
     margprob = c(prob_drugs, prob_events)
  ) 
  
  # add columns for the betas 
  for (beta_label in sprintf("beta%d", 0:max_in_degree)) {
    nodes <- tibble::add_column(nodes, !!(beta_label) := 0)
  }
  
  # add columns for the parents 
  for (parent_label in sprintf("parent%d", 1:max_in_degree)) { 
    nodes <- tibble::add_column(nodes, !!(parent_label) := NA)
  }
  
  # create also a tibble with all the edges
  edgelist <- igraph::as_edgelist(DAG)
  edges <- tibble(
    from = edgelist[,1], 
    to = edgelist[,2]
  )
  
  # walk through the edges and add the information to the nodes tibble
  for (i in 1:nrow(edges)) {
    to   <- edges[i,]$to 
    from <- edges[i,]$from 
    
    current_parent <- nodes[nodes$label == to,]$n_done + 1 # the current parent that needs to be filled in
    
    beta_label   <- sprintf("beta%d", current_parent)
    parent_label <- sprintf("parent%d", current_parent)

    # add the parent to the nodes tibble
    nodes[nodes$label == to,][[parent_label]] <- from
    
    # add the appropriate beta value for this edge
    if (startsWith(to, "drug")) { 
      nodes[nodes$label == to,][[beta_label]] <- beta_drugs
    } else { 
      nodes[nodes$label == to,][[beta_label]] <- beta
    }
    
    # one more done
    nodes[nodes$label == to,]$n_done <- nodes[nodes$label == to,]$n_done + 1
  }
  
  ### Set the intercepts, beta0, appropriately 
  nodes <- nodes %>% dplyr::mutate(
    beta0 = log(margprob / (1 - margprob))
  ) %>% dplyr::select(
    -n_done 
  )
  
  # walk through the nodes. If there are parents, update the beta0
  for (i in 1:nrow(nodes)) {
    n_parents <- nodes[i,]$in_degree
    # if there are parents, iterate over them
    if (n_parents != 0) {
      for (current_parent in 1:max_in_degree) { 
        beta_label   <- sprintf("beta%d", current_parent)
        parent_label <- sprintf("parent%d", current_parent)
        
        parent <- nodes[i,][[parent_label]]
        if (is.na(parent)) {
          break  
        }
        
        # get the beta and the marginal probability of the parent
        beta_parent <- nodes[nodes$label == parent,][[beta_label]] 
        margprob_parent <- nodes[nodes$label == parent,]$margprob
        
        # simply subtract their product to get the appropriate beta0
        nodes[i,]$beta0 <- nodes[i,]$beta0 - beta_parent*margprob_parent 
      }
    }
  }
  
  return(
    list(
      nodes = nodes,
      DAG = DAG, 
      prob_drugs = prob_drugs,
      prob_events = prob_events,
      max_in_degree = max_in_degree,
      adjacency_matrix = adjacency_matrix
    )
  )
}
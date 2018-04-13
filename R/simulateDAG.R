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
                        prob_drugs = rbeta(n_drugs, alpha_drugs, beta_drugs), 
                        alpha_events = 1.0,
                        beta_events = 5.0,
                        prob_events = rbeta(n_events, alpha_events, beta_events), 
                        method = 'er', 
                        exp_degree = 3, 
                        theta_drugs = 1.5,
                        n_correlated_pairs = 2,
                        theta = 2, 
                        valid_reports = TRUE, 
                        max_iter = 1000, 
                        seed = NULL,
                        verbose = TRUE) { 
  
  n <- n_drugs + n_events 
  beta_drugs <- log(theta_drugs) # the regression coefficient for the drug-drug pairs
  beta <- log(theta) # the regression coefficient for the drug-event pairs 
  
  # set the seed
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # create the labels used for the drug and the event nodes
  drug_labels <- sprintf("drug%d", 1:n_drugs)
  event_labels <- sprintf("event%d", 1:n_events)
  
  # create a random DAG for the drugs 
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
  
  # get the adjacency matrix with the weights
  adjacency_matrix <- as(DAG, "matrix")
  
  # transform to igraph
  DAG <- igraph::igraph.from.graphNEL(DAG)
  
  # add the attributes to the nodes of the graph
  DAG <- DAG %>% 
    igraph::set_vertex_attr("margprob", index = 1:n, value = c(prob_drugs, prob_events)) %>% 
    igraph::set_vertex_attr("done", index = 1:n, value = FALSE) 
  
  
  max_in_degree <- max(igraph::degree(DAG, mode = "in"))
  
  nodes <- dplyr::tibble(
     label = c(drug_labels, event_labels), 
     in_degree = as.factor(igraph::degree(DAG, mode = "in")), 
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
  
  edgelist <- igraph::as_edgelist(DAG)
  
  edges <- tibble(
    from = edgelist[,1], 
    to = edgelist[,2]
  )
  
  # walk through the edges and add the information to the nodes tibble
  for (i in 1:nrow(edges)) {
    to <- edges[i,]$to 
    from <- edges[i,]$from 
    
    current_parent <- nodes[nodes$label == to,]$n_done + 1
    beta_label <- sprintf("beta%d", current_parent)
    parent_label <- sprintf("parent%d", current_parent)

    if (startsWith(to, "drug")) { 
      nodes[nodes$label == to,][[beta_label]] <- beta_drugs
    } else { 
      nodes[nodes$label == to,][[beta_label]] <- beta
    }
    
    nodes[nodes$label == to,][[parent_label]] <- from
    
    nodes[nodes$label == to,]$n_done <- nodes[nodes$label == to,]$n_done + 1
  }
  
  # set the beta0 for the logistic regression
  nodes <- nodes %>% mutate(
    beta0 = log(margprob / (1 - margprob)),
    n_done = 0
  )
  
  print(nodes)
  
  for (i in 1:nrow(nodes)) {
    # update the beta0's in case there are parents
    n_parents <-nodes[i,]$in_degree
    print(n_parents)
    if (n_parents != 0) {
      print(1:(as.numeric(n_parents)-1))
      for (current_parent in 1:max_in_degree) { 
        beta_label <- sprintf("beta%d", current_parent)
        parent_label <- sprintf("parent%d", current_parent)
        
        parent <- nodes[i,][[parent_label]]
        
        if (is.na(parent)) {
         break  
        }
        
        beta_parent <- nodes[nodes$label == parent,][[beta_label]] 
        margprob_parent <- nodes[nodes$label == parent,]$margprob
        
        print(parent)
        print(beta_parent)
        print(margprob_parent)
        
        nodes[i,]$beta0 <- nodes[i,]$beta0 - beta_parent*margprob_parent 
      }
    }
  }
  
  ###### generating the reports 
  n_reports_done <- 0 # number of valid reports generated

  # TODO check whether report is valid 
  
  
  
  
  
  return(
    list(
      DAG = DAG, 
      adjacency_matrix = adjacency_matrix,
      nodes = nodes
    )
  )
}
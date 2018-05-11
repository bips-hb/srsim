library(SRSim)
library(ggplot2)

seed <- 2

res <-  simulateDAG(
  n_drugs = 3,
  n_events = 3,
  method = 'er',
  exp_degree = 0,
  theta_drugs = 1.5,
  n_correlated_pairs = 2,
  theta = 2,
  seed = NULL,
  verbose = TRUE
)

plot(res$DAG)

########################

n_drugs = 10
n_events = 10
alpha_drugs = 1.0
beta_drugs = 5.0
prob_drugs = rbeta(n_drugs, alpha_drugs, beta_drugs)
alpha_events = 1.0
beta_events = 5.0
prob_events = rbeta(n_events, alpha_events, beta_events)

n_reports <- 100000

res <- simulateSRDAG(
  n_reports = n_reports,
  prob_drugs = prob_drugs,
  prob_events = prob_events,
  method = 'er',
  exp_degree = 0,
  theta_drugs = 0,
  n_correlated_pairs = 0,
  theta = 10,
  seed = NULL,
  valid_reports = FALSE,
  verbose = TRUE
)

sr <- res$sr

tables <- SRSim::create2x2TablesDAG(res) %>% 
  mutate(est_or = (a * d) / (b * c), 
         difference_or = or - est_or)

ggplot2::ggplot(data = tables) + 
  ggplot2::geom_histogram(mapping = ggplot2::aes(x = difference_or))

hist(c(prob_drugs, prob_events) - colSums(sr) / n_reports, breaks=20)

sort(c(prob_drugs, prob_events) - colSums(sr) / n_reports)

# conclusion: the simulation is off! 



### Step-by step...

theta <- 3
theta_drugs <- 1.5
n_correlated_pairs <- 25
exp_degree <- 3

DAG <- SRSim::simulateDAG(
  n_drugs = n_drugs,
  n_events = n_events,
  prob_drugs = prob_drugs,
  prob_events = prob_events,
  method = 'er',
  exp_degree = exp_degree,
  theta_drugs = theta_drugs,
  n_correlated_pairs = n_correlated_pairs,
  theta = theta,
  seed = NULL
)

plot(DAG$DAG)
nodes <- DAG$nodes

igraph::as_edgelist(DAG$DAG)


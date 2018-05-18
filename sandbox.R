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





rep = 100

or = c()

for (i in 1:rep) { 
  res <-
    SRSim::simulateSR(
      n_reports = 10000,
      valid_reports = F,
      theta_drugs = 1000,
      theta = 5
    )
  
  a = sum(res$sr$drug1 * res$sr$drug6)
  b = sum((1 - res$sr$drug1) * res$sr$drug6)
  c = sum(res$sr$drug1 * (1 - res$sr$drug6))
  d = sum((1 - res$sr$drug1) * (1 - res$sr$drug6))

  or <- c(or, (a * d) / (b * c))

}


a = sum(res$sr$drug1 * res$sr$event9)
b = sum((1 - res$sr$drug1) * res$sr$event9)
c = sum(res$sr$drug1 * (1 - res$sr$event9))
d = sum((1 - res$sr$drug1) * (1 - res$sr$event9))


a = sum(res$sr$drug6 * res$sr$event9)
b = sum((1 - res$sr$drug6) * res$sr$event9)
c = sum(res$sr$drug6 * (1 - res$sr$event9))
d = sum((1 - res$sr$drug6) * (1 - res$sr$event9))

a = sum(res$sr$drug7 * res$sr$event5)
b = sum((1 - res$sr$drug7) * res$sr$event5)
c = sum(res$sr$drug7 * (1 - res$sr$event5))
d = sum((1 - res$sr$drug7) * (1 - res$sr$event5))





rep = 20

or_drug_events = c()
or_drug_bystanders = c()
or_bystanders_events = c()

n_innocent_bystanders <- 5
theta <- c(1.25,1)
theta_drugs <- 1000

drugs <- sprintf("drug%d", 1:n_innocent_bystanders)
bystanders <- sprintf("drug%d", (n_innocent_bystanders+1):10)



for (i in 1:rep) { 
  
  tryCatch({
    res <-
      SRSim::simulateSR(
        n_reports = 10000,
        beta_drugs = 20, 
        beta_events = 20,
        n_innocent_bystanders = n_innocent_bystanders,
        valid_reports = F,
        theta_drugs = theta_drugs,
        theta = theta
      )
    
    e <- res$nodes %>% filter(startsWith(label,"event"), parent_id != -1)
    events <- e$label 
    drug_events <- sprintf("drug%d", e$parent_id)
    
    # or drug events  
    for (i in 1:nrow(e)) {
      event <- e[i,]
      drug <- sprintf("drug%d", event$parent_id)
      
      bystander <- sprintf("drug%d", event$parent_id + 5)
      
      a = sum(res$sr[[drug]] * res$sr[[event$label]])
      b = sum((1 - res$sr[[drug]]) * res$sr[[event$label]])
      c = sum(res$sr[[drug]] * (1 - res$sr[[event$label]]))
      d = sum((1 - res$sr[[drug]]) * (1 - res$sr[[event$label]]))
      
      
      or_drug_events <- c(or_drug_events, (a * d) / (b * c))
      
      a = sum(res$sr[[bystander]] * res$sr[[event$label]])
      b = sum((1 - res$sr[[bystander]]) * res$sr[[event$label]])
      c = sum(res$sr[[bystander]] * (1 - res$sr[[event$label]]))
      d = sum((1 - res$sr[[bystander]]) * (1 - res$sr[[event$label]]))
      
      or_bystanders_events <- c(or_bystanders_events, (a * d) / (b * c)) 
    }
    
    # or drug bystanders
    for (k in 1:5) {
      l = k + 5
      i = sprintf("drug%d", k)
      j = sprintf("drug%d", l)
      a = sum(res$sr[[i]] * res$sr[[j]])
      b = sum((1 - res$sr[[i]]) * res$sr[[j]])
      c = sum(res$sr[[i]] * (1 - res$sr[[j]]))
      d = sum((1 - res$sr[[i]]) * (1 - res$sr[[j]]))
      
      or_drug_bystanders <- c(or_drug_bystanders, (a * d) / (b * c))
    }
  }, error=function(e){})
  
 
  
}

hist(or_drug_bystanders, breaks = 30)
hist(or_drug_events, breaks = 30)
hist(or_bystanders_events, breaks = 30)

hist(or_drug_events, breaks = 30, col=rgb(1,0,0,0.5))
hist(or_bystanders_events, breaks = 30, col=rgb(0,0,1,0.5), add = T)







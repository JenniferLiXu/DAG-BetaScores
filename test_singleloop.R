#Creating a small DAG model in R to 
#test the proposed method for obtaining beta values

install.packages("bnlearn")
install.packages("BiDAG")
install.packages("pcalg")
install.packages("igraph") 

library(BiDAG) 
library(pcalg)
library(igraph)
library(bnlearn)

# load the necessary functions
#source('./edgerevandstructure/newedgerevfns.R')
#source('./edgerevandstructure/newedgerevmove.R')
source('./orderMCMC_betas.R')
source('./orderscore_betas.R')
source('./samplefns-beta.R')
source('./scoring/scorefns.R')
source('./importance_sampling_beta.R')
source('./calculateBetaScoresArray.R')

# Example: Generating a random dataset
set.seed(123)

n <- 7
scoreParam <- BiDAG::scoreparameters("bge", BiDAG::Boston[1:100,1:n])
#itfit<-learnBN(scoreParam,algorithm="orderIter")
#maxEC<-getDAG(itfit,cp=FALSE)

#posterior probabilities of edges by averaging over a sample of DAGs obtained via an MCMC scheme.
samplefit<-sampleBN(scoreParam, "order")
edgesposterior<-edgep(samplefit, pdag=TRUE, burnin=0.2)

# DAG_Test <- list()
# DAG_Test[[1]] <- maxEC
# target_beta <- calculateBetaScoresArray(DAG_Test, k = length(DAG_Test), n)

# Move probabilities
prob1<-99
if(n>3){ prob1<-round(6*99*n/(n^2+10*n-24)) }
prob1<-prob1/100
moveprobs<-c(prob1,0.99-prob1,0.01)
moveprobs<-moveprobs/sum(moveprobs) # normalisation
if(!(length(moveprobs)==3)){print('Vector of move probabilities has the wrong length!')}

# Start with empty graph (OR user defined graph OR from other algo. and learn the beta matrix)
starting_dag <- list(matrix(c(0), nrow = n, ncol = n, byrow = TRUE))
betas_init <- calculateBetaScoresArray(starting_dag, k = 1 ,n)[,,1]
#base_score <- BiDAG::DAGscore(scoreParam, starting_dag[[1]])
base_score <- 0

# Initialization
iter <- 400
weighted_betas <- list(betas_init) # Starting beta_matrix
order <- list(c(1:n)) # Starting order

ess_DAGs <- numeric()
DAG_scores <- list()
DAG <- list() # Store the weighted DAG from each iteration based on the update weighted beta_matrix
DAG_total <- matrix(0, n, n)
#totallogscore<-list(sum(orderscore_betas(n,c(1:n), weighted_betas[[1]], order[[1]]))) #starting score

#Initialize for OrderMCMC
order_iter <-  100
order_stepsize <- 10

differences <- numeric()  # Initialize a vector to store the last 'num_steps' differences
diff_BiDAGs <- numeric() 

# Looping
for (i in 1:iter) {
  beta_prev <- weighted_betas[[i]]
  order_prev <- order[[i]]
  
  # Sampling the orders
  example <- orderMCMC_betas(n,startorder = order_prev ,iterations = order_iter, 
                             betas = beta_prev,
                             stepsave = order_stepsize, moveprobs)
  
  #Store the last order from the chain
  permy <- unlist(example[[4]][length(example[[4]])])
  order[[i+1]] <- permy

  # Initialize a list to store the sampled DAGs using the last order form OrderMCMC and their info
  sampled_DAGs <- vector("list", 10)
  incidence_matrices <- list()
  sampled_DAGs <- lapply(1:10, function(x) samplescore(n, beta_prev, permy, base_score))
  incidence_matrices <-lapply(sampled_DAGs, function(dag) dag$incidence) #store adjacency matrix of a sampled DAG each 'stepsave'
  incidence_logscore<-lapply(sampled_DAGs, function(dag) dag$logscore) #and log score of a sampled DAG

  ### Update beta matrix using importance sampling from the sampled DAGs
  is_results <- importance_DAG_prev(DAGs = incidence_matrices, s_betas = incidence_logscore)
  
  # Vectorized multiplication of each matrix by its corresponding weight and then summing up the matrices
  weights <- is_results$importance_weights
  #cat("weights",i, ":" , weights, "\n")
  beta_values <- calculateBetaScoresArray(incidence_matrices, k = length(incidence_matrices) ,n) #a list of matrices
  
  # Update the beta_metrix using the weights
  weighted_betas[[i + 1]]  <- Reduce("+", lapply(1:length(weights), 
                                                 function(k) beta_values[,,k] * weights[k]))
  ess_DAGs[i] <-is_results$ess_value
  
  ### Stopping criteria for beta matrix convergence
  # Calculate the Frobenius norm of the difference between current and previous matrices
  current_beta_matrix <- weighted_betas[[i + 1]]
  previous_beta_matrix <- weighted_betas[[i]]
  
  difference <- norm(current_beta_matrix - previous_beta_matrix, type = "F")
  # Update the differences vector
  differences <- c(differences, difference)
  
  # The weight based on the current/proposed(start/end) orderlogscore using the current/proposed beta_matrix
  order_s_beta_s <- sum(orderscore_betas(n,c(1:n), weighted_betas[[i]], order[[i]]))
  order_s_beta_e <- sum(orderscore_betas(n,c(1:n), weighted_betas[[i+1]], order[[i]]))
  order_e_beta_s <- sum(orderscore_betas(n,c(1:n), weighted_betas[[i]], order[[i+1]]))
  order_e_beta_e <- sum(orderscore_betas(n,c(1:n), weighted_betas[[i+1]], order[[i+1]]))
  #cat("order_beta",i, ":" , order_s_beta_s ,order_s_beta_e ,order_e_beta_s,order_e_beta_e,"\n")
  
  ratio <- (order_s_beta_s + order_e_beta_s) / (order_s_beta_e + order_e_beta_e)
  #cat("ratio",i, ":" , ratio,"\n")
  DAG[[i]] <- ratio*samplescore(n, weighted_betas[[i + 1]], order[[i+1]], base_score)$incidence
  
  DAG_total <- DAG_total + DAG[[i]]
  diff_BiDAG <- norm(as.matrix(DAG_total/i - edgesposterior), type = "F")
  # if(i < 7){
  #   DAG_total <- DAG_total + DAG[[i]]
  #   diff_BiDAG <- norm(as.matrix(DAG_total/i - edgesposterior), type = "F")
  # }else if(i == 7){
  #   DAG_total <- DAG[[i]]
  #   diff_BiDAG <- norm(as.matrix(DAG_total - edgesposterior), type = "F")
  # }else{
  #   DAG_total <- DAG_total + DAG[[i]]
  #   diff_BiDAG <- norm(as.matrix(DAG_total/(i-6) - edgesposterior), type = "F")
  # }

  
  diff_BiDAGs <- c(diff_BiDAGs, diff_BiDAG)

}

print(ess_DAGs)
plot(diff_BiDAGs)
plot(diff_BiDAGs[-c(1:100)])

# Plotting the differences for the beta matrices
plot(#differences[-c(1:3)], 
  differences, 
  col = "blue",
  type = "b", main = "Differences for Beta_Matrix Per Iteration", 
  xlab = "Iteration", ylab = "Difference")












########### Graphing part ############

### Build a Consensus Graph from sampled DAGs
# Include an edge in the consensus graph if it appears in a significant number of DAGs
# Initialize matrices to store edge frequencies and cumulative beta values


edge_freq_perc <- DAG_total/length(DAG)

consensus_graph <- edge_freq_perc >= 0.22
consensus_adj_mat <- ifelse(consensus_graph, 1, 0)

# Create an igraph graph from the adjacency matrix
graph <- graph_from_adjacency_matrix(consensus_adj_mat, mode = "undirected", diag = FALSE)

# Plot the graph
plot(graph, 
     main = "Consensus Graph",
     edge.arrow.size = 0.5,
     vertex.color = "lightblue",
     vertex.size = 15,
     vertex.label.color = "black",
     vertex.label.cex = 0.8)


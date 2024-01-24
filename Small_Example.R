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

seedset<-1 # set a seed?
seednumber<-101 # an example used


# load the necessary functions
#source('./edgerevandstructure/newedgerevfns.R')
#source('./edgerevandstructure/newedgerevmove.R')
source('./orderMCMC_betas.R')
source('./orderscore_betas.R')
source('./samplefns-beta.R')
source('./scoring/scorefns.R')

#Use the function to calculate beta values
source('./calculateBetaScoresArray.R')


# Example: Generating a random dataset
set.seed(123)
data <- matrix(rnorm(1400), ncol = 14)  # A dataset with 14 variables

# Specify the conditional independence test
# gaussCItest is used for continuous data (assuming Gaussian distribution)
ciTest <- gaussCItest

# Apply the PC algorithm
pcFit <- pc(suffStat = list(C = cor(data), n = nrow(data)), 
            indepTest = ciTest, p = 14, alpha = 0.05)  
# Convert the output to an adjacency matrix
adjMatrix <- as(pcFit, "amat")

#DAG <-samplescore(n,betas,permy)
#calculateBetaScoresArray(sampledDAGs = DAG, n)

#Score of a DAG with given beta matrix
DAGscore_under_betas <- function(incidence, betas){
  diag(betas) <- 0
  sampledscore<-0
  for (child in 1:n){
      for (parent in 1:n){
        sampledscore <- sampledscore + betas[parent, child]*incidence[parent, child]   
      }
  }
  return(sampledscore)
}

#Importance Weights for the sampled DAGs
importance_DAG <- function(DAGs, betas){
  differ_score <- numeric()
  
  #betas = beta_values
  for (i in 1:length(DAGs)){
    target_score <- BiDAG::DAGscore(scoreParam, DAGs[[i]])
    #cat("target_score",i, ":" , target_score, "\n")
    beta_score <- DAGscore_under_betas(DAGs[[i]], betas[,,i])
    #cat("beta_score",i, ":" , beta_score, "\n")
    differ_score[i] <- target_score - beta_score
  }
  # To avoid numerical issues, subtract the max score before exponentiating
  max_score <- max(differ_score)
  exp_scores <- exp(differ_score - max_score)
  
  # Normalize the exponentiated scores to get the importance weights
  importance_weights <- exp_scores / sum(exp_scores)
  #cat("importance_weights", ":" , importance_weights, "\n")
  ess_value <- 1 / sum(importance_weights^2)
  
  #cat("ess_value", ":" , ess_value, "\n")
  
  return(list(importance_weights = importance_weights, ess_value = ess_value))
}

importance_DAG_prev <- function(DAGs, s_betas){
  differ_score <- numeric()
  target_scores <- numeric()
  s_betas = DAG_scores[[i]]
  for (i in 1:length(DAGs)){
    target_score <- BiDAG::DAGscore(scoreParam, DAGs[[i]])
    beta_score <- s_betas[[i]]
    #cat("beta_score",i, ":" , beta_score, "\n")
    differ_score[i] <- target_score - beta_score
  }
  # To avoid numerical issues, subtract the max score before exponentiating
  max_score <- max(differ_score)
  exp_scores <- exp(differ_score - max_score)
  
  # Normalize the exponentiated scores to get the importance weights
  importance_weights <- exp_scores / sum(exp_scores)
  #cat("importance_weights", ":" , importance_weights, "\n")
  ess_value <- 1 / sum(importance_weights^2)
  #cat("ess_value", ":" , ess_value, "\n")
  return(list(importance_weights = importance_weights, ess_value = ess_value))
}


# get the probability of a parent node being the parent of a child node
get_probability <- function(betas, parentNode, childNode, DAG_Test){  
  incidence <- DAG_Test[[1]]
  get_probability <- numeric()
  if (parentNode != childNode) {
    # Check if parentNode is actually a parent in this DAG
    exp_beta <- exp(betas[parentNode,childNode, ])
    if (incidence[parentNode, childNode] == 1){
      get_probability <- exp_beta/(1 + exp_beta)
    }else{
      get_probability <- 1/(1 + exp_beta)
    }
  }else{
    # weight is not applicable for self (node cannot be its own parent)
    get_probability <- NA
  }
  return(get_probability)
}


####### Example 2: Generate random DAGs with 4 nodes

binary <- TRUE

if(binary){ 
  # For binary data
  n <- 8
  #scoreParam <- BiDAG::scoreparameters("bde", BiDAG::Asia[1:500,])
  scoreParam <- BiDAG::scoreparameters("bde", BiDAG::Asia)
  #BiDAG:::DAGcorescore(2, c(3,5), n = 8, param = scoreParam)
  DAG_Test <- list()
  DAG_Test[[1]] <- BiDAG::Asiamat
  targetDAG_score <- BiDAG::DAGscore(scoreParam, DAG_Test[[1]]) #-11105.32
}else{
  # For continuous data
  n <- 14 
  scoreParam <- BiDAG::scoreparameters("bge", BiDAG::Boston)
  itfit<-learnBN(scoreParam,algorithm="orderIter")
  maxEC<-getDAG(itfit,cp=FALSE)
  
  DAG_Test <- list()
  DAG_Test[[1]] <- maxEC
  target_beta <- calculateBetaScoresArray(DAG_Test, k = length(DAG_Test), n)
  targetDAG_score <- BiDAG::DAGscore(scoreParam, DAG_Test[[1]])  #> targetDAG_score[1] -20157.77
  #targetDAG_score <- 0.89 # Selected test value
  }

# Calculate move probabilities
prob1<-99
if(n>3){ prob1<-round(6*99*n/(n^2+10*n-24)) }
prob1<-prob1/100
moveprobs<-c(prob1,0.99-prob1,0.01)
moveprobs<-moveprobs/sum(moveprobs) # normalisation
if(!(length(moveprobs)==3)){print('Vector of move probabilities has the wrong length!')}
   
#start with empty graph or user defined graph  or from other algo. and learn the beta matrix
starting_dag <- list(matrix(c(0), nrow = n, ncol = n, byrow = TRUE))
# Replace diagonal elements with NA
betas_init <- calculateBetaScoresArray(starting_dag, k = 1 ,n)[,,1]
#base_score <- BiDAG::DAGscore(scoreParam, starting_dag[[1]])
base_score <- 0

# Threshold
# Define a threshold for individual differences
epsilon <- 0.8  # Boston 400, Asia 0.8
num_steps <- 5  # Number of steps to consider for fluctuation
# Define a threshold for the standard deviation of differences 
difference_threshold <- 0.01  # Boston 40, Asia 0.1
differences <- numeric()  # Initialize a vector to store the last 'num_steps' differences

# Define the starting order
permy <- c(1:n)

# Number of iterations
iter <- 100
weighted_betas <- list(betas_init)
averaged_beta_ess <- numeric()
ess_DAGs <- numeric()
DAG_scores <- list()
order_score <- numeric()

#Initialize for OrderMCMC
order_iter <-  100
order_stepsize <- 5

# Looping
for (i in 1:iter) {
  cat("permy for iter",i, ":" , permy, "\n")
  example <- orderMCMC_betas(n,startorder = permy,iterations = order_iter, 
                             betas = weighted_betas[[i]],
                             stepsave = order_stepsize, moveprobs)
  #print(permy)
  DAGs <- example[[1]]
  DAG_scores[[i]] <-example[[2]]
  
  ### Update beta matrix using importance sampling
  #is_results <- importance_DAG(DAGs = DAGs, betas = beta_values)
  is_results <- importance_DAG_prev(DAGs = DAGs, s_betas = unlist(DAG_scores[[i]]))
  
  # Vectorized multiplication of each matrix by its corresponding weight and then summing up the matrices
  weights <- is_results$importance_weights
  #cat("weights",i, ":" , weights, "\n")
  beta_values <- calculateBetaScoresArray(DAGs, k = length(DAGs) ,n) #a list of matrices
  
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
  
  if((i > num_steps+1)&& 
     (difference < epsilon) && 
     sd(differences[(i-num_steps): i]) < difference_threshold) {
    cat("Convergence achieved based on fluctuation criterion in iteration ", i, "\n")
    cat("Current order is ", unlist(example[[4]][length(example[[4]])]), "\n")
    break
  }
  

  #if((i > 1) && (ess_DAGs[i] > ess_DAGs[i-1])){
  if(difference > 300){
    #cat("ess_DAGs[i] < ess_DAGs[i]", permy, "\n")
    #permy <- unlist(example[[4]][1])
    #permy <- unlist(example[[4]][length(example[[4]])])
    #permy <- unlist(example[[4]][which.max(DAG_scores[[i]])[1]])
    permy <- permy
    order_score[i] <- example[[3]][length(example[[3]])]
  }else{
    ### Update the order for the next iteration
    #permy <- unlist(example[[4]][length(example[[4]])])
    #permy <- unlist(example[[4]][which.max(example[[3]])[1]])
    permy <- unlist(example[[4]][which.max(weights)[1]])
    order_score[i] <- example[[3]][which.max(weights)]
  }
}

print(ess_DAGs)

# Plotting the differences for the beta matrices
plot(#differences[-c(1:3)], 
     differences, 
     col = "blue",
     type = "b", main = "Differences for Beta_Matrix Per Iteration", 
     xlab = "Iteration", ylab = "Difference")
lines(c(1:iter), ess_DAGs, type = "o", col = "red")

final_betas <- weighted_betas[[i + 1]]


### Build a Consensus Graph from sampled DAGs
# Include an edge in the consensus graph if it appears in a significant number of DAGs
# Initialize matrices to store edge frequencies and cumulative beta values
edge_freq <- matrix(0, n, n)
#cumulative_beta <- matrix(0, n, n)

# Aggregate information from each DAG
for (i in 1:length(DAGs)) {
  dag <- DAGs[[i]]
  #beta_matrix <- beta_values[,,i]
  for (row in 1:n) {
    for (col in 1:n) {
      if (dag[row, col] == 1) {  # Check if there's an edge
        edge_freq[row, col] <- edge_freq[row, col] + 1
        #cumulative_beta[row, col] <- cumulative_beta[row, col] + beta_matrix[row, col]
      }
    }
  }
}

# Calculate average beta values
#average_beta <- cumulative_beta / edge_freq # beta values for each edge
#average_beta[is.nan(average_beta)] <- 0  # Handle division by zero
threshold_percentage <- seq(0, 1, by = 0.01)  # steps of 0.01

#Calculate TP and FP Rates
true_adj_mat <- BiDAG::Asiamat# the adjacency matrix of the true DAG

TP_rates <- numeric(length(threshold_percentage))
FP_rates <- numeric(length(threshold_percentage))

for (i in 1:length(threshold_percentage)) {
  threshold_edge <- length(DAGs) * threshold_percentage[i]
  consensus_graph <- edge_freq >= threshold_edge
  consensus_adj_mat <- ifelse(consensus_graph, 1, 0)
  
  TP <- sum(consensus_adj_mat == 1 & true_adj_mat == 1)
  FP <- sum(consensus_adj_mat == 1 & true_adj_mat == 0)
  TN <- sum(consensus_adj_mat == 0 & true_adj_mat == 0)
  FN <- sum(consensus_adj_mat == 0 & true_adj_mat == 1)
  
  TP_rates[i] <- TP / (TP + FN)
  FP_rates[i] <- FP / (FP + TN)
}

plot(FP_rates, TP_rates, type = "l", xlab = "False Positive Rate", ylab = "True Positive Rate", main = "Consensus Graph Performance")
abline(0, 1, col = "red", lty = 2)  # Diagonal line for reference

# Create an igraph graph from the adjacency matrix
graph <- graph_from_adjacency_matrix(consensus_adj_mat, mode = "directed", diag = FALSE)

# Plot the graph
plot(graph, 
     main = "Consensus Graph",
     edge.arrow.size = 0.5,
     vertex.color = "lightblue",
     vertex.size = 15,
     vertex.label.color = "black",
     vertex.label.cex = 0.8)


# To the skeleton space (from DAG to undirected graph)
skeleton_matrix <- (consensus_adj_mat | t(consensus_adj_mat))  # Union of the DAG and its transpose
diag(skeleton_matrix) <- 0  # Remove self-loops if present

library(igraph)
graph_skeleton <- graph_from_adjacency_matrix(skeleton_matrix, mode = "undirected", diag = FALSE)
#skeleton_object <- as.undirected(graph_skeleton, mode = "mutual")
# Plot the graph
plot(graph_skeleton, 
     main = "Skeleton Graph",
     edge.arrow.size = 0.5, vertex.color = "lightgreen",
     vertex.size = 15, vertex.label.color = "black", vertex.label.cex = 0.8)



#pattern space
#to benchpress

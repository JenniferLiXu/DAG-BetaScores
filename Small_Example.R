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
  sampledscore<-0
  for (child in 1:n){
      for (parent in parentnodes){
        sampledscore <- sampledscore + betas[parent, child]*incidence[parent, child]                      
      }
  }
  return(sampledscore)
}

#Importance Weights for the sampled DAGs

importance_DAG <- function(scores, targetDAG_score){
  # To avoid numerical issues, subtract the max score before exponentiating
  scores <- unlist(scores)
  #max_score <- max(scores)
  #exp_scores <- exp(scores - max_score)
  unormalized_weights  <-  targetDAG_score/scores
  
  # Normalize the exponentiated scores to get the importance weights
  importance_weights <- unormalized_weights / sum(unormalized_weights)
  ess_value <- 1 / sum(importance_weights^2)
  
  return(list(importance_weights = importance_weights, ess_value = ess_value))
}


####### Example 2: Generate random DAGs with 4 nodes

# For binary data
#scoreParam <- BiDAG::scoreparameters("bde", BiDAG::Asia)
#BiDAG:::DAGcorescore(2, c(3,5), n = 8, param = scoreParam)

n <- 14 
# For continuous data
scoreParam <- BiDAG::scoreparameters("bge", BiDAG::Boston)
#BiDAG:::DAGcorescore(2, c(3,5), n = 14, param = scoreParam)
itfit<-learnBN(scoreParam,algorithm="orderIter")
maxEC<-getDAG(itfit,cp=TRUE)

DAG_Test <- list()
DAG_Test[[1]] <- maxEC
#target_beta <- calculateBetaScoresArray(DAG_Test, k = length(DAG_Test), n)
targetDAG_score <- BiDAG::DAGscore(scoreParam, DAG_Test[[1]])  #> targetDAG_score[1] -20157.77
#targetDAG_score <- 0.89 # Selected test value

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


# Initialize

# Test Case 1 : Create a zero-matrix as our starting betas
#beta_matrix_init <- matrix(c(0), nrow = 14, ncol = 14, byrow = TRUE)

# Test Case 2 : Create a random matrix
# Assuming you want values to be drawn from a normal distribution
beta_matrix_init <- matrix(rnorm(n * n), nrow = n, ncol = n)

# Replace diagonal elements with NA
diag(beta_matrix_init) <- NA
# Convert the matrix to an array
betas_init <- array(beta_matrix_init, dim = c(14, 14))

# Define the starting order
permy <- c(1:14)

# Calculate move probabilities
prob1<-99
if(n>3){ prob1<-round(6*99*n/(n^2+10*n-24)) }
prob1<-prob1/100
moveprobs<-c(prob1,0.99-prob1,0.01)
moveprobs<-moveprobs/sum(moveprobs) # normalisation

# Number of iterations
iter <- 100

weighted_betas <- list(betas_init)
averaged_beta_ess <- numeric(iter)
ess_DAGs <- numeric(iter)
DAG_scores <- list()

# Threshold
epsilon <- 10  # Define a threshold for individual differences
num_steps <- 5  # Number of steps to consider for fluctuation
difference_threshold <- 1  # Define a threshold for the standard deviation of differences
differences <- numeric()  # Initialize a vector to store the last 'num_steps' differences

for (i in 1:iter) {
  example <- orderMCMC_betas(n,startorder = permy,iterations = 20, 
                             betas = weighted_betas[[i]],
                             stepsave = 1, moveprobs)
  DAGs <- example[[1]]
  DAG_scores[[i]] <-example[[2]]
  
  beta_values <- calculateBetaScoresArray(DAGs, k = length(DAGs) ,n) #a list of matrices
  
  ### Update beta matrix using importance sampling
  is_results <- importance_DAG(DAG_scores[[i]], targetDAG_score)
  
  # Vectorized multiplication of each matrix by its corresponding weight
  # and then summing up the matrices
  weights <- is_results$importance_weights
  #cat("weights",i, ":" , weights, "\n")
  weighted_betas[[i + 1]]  <- Reduce("+", lapply(1:length(weights), 
                                                 function(k) beta_values[,,k] * weights[k]))
  ess_DAGs[i] <-is_results$ess_value
  
  ### Stopping criteria for beta matrix convergence
  # Calculate the Frobenius norm of the difference between current and previous matrices
  current_beta_matrix <- weighted_betas[[i + 1]]
  previous_beta_matrix <- weighted_betas[[i]]
  # # Calculate difference only for non-NA pairs. Replace NA with 0 or other appropriate value
  current_beta_matrix[is.na(current_beta_matrix)] <- 0
  previous_beta_matrix[is.na(previous_beta_matrix)] <- 0
  
  # Calculate difference using Frobenius norm
  difference <- norm(current_beta_matrix - previous_beta_matrix, type = "F")
  # Update the differences vector
  differences <- c(differences, difference)
  
  if(all(differences < epsilon) && sd(differences[(i-num_steps): i]) < difference_threshold) {
    cat("Convergence achieved based on fluctuation criterion.\n")
    break
  }
  
  ### Update the order for the next iteration
  permy <- unlist(example[[4]][length(example[[4]])])
}


print(ess_DAGs)

# Plotting the differences for the beta matrices
plot(differences, type = "b", main = "Differences Per Iteration", 
     xlab = "Iteration", ylab = "Difference", col = "blue")

final_betas <- weighted_betas[[iter + 1]]


### Build a Consensus Graph from sampled DAGs

# Include an edge in the consensus graph if it appears in a significant number of DAGs
# Initialize matrices to store edge frequencies and cumulative beta values
edge_freq <- matrix(0, n, n)
cumulative_beta <- matrix(0, n, n)

# Aggregate information from each DAG
for (i in 1:length(DAGs)) {
  dag <- DAGs[[i]]
  beta_matrix <- beta_values[,,i]
  
  for (row in 1:n) {
    for (col in 1:n) {
      if (dag[row, col] == 1) {  # Check if there's an edge
        edge_freq[row, col] <- edge_freq[row, col] + 1
        cumulative_beta[row, col] <- cumulative_beta[row, col] + beta_matrix[row, col]
      }
    }
  }
}

# Calculate average beta values
average_beta <- cumulative_beta / edge_freq # beta values for each edge
average_beta[is.nan(average_beta)] <- 0  # Handle division by zero

# Define a threshold for including edges in the consensus graph
threshold <- length(DAGs) * 0.5  # Example: edge must appear in more than 50% of DAGs

# Build consensus graph
consensus_graph <- edge_freq >= threshold # The consensus graph 

# Convert the consensus graph to an adjacency matrix format expected by igraph
consensus_adj_mat <- ifelse(consensus_graph, 1, 0)

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




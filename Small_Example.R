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

# Choose whether to save the plots to png or pdf
saveoutput<-1 # 1 means save
filetypes<-c(".png",".pdf") # choose png or pdf for the figures
filetypey<-1 # we use png here  

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


#Importance Sampling
importance_beta <- function(target_beta, beta_values){
  if (is.null(dim(target_beta)) || is.null(dim(beta_values))) {
    stop("target_beta and beta_values must be matrices or arrays.")
  }
  
  n <- dim(target_beta)[1]
  
  # Initialize the matrix for averaged beta values 
  averaged_beta_matrix <- matrix(NA, n, n)  
  ess_matrix <- matrix(NA, n, n)  
  
  total_ess <- 0
  total_entries <- 0

  # Iterate over each matrix position
  for (row in 1:n) {
    for (col in 1:n) {
      if (row != col) {  # Exclude diagonal elements
        values <- beta_values[row, col, ]  # Extract values from the same position in all matrices
        weights <- (target_beta[row, col, 1])/values  # Define or compute weights for these values
        
        normalized_weights <- weights / sum(weights) # Normalized weights
        
        ess_value <- 1 / sum(normalized_weights^2)
        ess_matrix[row, col] <- ess_value
        total_ess <- total_ess + ess_value
        total_entries <- total_entries + 1
        
        # Compute the weighted average and assign it to the averaged matrix
        averaged_beta_matrix[row, col] <- sum(values * normalized_weights)
      }
    }
  }
  average_ess <- total_ess / total_entries
  return(list(averaged_beta_matrix = averaged_beta_matrix, ess_matrix = ess_matrix, average_ess = average_ess))
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
target_beta <- calculateBetaScoresArray(DAG_Test, k = length(DAG_Test), n)

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

# with weights defined by probabilities 
is_with_get_weight <- function(beta_values, DAG_Test){
  # Initialize the matrix for averaged beta values 
  averaged_beta_matrix <- matrix(0, n, n)  
  # Replace diagonal elements with NA
  diag(averaged_beta_matrix) <- NA
  values <-  numeric()
  weights <-  numeric()
  
  # Iterate over each matrix position
  for (row in 1:n) {
    for (col in 1:n) {
      if (row != col) {  # Exclude diagonal elements
        values <- beta_values[row, col, ]  # Extract values from the same position in all matrices
        weights <- get_probability(beta_values, row, col, DAG_Test)  # Define or compute weights for these values
        
        # Compute the weighted average and assign it to the averaged matrix
        averaged_beta_matrix[row, col] <- sum(values * weights)
      }
    }
  }
  return(averaged_beta_matrix)
}

# Initialize

# Create a zero-matrix as our starting betas
beta_matrix <- matrix(c(0), nrow = 14, ncol = 14, byrow = TRUE)
# Replace diagonal elements with NA
diag(beta_matrix) <- NA

# Convert the matrix to an array
betas <- array(beta_matrix, dim = c(14, 14))

# Define the starting order
permy <- c(1:14)

# Calculate move probabilities
prob1<-99
if(n>3){ prob1<-round(6*99*n/(n^2+10*n-24)) }
prob1<-prob1/100
moveprobs<-c(prob1,0.99-prob1,0.01)
moveprobs<-moveprobs/sum(moveprobs) # normalisation

# Number of iterations
iter <- 10

weighted_betas <- list(betas)
averaged_ess <- numeric(iter)

for (i in 1:iter) {
  example <- orderMCMC_betas(n,startorder = permy,iterations = 10, 
                             betas = weighted_betas[[i]],
                             stepsave = 1, moveprobs)
  DAGs <- example[[1]]
  
  # beta_values is a list of matrices from calculateBetaScoresArray
  beta_values <- calculateBetaScoresArray(DAGs, k = length(DAGs) ,n)
  
  # Update beta values using importance sampling
  is_results <- importance_beta(target_beta, beta_values)
  
  weighted_betas[[i + 1]] <- is_results$averaged_beta_matrix
  averaged_ess[i] <- is_results$average_ess
  
  #weighted_betas[[2]] <- is_with_get_weight(beta_values, DAG_Test)
  
  # Update the order for the next iteration
  permy <- unlist(example[[4]][length(example[[4]])])
}


print(averaged_ess)




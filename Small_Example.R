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

n <- 14
# Create a matrix
beta_matrix <- matrix(c(0), 
                      nrow = 14, ncol = 14, byrow = TRUE)

# Replace diagonal elements with NA
diag(beta_matrix) <- NA

# Convert the matrix to an array
betas <- array(beta_matrix, dim = c(14, 14))

permy <- c(1:14)

#DAG <-samplescore(n,betas,permy)
#calculateBetaScoresArray(sampledDAGs = DAG, n)


#Importance Sampling
importance_beta <- function(beta_values){
  # Initialize the matrix for averaged beta values 
  averaged_beta_matrix <- matrix(0, n, n)  
  # Replace diagonal elements with NA
  diag(averaged_beta_matrix) <- NA
  
  values = numeric()
  # Iterate over each matrix position
  for (row in 1:n) {
    for (col in 1:n) {
      if (row != col) {  # Exclude diagonal elements
        values <- beta_values[row, col, ]  # Extract values from the same position in all matrices
        weights <- (target_beta[row, col, 1])/values  # Define or compute weights for these values
        
        # Compute the weighted average and assign it to the averaged matrix
        averaged_beta_matrix[row, col] <- sum(values * weights / sum(weights))
      }
    }
  }
  return(averaged_beta_matrix)
}


####### Example 2: Generate random DAGs with 4 nodes

# For binary data
#scoreParam <- BiDAG::scoreparameters("bde", BiDAG::Asia)
#BiDAG:::DAGcorescore(2, c(3,5), n = 8, param = scoreParam)

# For continuous data
scoreParam <- BiDAG::scoreparameters("bge", BiDAG::Boston)
#BiDAG:::DAGcorescore(2, c(3,5), n = 14, param = scoreParam)
itfit<-learnBN(scoreParam,algorithm="orderIter")
maxEC<-getDAG(itfit,cp=TRUE)

DAG_Test <- list()
DAG_Test[[1]] <- maxEC
target_beta <- calculateBetaScoresArray(DAG_Test, k = length(DAG_Test), n)

prob1<-99
if(n>3){ prob1<-round(6*99*n/(n^2+10*n-24)) }
prob1<-prob1/100
moveprobs<-c(prob1,0.99-prob1,0.01)
moveprobs<-moveprobs/sum(moveprobs) # normalisation

# Number of iterations
iter <- 1

weighted_betas <- list()

for (round in 1:iter) {
  example <- orderMCMC_betas(n,startorder = permy,iterations = 10, betas,stepsave = 1, moveprobs)
  DAGs <- example[[1]]
  
  # beta_values is a list of matrices from calculateBetaScoresArray
  beta_values <- calculateBetaScoresArray(DAGs, k = length(DAGs) ,n)
  
  # use Importance Sampling to update values of betas
  betas <- importance_beta(beta_values)
  
  # store the weighted beta values
  weighted_betas[[round]] <- betas
}







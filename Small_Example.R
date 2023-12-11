#Creating a small DAG model in R to 
#test the proposed method for obtaining beta values

#install.packages("bnlearn")
#install.packages("BiDAG")
install.packages("gRain")

#Create a Small DAG
#library(bnlearn)
# Manually define a small DAG
#dag <- model2network("[A][B|A][C|A][D|B:C]")

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


install.packages("BiDAG")
install.packages("pcalg")
install.packages("igraph") 

library(BiDAG) 
library(pcalg)
library(igraph)

n <- 4
# Create a matrix
#beta_matrix <- matrix(c(0.5, 0.2, 0.3, 0.4), 
#                      nrow = 4, ncol = 4, byrow = TRUE)

# Replace diagonal elements with NA
#diag(beta_matrix) <- NA

# Convert the matrix to an array
#betas <- array(beta_matrix, dim = c(4, 4))

# Print the array
#print(betas)

permy <- c(3, 1, 4, 2)


#DAG <-samplescore(n,betas,permy)

####### Example1: Generate a random DAG with 4 nodes
#dag <- randDAG(4,2)

# Convert graphNEL object to igraph object
#igraphObj <- igraph::graph_from_graphnel(dag)

# Convert to adjacency matrix
#adjMatrix <- as_adjacency_matrix(igraphObj, sparse = FALSE)
#DAG <- array(adjMatrix,dim = dim(adjMatrix))

#calculateBetaScoresArray(sampledDAGs = DAG, n)


####### Example 2: Generate random DAGs with 4 nodes
prob1<-99
if(n>3){ prob1<-round(6*99*n/(n^2+10*n-24)) }
prob1<-prob1/100
moveprobs<-c(prob1,0.99-prob1,0.01)
moveprobs<-moveprobs/sum(moveprobs) # normalisation

example <- orderMCMC_betas(n,startorder = permy,iterations = 10, betas,stepsave = 1, moveprobs)
DAGs <- example[[1]]

beta_values <- calculateBetaScoresArray(DAGs,n)




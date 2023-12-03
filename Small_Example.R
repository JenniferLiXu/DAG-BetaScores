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

source('./edgerevandstructure/structurefns.R')
source('./edgerevandstructure/structureMCMC.R')
#source('./edgerevandstructure/newedgerevfns.R')
#source('./edgerevandstructure/newedgerevmove.R')

source('./orderandpartition/samplefns.R')
source('./scoring/combinations.R')
source('./scoring/scorefns.R')

# load a simple score proportional to the number of edges in the DAG
source('./scoring/numedgescore.R')

#Use the function to calculate beta values
source('./calculateBetaScoresArray.R')


install.packages("BiDAG")
install.packages("pcalg")
install.packages("igraph") 

library(BiDAG) 
library(pcalg)
library(igraph)

# Example: Generate a random DAG with 4 nodes
dag <- randDAG(4,2)

# Convert graphNEL object to igraph object
igraphObj <- igraph::graph_from_graphnel(dag)

# Convert to adjacency matrix
adjMatrix <- as_adjacency_matrix(igraphObj, sparse = FALSE)
DAG <- array(adjMatrix,dim = dim(adjMatrix))


calculateBetaScoresArray(sampledDAGs = DAG, n)


#Use an existing MCMC method to sample new DAGs and update β values.
# This step will involve resampling DAGs, updating beta values, and using MCMC methods


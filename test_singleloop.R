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
source('./CompareDAG_skeleton.R')
source('./BetaOrderSampler.R')

# Example: Generating a random dataset
set.seed(123)

n <- 7 # Define the number of nodes
scoreParam <- BiDAG::scoreparameters("bge", BiDAG::Boston[1:100,1:n])
#itfit<-learnBN(scoreParam,algorithm="orderIter")
#maxEC<-getDAG(itfit,cp=FALSE)

#posterior probabilities of edges by averaging over a sample of DAGs obtained via an MCMC scheme.
# In skeleton space
samplefit<-sampleBN(scoreParam, "order")
edgesposterior<-edgep(samplefit, pdag=FALSE, burnin=0.2)


# Move probabilities
prob1<-99
if(n>3){ prob1<-round(6*99*n/(n^2+10*n-24)) }
prob1<-prob1/100
moveprobs<-c(prob1,0.99-prob1,0.01)
moveprobs<-moveprobs/sum(moveprobs) # normalisation
if(!(length(moveprobs)==3)){print('Vector of move probabilities has the wrong length!')}
 
# Start with empty graph (OR user defined graph OR from other algo. and learn the beta matrix)
base_score <- 0 # Initialize the base score

# Initialization Parameters
num_iterations <- 200 # Total iterations

# Example 
results <- BetaOrderSampler(n = n, iter = num_iterations, order_iter = 100, 
                            order_stepsize = 10, moveprobs = moveprobs, 
                            edgesposterior = edgesposterior)


sum(results$acceptCount)

########## Plotting the differences for the beta matrices
plot(#differences[-c(1:3)],  
  results$diffBiDAGs[seq(1, length(results$diffBiDAGs), by = 10)],
  #diff_BiDAGs[seq(5, length(diff_BiDAGs), by = 5)], 
  col = "blue",
  type = "b", main = "Differences for Beta_Matrix Per Iteration", 
  xlab = "Iteration", ylab = "Difference")

########## For each edge, we plot a graph to see its changes in each iteration
# Setting up the plot
plot(NULL, xlim = c(1, num_iterations), ylim = range(results$edgeDifferences, na.rm = TRUE), #c(0,1),
     xlab = 'Iteration', ylab = 'Difference in Edge Probability', 
     main = 'Change in Edge Differences Over Iterations')

colors <- rainbow((n * (n - 1)) / 2)
legend_labels <- c()
color_index <- 1

# Plotting each edge difference over time
for (row in 1:(n - 1)) {
  for (col in (row + 1):n) {
    #lines(1:num_iterations, results$edgeDifferences[row, col, ], col = colors[color_index], type = 'l')
    lines(seq(1, num_iterations, by = 10), results$edgeDifferences[row, col, seq(1, num_iterations, by = 10)], col = colors[color_index], type = 'l')
    legend_labels <- c(legend_labels, paste('Edge', row, '-', col))
    color_index <- color_index + 1
  }
}

# Add a legend if needed
legend('topright', legend = legend_labels, col = colors, lty = 1, cex = 0.3)





# These graphs should be more or less similar with the BiDAG graph

# Also compare with partition MCMC


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


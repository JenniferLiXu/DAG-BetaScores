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
source('./orderandpartition_beta/orderMCMC_betas.R')
source('./orderandpartition_beta/orderscore_betas.R')
source('./orderandpartition_beta/partitionMCMC_betas.R')
source('./orderandpartition_beta/partitionscore_betas.R')
source('./orderandpartition_beta/partitionmoves_betas.R')
source('./orderandpartition_beta/samplefns-beta.R')
source('./orderandpartition_beta/samplefns_party-beta.R')

source('./importance_sampling_beta.R')
source('./calculateBetaScoresArray.R')
source('./CompareDAG_skeleton.R')
source('./BetaOrderSampler.R') #our method

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

# First choose the MCMC scheme
MCMCtype<-4 # 1 means standard structure, 2 with new edge reversal
# 3 means order MCMC, 4 means partition MCMC
# 5 means partition MCMC with new edge reversal


switch(as.character(MCMCtype),
       "1"={ # standard structure MCMC
         iterations<-100 #number of iterations in the chain
         moveprobs<-c(1) # having length 1 disallows the new edge reversal move
         if(!(length(moveprobs)==1)){print('Vector of move probabilities has the wrong length!')}
       },
       "2"={ # with new edge reversal
         iterations<-100 #number of iterations in the chain
         # Choose the probability of the different moves
         # 1 is standard structure MCMC (including the possibility to stay still [officially needed for convergence])
         # 2 is new edge reversal move
         moveprobs<-c(0.93,0.07)
         moveprobs<-moveprobs/sum(moveprobs) # normalisation
         if(!(length(moveprobs)==2)){print('Vector of move probabilities has the wrong length!')}
       },
       "3"={ # order MCMC
         iterations<-100 #number of iterations in the chain
         # Choose the probability of the different moves
         # 1 is swap any two elements
         # 2 is to only swap adjacent elements
         # 3 is to stay still (officially needed for convergence)
         prob1<-99
         if(n>3){ prob1<-round(6*99*n/(n^2+10*n-24)) }
         prob1<-prob1/100
         moveprobs<-c(prob1,0.99-prob1,0.01)
         moveprobs<-moveprobs/sum(moveprobs) # normalisation
         if(!(length(moveprobs)==3)){print('Vector of move probabilities has the wrong length!')}
       },
       "4"={ # partition MCMC
         iterations<-100 #number of iterations in the chain
         # Choose the probability of the different moves
         # 1 is swap two nodes from different partition elements
         # 2 is to only swap nodes from adjacent elements
         # 3 is to split or join partition elements
         # 4 is to move a single node
         # 5 is to stay still (officially needed for convergence)
         prob1start<-40/100
         prob1<-prob1start*100
         if(n>3){ prob1<-round(6*prob1*n/(n^2+10*n-24)) }
         prob1<-prob1/100
         prob2start<-99/100-prob1start
         prob2<-prob2start*100
         if(n>3){ prob2<-round(6*prob2*n/(n^2+10*n-24)) }
         prob2<-prob2/100
         moveprobs<-c(prob1,prob1start-prob1,prob2start-prob2,prob2,0.01)
         moveprobs<-moveprobs/sum(moveprobs) # normalisation
         if(!(length(moveprobs)==5)){print('Vector of move probabilities has the wrong length!')}
       },
       "5"={ # partition MCMC with edge reversal
         iterations<-100 #number of iterations in the chain
         # Choose the probability of the different moves
         # 1 is swap two nodes from different partition elements
         # 2 is to only swap nodes from adjacent elements
         # 3 is to split or join partition elements
         # 4 is to move a single node
         # 5 is to stay still (officially needed for convergence)
         # 6 is the new edge reversal
         prob1start<-37/100
         prob1<-prob1start*100
         if(n>3){ prob1<-round(6*prob1*n/(n^2+10*n-24)) }
         prob1<-prob1/100
         prob2start<-92/100-prob1start
         prob2<-prob2start*100
         if(n>3){ prob2<-round(6*prob2*n/(n^2+10*n-24)) }
         prob2<-prob2/100
         moveprobs<-c(prob1,prob1start-prob1,prob2start-prob2,prob2,0.01,0.07)
         moveprobs<-moveprobs/sum(moveprobs) # normalisation
         if(!(length(moveprobs)==6)){print('Vector of move probabilities has the wrong length!')}
       },
       {# if none is chosen, we have a problem
         print('Not implemented')
       })


 
# Start with empty graph (OR user defined graph OR from other algo. and learn the beta matrix)
base_score <- 0 # Initialize the base score

# Initialization Parameters
num_iterations <- 500 # Total iterations

# Example 
switch(as.character(MCMCtype),
       "3"={ # # order MCMC
         results <- BetaOrderSampler(n = n, iter = num_iterations, order_iter = 100, 
                                     order_stepsize = 10, moveprobs = moveprobs, 
                                     edgesposterior = edgesposterior )
       },
       "4"={ # partition MCMC
         results <- BetaPartitionSampler(n = n, iter = num_iterations, party_iter = 100, 
                                         party_stepsize = 10, moveprobs = moveprobs, 
                                         edgesposterior = edgesposterior )
       }
       )


sum(results$acceptCount)

########## Plotting the differences between our matrix and the one from BiDAG
plot(#differences[-c(1:3)],  
  results$diffBiDAGs[seq(1, length(results$diffBiDAGs), by = 1)],
  #diff_BiDAGs[seq(5, length(diff_BiDAGs), by = 5)], 
  col = "blue",
  type = "b", main = "Differences between Matrices Per Iteration", 
  xlab = "Iteration", ylab = "Difference")

########## For each edge, we plot a graph to see its changes in each iteration
# # Setting up the plot
# plot(NULL, xlim = c(1, length(results$DAGs)), ylim = c(0,1),#range(results$edge_prob, na.rm = TRUE), #c(0,1),#
#      xlab = 'Iteration', ylab = 'Difference in Edge Probability', 
#      main = 'Change in Edge Differences Over Iterations')
# 
# colors <- rainbow((n * (n - 1)) / 2)
# legend_labels <- c()
# color_index <- 1
# 
# # Plotting each edge difference over time
# for (row in 1:(n - 1)) {
#   for (col in (row + 1):n) {
#     #lines(1:num_iterations, results$edgeDifferences[row, col, ], col = colors[color_index], type = 'l')
#     lines(seq(1, length(results$diffBiDAGs), by = 1), 
#           results$edgeDifferences[row, col, seq(1, length(results$diffBiDAGs), by = 1)], 
#           col = colors[color_index], type = 'l')
#     legend_labels <- c(legend_labels, paste('Edge', row, '-', col))
#     color_index <- color_index + 1
#   }
# }
# # Add a legend if needed
# legend('topright', legend = legend_labels, col = colors, lty = 1, cex = 0.3)



plot(NULL, xlim = c(1, length(results$diffBiDAGs)), ylim = c(0, 1),  # Adjust ylim based on actual range of your data if needed
     xlab = 'Iteration', ylab = 'Edge Probability', 
     main = 'Edge Probability Over Iterations')

# Generate enough colors for all edges, considering a fully connected directed graph without self-loops
colors <- rainbow(n * (n - 1))
legend_labels <- vector("character", length = n * (n - 1))
color_index <- 1

# Plotting each edge over time, skipping diagonals
for (row in 1:n) {
  for (col in 1:n) {
    if (row != col) {  # Skip diagonals
      lines(seq(1, length(results$diffBiDAGs), by = 1), 
            results$edge_prob[row, col, seq(1, length(results$diffBiDAGs), by = 1)], 
            col = colors[color_index], type = 'l')
      legend_labels[color_index] <- paste('Edge', row, '->', col)
      color_index <- color_index + 1
    }
  }
}

# Adding legend
# Note: Displaying a legend for a large number of edges might not be practical
# Consider using a subset or interactive plotting for detailed inspection
if (color_index <= 10) {  # Arbitrary threshold to avoid clutter
  legend("topright", legend = legend_labels, col = colors, lty = 1, cex = 0.5)
}



orderfitBoston100<-orderMCMC(scoreParam,plus1=TRUE,chainout=TRUE)
plotpedges(orderfitBoston100, cutoff = 0, pdag=FALSE)

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


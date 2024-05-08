#Creating a small DAG model in R to 
#test the proposed method for obtaining beta values

#install.packages("bnlearn")
# install.packages("BiDAG")
# install.packages("pcalg")
# install.packages("igraph") 
# install.packages("combinat")

.libPaths("/cluster/home/xuwli/R/x86_64-pc-linux-gnu-library/4.3") 
library(BiDAG) 
library(pcalg)
library(igraph)
library(bnlearn)
library(ggplot2)
library(combinat)

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

source('./calculateBetaScoresArray.R')
source('./compareDAG_skeleton.R')
source('./BetaOrderSampler.R') #our method
source('./BetaOrderSampler_OneDAG.R') 
source('./BetaOrderSampler_gibbs.R') 


# Example: Generating a random dataset
n <- 4 # Define the number of nodes,
scoreParam <- BiDAG::scoreparameters("bge", BiDAG::Boston[1:250,1:n])
#itfit<-learnBN(scoreParam,algorithm="orderIter")
#maxEC<-getDAG(itfit,cp=FALSE)

#posterior probabilities of edges by averaging over a sample of DAGs obtained via an MCMC scheme.
# In skeleton space
samplefit<-sampleBN(scoreParam, "orderIter")
edgesposterior<-edgep(samplefit, pdag=FALSE, burnin=0.2)

edgesposterior <-  matrix(c(0,0.003747,0.37539,0.536539,
                      0.0006246,0, 0.179262,0.004372,
                      0.624609,0.820737,0,0.06121,
                      0.131792,0.003123,0.019987,0), nrow = 4, ncol = 4, byrow = TRUE)

#          crim          zn      indus        chas
# crim  0.0000000000 0.003747658 0.37539038 0.536539663
# zn    0.0006246096 0.000000000 0.17926296 0.004372267
# indus 0.6246096190 0.820737039 0.00000000 0.061211743
# chas  0.1317926296 0.003123048 0.01998751 0.000000000

MCMCtype<-3 # 1 means standard structure, 2 with new edge reversal
# 3 means order MCMC, 4 means partition MCMC, 5 means partition MCMC with new edge reversal
 
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
       {# if none is chosen, we have a problem
         print('Not implemented')
       })


# Initialization Parameters
# num_iterations <- 1000# Total iterations

# Example 
starttime_model<-proc.time() 

switch(as.character(MCMCtype),
       "3"={ # # order MCMC
         num_iterations <- 3000 
         set.seed(100) 
         results_seed123 <- BetaOrderSampler(n = n, iteration = num_iterations, order_iter = 100, 
                                                    # order = list(c(4,2,1,3)),
                                                    order_stepsize = 100, moveprobs = moveprobs,
                                                    edgesposterior = edgesposterior )
         num_iterations <- 3000
         set.seed(100) 
         results_seed300 <- BetaOrderSampler_gibbs_ver2(n = n, iteration = num_iterations, order_iter = 100, 
                                                    # order = list(c(4,2,1,3)),
                                                    order_stepsize = 10, moveprobs = moveprobs, 
                                                    edgesposterior = edgesposterior )
         num_iterations <- 3000
         set.seed(100) 
         results_seed500 <- BetaOrderSampler_gibbs_ver3(n = n, iteration = num_iterations, order_iter = 100, 
                                             # order = list(c(4,2,1,3)),
                                             order_stepsize = 100, moveprobs = moveprobs, 
                                             edgesposterior = edgesposterior )
         # num_iterations <- 10000
         # set.seed(100) 
         # results_seed1000 <- BetaOrderSampler(n = n, iteration = num_iterations, order_iter = 100, 
         #                                     # order = list(c(4,2,1,3)),
         #                                     order_stepsize = 100, moveprobs = moveprobs, 
         #                                     edgesposterior = edgesposterior )
       },
       "4"={ # partition MCMC
         set.seed(123) 
         results_seed123 <- BetaPartitionSampler(n = n, iter = num_iterations, party_iter = 100, 
                                         party_stepsize = 100, moveprobs = moveprobs, 
                                         edgesposterior = edgesposterior )
         set.seed(100)
         results_seed100<- BetaPartitionSampler(n = n, iter = num_iterations, party_iter = 100, 
                                 party_stepsize = 100, moveprobs = moveprobs, 
                                 edgesposterior = edgesposterior )
       }
       )

endtime_model<-proc.time()
endtime_model<-endtime_model-starttime_model
print(endtime_model)

# sum(results_seed123$acceptCount)/num_iterations
# sum(results_seed100$acceptCount)/num_iterations
# num_iterations <- 1e4 
# user   system  elapsed 
# 1307.209    5.505 1319.328 (calculateBetaScoresArray)
# 767.520   8.882 786.310    (calculateBetaScoresArray_hash)

# num_iterations <- 1e3 
# user  system elapsed 
#  143.562   0.833 147.043 (calculateBetaScoresArray)
#  72.534   0.695  72.791 (calculateBetaScoresArray_hash)
#  33.980   0.408  35.066  (calculateBetaScoresArray_hash new version)


# ########## Plotting the differences between our matrix and the one from BiDAG
# plot(#differences[-c(1:3)],  
#   results_seed123$diffBiDAGs[seq(1, length(results_seed123$diffBiDAGs), by = 1)],
#   #diff_BiDAGs[seq(5, length(diff_BiDAGs), by = 5)], 
#   col = "blue",
#   type = "b", main = "Differences between Matrices Per Iteration", 
#   xlab = "Iteration", ylab = "Difference")
# 
# 
# plot(NULL, xlim = c(1, length(results_seed123$diffBiDAGs)), ylim = c(0, 1),  # Adjust ylim based on actual range of your data if needed
#      xlab = 'Iteration', ylab = 'Edge Probability', 
#      main = 'Edge Probability Over Iterations')
# 
# # Generate enough colors for all edges, considering a fully connected directed graph without self-loops
# colors <- rainbow(n * (n - 1))
# legend_labels <- vector("character", length = n * (n - 1))
# color_index <- 1
# 
# # Plotting each edge over time, skipping diagonals
# for (row in 1:n) {
#   for (col in 1:n) {
#     if (row != col) {  # Skip diagonals
#       lines(seq(1, length(results_seed123$diffBiDAGs), by = 1), 
#             results_seed123$edge_prob[row, col, seq(1, length(results_seed123$diffBiDAGs), by = 1)], 
#             col = colors[color_index], type = 'l')
#       legend_labels[color_index] <- paste('Edge', row, '->', col)
#       color_index <- color_index + 1
#     }
#   }
# }
# 
# # Adding legend
# # Note: Displaying a legend for a large number of edges might not be practical
# # Consider using a subset or interactive plotting for detailed inspection
# if (color_index <= 200) {  # Arbitrary threshold to avoid clutter
#   legend("topright", legend = legend_labels, col = colors, lty = 1, cex = 0.5)
# }


starting_mat <- matrix(1, nrow = n, ncol = n)
diag(starting_mat) <- 0
set.seed(100)
orderfitBoston100<-orderMCMC(scoreParam, iterations = num_iterations , MAP = FALSE,chainout=TRUE, startspace = starting_mat)
set.seed(123)
orderfitBoston123<-orderMCMC(scoreParam, iterations = num_iterations , MAP = FALSE,chainout=TRUE, startspace = starting_mat)
plotpedges(orderfitBoston123, cutoff = 0, pdag=FALSE)

# These graphs should be more or less similar with the BiDAG graph

# Also compare with partition MCMC

#results_seed1 and 2 for 1e5 iterration
pedges <-  list()
pedges[[1]] <-  edgep(orderfitBoston100, pdag=FALSE)
pedges[[2]] <- edgep(orderfitBoston123, pdag=FALSE)
# pdf("0506plot_order_betaOrder_OneDAGs.pdf")
plot_order_Order <- plotpcor(pedges, xlab="run1", ylab="run2",printedges=TRUE, main = paste("Iteration", num_iterations) )
cat("order_Order: ",plot_order_Order$MSE, plot_order_Order$R2 , "\n")

pedges_comp123 <-  list()
pedges_comp123[[1]] <- edgep(orderfitBoston100, pdag=FALSE)
pedges_comp123[[2]] <- results_seed123$edge_prob[,,length(results_seed123$edge_prob[1,1,])]
dimnames(pedges_comp123[[2]]) <- dimnames(pedges_comp123[[1]])
plot_order_betaOrder <- plotpcor(pedges_comp123, xlab="run1", ylab="run2",printedges=TRUE, main = "OrderMCMC and BetaSampler-5000iter-weigthed")
cat("order_betaOrder123: ",plot_order_betaOrder$MSE, plot_order_betaOrder$R2 , "\n")

pedges_comp300 <-  list()
pedges_comp300[[1]] <- edgep(orderfitBoston100, pdag=FALSE)
pedges_comp300[[2]] <- results_seed300$edge_prob[,,length(results_seed300$edge_prob[1,1,])]
dimnames(pedges_comp300[[2]]) <- dimnames(pedges_comp300[[1]])
plot_order_betaOrder <- plotpcor(pedges_comp300, xlab="run1", ylab="run2",printedges=TRUE, main = "OrderMCMC and GibbBetaSampler-5000iter")
cat("order_betaOrder300: ",plot_order_betaOrder$MSE, plot_order_betaOrder$R2 , "\n")

pedges_comp500 <-  list()
pedges_comp500[[1]] <- edgep(orderfitBoston100, pdag=FALSE)
pedges_comp500[[2]] <- results_seed500$edge_prob[,,length(results_seed500$edge_prob[1,1,])]
dimnames(pedges_comp500[[2]]) <- dimnames(pedges_comp500[[1]])
plot_order_betaOrder <- plotpcor(pedges_comp500, xlab="run1", ylab="run2",printedges=TRUE, main = "OrderMCMC and GibbBetaSampler-5000iter-ver3")
cat("order_betaOrder500: ",plot_order_betaOrder$MSE, plot_order_betaOrder$R2 , "\n")

# pedges_comp1000 <-  list()
# pedges_comp1000[[1]] <- edgep(orderfitBoston100, pdag=FALSE)
# pedges_comp1000[[2]] <- results_seed1000$edge_prob[,,length(results_seed1000$edge_prob[1,1,])]
# dimnames(pedges_comp1000[[2]]) <- dimnames(pedges_comp1000[[1]])
# plot_order_betaOrder <- plotpcor(pedges_comp1000, xlab="run1", ylab="run2",printedges=TRUE, main = "OrderMCMC and BetaSampler-1e4iter")
# cat("order_betaOrder1000: ",plot_order_betaOrder$MSE, plot_order_betaOrder$R2 , "\n")

pedges_seed <-  list()
pedges_seed[[1]] <- results_seed123$edge_prob[,,length(results_seed123$edge_prob[1,1,])]
pedges_seed[[2]] <- results_seed300$edge_prob[,,length(results_seed300$edge_prob[1,1,])]
dimnames(pedges_seed[[1]]) <- dimnames(pedges_comp[[1]])
dimnames(pedges_seed[[2]]) <- dimnames(pedges_comp[[1]])
plot_order_betaOrder_seed <- plotpcor(pedges_seed, xlab="run1", ylab="run2",printedges=TRUE, main = "Comparison betw. BetaSamplers")
cat("order_betaOrder_seed: ",plot_order_betaOrder_seed$MSE, plot_order_betaOrder_seed$R2 , "\n")

dev.off()

########### Graphing part ############

### Build a Consensus Graph from sampled DAGs
# Include an edge in the consensus graph if it appears in a significant number of DAGs
# Initialize matrices to store edge frequencies and cumulative beta values
# 
# 
# edge_freq_perc <- DAG_total/length(DAG)
# 
# consensus_graph <- edge_freq_perc >= 0.22
# consensus_adj_mat <- ifelse(consensus_graph, 1, 0)
# 
# # Create an igraph graph from the adjacency matrix
# graph <- graph_from_adjacency_matrix(consensus_adj_mat, mode = "undirected", diag = FALSE)
# 
# # Plot the graph
# plot(graph, 
#      main = "Consensus Graph",
#      edge.arrow.size = 0.5,
#      vertex.color = "lightblue",
#      vertex.size = 15,
#      vertex.label.color = "black",
#      vertex.label.cex = 0.8)




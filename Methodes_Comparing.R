install.packages("pcalg")
install.packages("MASS")
library(pcalg)
library(MASS)

source('./SNOrderSampler_OneDAG.R')

p=30
n=100

set.seed(123)
generate_data <- function(p, n) {
  # Generate a random DAG
  g <- randomDAG(p, prob = 0.2)  # Adjust edge probability as needed
  
  # Get adjacency matrix from graphNEL object
  adj_matrix <- as(g, "matrix")
  
  # Assign random weights to edges, avoiding values close to zero
  weights <- runif(sum(adj_matrix != 0), min = 0.25, max = 1)
  weights <- ifelse(runif(length(weights)) > 0.5, weights, -weights)  # Randomly assign signs
  
  # Apply weights to the adjacency matrix
  w_adjMat <- matrix(0, nrow = p, ncol = p)
  w_adjMat[adj_matrix != 0] <- weights
  
  # Simulate data from the SEM with additive Gaussian noise
  error_terms <- matrix(rnorm(n * p), nrow = n)
  X <- error_terms  # Start with the noise
  
  for (i in 1:p) {
    # Update the data based on the SEM
    X[i, ] <- X[i, ] + X[i, ] %*% w_adjMat
  }
  
  return(list(data = X, graph = g, weights = w_adjMat))
}


data_100 <- generate_data(p = 30, n = 100)

data_1000 <- generate_data(p = 30, n = 1000)

n <- 30 # Define the number of nodes,

n <- 5
scoreParam <- BiDAG::scoreparameters("bge", data_1000$data[,1:n])
samplefit<-sampleBN(scoreParam, "order")
edgesposterior<-edgep(samplefit, pdag=FALSE, burnin=0.2)

edgesposterior <- as(data_100$graph, "matrix")

num_iterations <- 1000
set.seed(123)
results_seed_test_100 <- BetaOrderSampler_OneDAG_test(n = n, iteration = num_iterations, order_iter = 100,
                                                # order = list(c(4,2,1,3)),
                                                order_stepsize = 10, moveprobs = moveprobs,
                                                edgesposterior = edgesposterior )

num_iterations <- 500
set.seed(100)
results_seed_test_100_1 <- BetaOrderSampler_OneDAG_test(n = n, iteration = num_iterations, order_iter = 100,
                                                      # order = list(c(4,2,1,3)),
                                                      order_stepsize = 10, moveprobs = moveprobs,
                                                      edgesposterior = edgesposterior )

pedges_comp_test <-  list() 
# pedges_comp_test[[1]] <- edgep(orderfitBoston100, pdag=FALSE)
pedges_comp_test[[1]] <- edgesposterior
pedges_comp_test[[2]] <- results_seed_test_100$edge_prob[,,length(results_seed_test_100$edge_prob[1,1,])]
dimnames(pedges_comp_test[[2]]) <- dimnames(pedges_comp_test[[1]])
plot_order_betaOrder <- plotpcor(pedges_comp_test, xlab="run1", ylab="run2",printedges=TRUE, main = "OrderMCMC & SN -1 DAG")
cat("order_betaOrder123: ",plot_order_betaOrder$MSE, plot_order_betaOrder$R2 , "\n")
plot(results_seed_test_100$edge_prob[3,2,])
sum(results_seed_test_100$acceptCount)

pedges_comp_test_1 <-  list() 
# pedges_comp_test[[1]] <- edgep(orderfitBoston100, pdag=FALSE)
pedges_comp_test_1[[1]] <- edgesposterior
pedges_comp_test_1[[2]] <- results_seed_test_100_1$edge_prob[,,length(results_seed_test_100_1$edge_prob[1,1,])]
# dimnames(pedges_comp123[[2]]) <- dimnames(pedges_comp123[[1]])
plot_order_betaOrder_1 <- plotpcor(pedges_comp_test_1, xlab="run1", ylab="run2",printedges=TRUE, main = "OrderMCMC & SN -1 DAG")
cat("order_betaOrder123: ",plot_order_betaOrder_1$MSE, plot_order_betaOrder_1$R2 , "\n")
plot(results_seed_test_100$edge_prob[3,2,])
sum(results_seed_test_100$acceptCount)


pedges_SN_seed <- list() 
pedges_SN_seed[[1]] <- results_seed_test_100$edge_prob[,,length(results_seed_test_100$edge_prob[1,1,])]
pedges_SN_seed[[2]] <- results_seed_test_100_1$edge_prob[,,length(results_seed_test_100_1$edge_prob[1,1,])]
dimnames(pedges_SN_seed[[1]]) <- dimnames(pedges_comp_test[[1]])
plot_SN_seed <- plotpcor(pedges_SN_seed, xlab="run1", ylab="run2",printedges=TRUE, main = "OrderMCMC & SN -1 DAG")
cat("SN_seed: ",plot_SN_seed$MSE, plot_SN_seed$R2 , "\n")
plot(results_seed_test_100_1$edge_prob[3,2,])
sum(results_seed_test_100$acceptCount)

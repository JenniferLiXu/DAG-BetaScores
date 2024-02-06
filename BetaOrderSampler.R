# n: number of nodes in the DAG.
# iter : the total number of iterations to run.
# order_iter and order_stepsize are parameters for the Order MCMC process.
# moveprobs : a vector of probabilities for moves in the MCMC process.
# base_score : the initial score, set to 0 by default.
# starting_dag : an optional parameter to start with a user-defined DAG. If not provided, it starts with an empty DAG.
# edgesposterior : the matrix of edge probabilities in the target (true) DAG.
# burrnin : burn in iteration, by default 20 percent 

# Output: a list containing the final DAGs, edge differences over each iteration, 
#         ESS values, acceptance counts, and differences from the target DAG.

BetaOrderSampler <- function(n, iter, order_iter, order = list(seq_len(n)), 
                             order_stepsize, moveprobs, base_score = 0, 
                             starting_dag = NULL, 
                             edgesposterior, burnin = 0.2 ) {
  # Initialize starting DAG if not provided
  if (is.null(starting_dag)) {
    starting_dag <- list(matrix(0, nrow = n, ncol = n))
  }
  
  # Initialize variables
  betas_init <- calculateBetaScoresArray(starting_dag, k = 1, n)[,,1]
  weighted_betas <- list(betas_init)
  weights <- list(1)
  ess_DAGs <- numeric()
  DAG <- starting_dag
  skeleton_DAG <- starting_dag
  total_DAG <- matrix(0, nrow = n, ncol = n)
  single_DAG <- starting_dag
  count_accept <- numeric()
  diff_BiDAGs <- numeric()
  edge_diff_over_time <- array(0, dim = c(n, n, iter))
  burin_iter <- floor(burnin*iter)
  
  # Looping through iterations
  for (i in 1:iter) {
    beta_prev <- weighted_betas[[i]]
    order_prev <- order[[i]]
    
    # Sampling orders with OrderMCMC
    example <- orderMCMC_betas(n,startorder = order_prev ,iterations = order_iter, 
                               betas = beta_prev,
                               stepsave = order_stepsize, moveprobs) # run the Order MCMC code
    
    #Store the last order from the chain
    #permy <- unlist(example[[4]][length(example[[4]])])
    
    max_order_score <- max(unlist(example[[3]]))
    permy <- example[[4]][which(example[[3]] == max_order_score)][[1]]
    
    #  Sample 10 DAGs using the order with highest score from OrderMCMC
    sampled_DAGs <- lapply(1:20, function(x) samplescore(n, beta_prev, permy, base_score))
    
    # Extracting incidence matrices and log scores
    incidence_matrices <-lapply(sampled_DAGs, function(dag) dag$incidence) 
    incidence_logscore<-lapply(sampled_DAGs, function(dag) dag$logscore) 
    
    # Update beta matrix using importance sampling
    is_results <- importance_DAG(DAGs = incidence_matrices, score_under_betas = incidence_logscore)
    weights_proposed <- is_results$importance_weights
    ess_DAGs[i] <-is_results$ess_value
    single_DAG[[i+1]] <- is_results$compress_dag
    
    # Update beta matrix and calculate proposed log score
    beta_values <- calculateBetaScoresArray(incidence_matrices, k = length(incidence_matrices) ,n) 
    weighted_betas_proposed <- Reduce("+", lapply(1:length(weights_proposed), 
                                                  function(k) beta_values[,,k] * weights_proposed[k]))
    
    # Log score of Proposed DAG set under the previous beta_matrix 
    proposed_logscore <- Reduce("+", lapply(1:length(weights_proposed), 
                                            function(k) incidence_logscore[[k]] * weights_proposed[k]))
    
    # Calculate current log score
    #current_logscore <- calculate_DAG_score(DAG = DAGs[[i]], permy = order[[i]], 
    #                                        weights = weights[[i]],betas = weighted_betas_proposed)
    # Log score of previous DAG set under the proposed beta matrix
    current_logscore <- calculate_singleDAG_score(DAG = single_DAG[[i]], permy = order[[i]],
                                                  betas = weighted_betas_proposed)
    
    # Acceptance ratio
    ratio <-  exp(proposed_logscore - current_logscore)
    
    # Metropolis-Hastings acceptance step
    if((i == 1 )|| (runif(1) < ratio)){ 
      DAG[[i]] <- single_DAG[[i+1]]
      skeleton_DAG[[i]] <- combineProbabilities(single_DAG[[i+1]])
      weighted_betas[[i+1]] <- weighted_betas_proposed
      order[[i+1]] <- permy
      weights[[i+1]] <- weights_proposed
      count_accept[i] <- 1 # Accept 
    }else{
      DAG[[i]] <- DAG[[i-1]]
      skeleton_DAG[[i]] <- skeleton_DAG[[i-1]]
      weighted_betas[[i+1]] <- weighted_betas[[i]]
      order[[i+1]] <- order[[i]]
      weights[[i+1]] <- weights[[i]]
      count_accept[i] <- 0 # Reject
    }
    
    if (length(DAG) > burin_iter) {
      # Get the last few matrices
      total_DAG <- total_DAG + skeleton_DAG[[i]]
      #sum_matrix <- Reduce("+", last_matrices) # Sum the matrices
      # Update diff_BiDAG
      diff_mat <- CompareDAG_skeleton(total_DAG/(i-burin_iter), edgesposterior)
    }else{
      sum_matrix <- Reduce("+", skeleton_DAG[1:length(DAG)])
      # Update diff_BiDAG
      diff_mat <- CompareDAG_skeleton(sum_matrix/i, edgesposterior)
    }
    
    edge_diff_over_time[,,i] <- diff_mat # Store the difference per edge
    diff_BiDAG <-norm(diff_mat,type = "F") # Store the difference of the DAG matrix
    diff_BiDAGs <- c(diff_BiDAGs, diff_BiDAG)
    
  }
  
  # Return the results
  return(list(DAGs = DAG, edgeDifferences = edge_diff_over_time, essValues = ess_DAGs, acceptCount = count_accept, diffBiDAGs = diff_BiDAGs))
}

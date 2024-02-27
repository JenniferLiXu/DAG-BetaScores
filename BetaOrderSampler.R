# n: number of nodes in the DAG.
# iter : the total number of iterations to run.
# order_iter and order_stepsize are parameters for the Order MCMC process.
# moveprobs : a vector of probabilities for moves in the MCMC process.
# base_score : the initial score, set to 0 by default.
# starting_dag : an optional parameter to start with a user-defined DAG. If not provided, it starts with an empty DAG.
# edgesposterior : the matrix of edge probabilities in the target (true) DAG.
# burrnin : burn in iteration, by default 10 percent 

# Output: a list containing the final DAGs, edge differences over each iteration, 
#         ESS values, acceptance counts, and differences from the target DAG.

BetaOrderSampler <- function(n, iter, order_iter, order = NULL, 
                             order_stepsize, moveprobs, base_score = 0, 
                             starting_dag = NULL, betas_init = NULL, skeleton = FALSE,
                             edgesposterior, burnin = 0.2 ) {
  # Initialize starting DAG if not provided
  if (is.null(starting_dag)) {
    starting_dag <- list(matrix(0, nrow = n, ncol = n))
  }
  # Initialize beta matrix
  if (is.null(betas_init)) {
    betas_init <- calculateBetaScoresArray(starting_dag, k = 1, n)[,,1]
  }
  # Initialize order
  if (is.null(order)) {
    order <-  list(seq_len(n))
  }
  
  # Initialize variables
  weighted_betas <- list(betas_init)
  ess_DAGs <- numeric()
  count_accept <- numeric()
  diff_BiDAGs <- numeric()
  
  DAG <- starting_dag
  single_DAG <- starting_dag
  total_DAG <- matrix(0, nrow = n, ncol = n)
  edge_over_time <- array(0, dim = c(n, n, iter))
  edge_diff_over_time <- array(0, dim = c(n, n, iter))
  burin_iter <- floor(burnin*iter)
  prev_weight <- 1
  
  # Looping through iterations
  for (i in 1:iter) {
    beta_prev <- weighted_betas[[i]]
    order_prev <- order[[i]]
    
    # Sampling orders with OrderMCMC
    example <- orderMCMC_betas(n,startorder = order_prev ,iterations = order_iter, 
                               betas = beta_prev,
                               stepsave = order_stepsize, moveprobs) # run the Order MCMC code
    
    #Store the last order from the chain
    permy <- unlist(example[[4]][length(example[[4]])])
    
    #  Sample 30 DAGs using the last sampled order from OrderMCMC
    sampled_DAGs <- lapply(1:10, function(x) samplescore(n, beta_prev, permy, base_score))
    
    # Extracting incidence matrices and log scores
    incidence_matrices <-lapply(sampled_DAGs, function(dag) dag$incidence) 
    incidence_logscore<-lapply(sampled_DAGs, function(dag) dag$logscore) 
    
    # Update beta matrix using importance sampling
    is_results <- importance_DAG(DAGs = incidence_matrices, score_under_betas = incidence_logscore)
    weights_proposed <- is_results$importance_weights
    ess_DAGs[i] <-is_results$ess_value
    
    # Sample one DAG from our sampled DAGs, using normalised weights as the probability
    represent_sample <- sample(c(1:length(sampled_DAGs)),size = 1, prob = weights_proposed)
    represent_DAG <- incidence_matrices[[represent_sample]]
    represent_weight <- weights_proposed[represent_sample]
    
    # Update beta matrix using the weights from sampled DAGs
    beta_values <- calculateBetaScoresArray(incidence_matrices, k = length(incidence_matrices) ,n) 
    weighted_betas_proposed <- Reduce("+", lapply(1:length(weights_proposed), 
                                                  function(k) beta_values[,,k] * weights_proposed[k]))
    
    # Log score of Proposed DAG set under the previous beta_matrix 
    #proposed_logscore <- Reduce("+", lapply(1:length(weights_proposed), function(k) incidence_logscore[[k]] * weights_proposed[k]))
    proposed_logscore <- calculate_DAG_score(DAG_list = list(represent_DAG),permy = permy, weights = c(1), betas =  beta_prev)
    # Calculate current log score(DAG from last iteration under the new beta)
    #current_logs3core <- calculate_DAG_score(DAG_list = list(single_DAG[[i]]),permy = order_prev, weights = c(1) ,betas = weighted_betas_proposed)
    current_logscore <- calculate_DAG_score(DAG_list = list(DAG[[i]]),permy = order_prev, weights = c(1), betas = weighted_betas_proposed)
    
    # Acceptance ratio
    ratio <-  exp(proposed_logscore - current_logscore)
    
    # Metropolis-Hastings acceptance step
    if(runif(1) < ratio){ 
      DAG[[i+1]] <- represent_DAG
      single_DAG[[i+1]] <- is_results$compress_dag
      weighted_betas[[i+1]] <- weighted_betas_proposed
      order[[i+1]] <- permy
      count_accept[i] <- 1 # Accept 
      prev_weight <- represent_weight
      cat("ratio:", ratio, "\n")
    }else{
      DAG[[i+1]] <- DAG[[i]]
      single_DAG[[i+1]] <- single_DAG[[i]]
      weighted_betas[[i+1]] <- weighted_betas[[i]]
      order[[i+1]] <- order[[i]]
      count_accept[i] <- 0 # Reject
    }
    
    if (length(single_DAG)-1 > burin_iter) {
      total_DAG <- total_DAG + single_DAG[[i+1]]
      current_mat <- total_DAG/(i - burin_iter) # Average the edges of DAGs after burn in part
      diff_mat <- CompareDAG(current_mat, edgesposterior)
    }else{
      sum_matrix <- Reduce("+", single_DAG[1:length(single_DAG)])
      current_mat <- sum_matrix/i
      diff_mat <- CompareDAG(current_mat, edgesposterior)
    }
    
    edge_over_time[,,i] <- current_mat
    edge_diff_over_time[,,i] <- diff_mat # Store the difference per edge
    diff_BiDAG <-norm(diff_mat,type = "F") # Store the difference of the DAG matrix
    diff_BiDAGs <- c(diff_BiDAGs, diff_BiDAG)
    
  }
  
  # Return the results
  return(list(DAGs = DAG[-c(1:burin_iter)], 
              edgeDifferences = edge_diff_over_time[,,-c(1:burin_iter)], 
              edge_prob = edge_over_time[,,-c(1:burin_iter)], 
              essValues = ess_DAGs[-c(1:burin_iter)], 
              acceptCount = count_accept[-c(1:burin_iter)], 
              diffBiDAGs = diff_BiDAGs[-c(1:burin_iter)]))
}


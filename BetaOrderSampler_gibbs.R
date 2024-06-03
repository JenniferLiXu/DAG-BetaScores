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


# In this version: samples from each iteration are treated with equal weights
BetaOrderSampler_gibbs <- function(n, iteration, order_iter , order = NULL, 
                                   order_stepsize , moveprobs, base_score = 0, 
                                   starting_dag = NULL, betas_init = NULL, skeleton = FALSE,
                                   edgesposterior, burnin = 0.3 ) {
  elements <- c(1:n) # Define set of elements
  # Initialize order
  if (is.null(order)) {
    order <-  list(seq_len(n))
  }
  # Initialize starting DAG if not provided
  if (is.null(starting_dag)) {
    starting_dag <-list(samplescore(n, betas = matrix(0, nrow = n, ncol = n), order[[1]])$incidence)
  }
  
  # Initialize beta matrix
  if (is.null(betas_init)) {
    calcultion_betas_init <- calculateBetaScoresArray_hash(starting_dag, k = 1, n, base_score = 0)
    betas_init <- calcultion_betas_init$allBetaScores[,,1]
    # base_score <- calcultion_betas_init$target_DAG_score
    base_score <- 0
  }
  
  # Initialize variables
  weighted_betas <- list(betas_init)
  order_prev <- order
  ess_DAGs <- numeric()
  count_accept <- 0
  diff_BiDAGs <- numeric()
  burin_iter <- floor(burnin*iteration)
  iter <- iteration + burin_iter
  
  compress_DAG <- starting_dag
  total_DAG <- matrix(0, nrow = n, ncol = n)
  edge_over_time <- array(0, dim = c(n, n, iter))
  edge_diff_over_time <- array(0, dim = c(n, n, iter))

  # Looping through iterations
  for (i in 1:iter) {
    beta_prev <- weighted_betas[[i]]
    # Sampling orders with OrderMCMC
    example <- orderMCMC_betas(n,startorder = order_prev[[length(order_prev)]] ,iterations = order_iter, betas = beta_prev,
                               stepsave = order_stepsize, moveprobs) # run the Order MCMC code
    
    proposed_orders <- example[[4]][-1]
    # print(proposed_orders)
    
    sampled_DAGs_fromOrder <- DAGs_from_order(order_list = proposed_orders, nr_sample = 10, beta_matrix = beta_prev)
    incidence_matrices <- sampled_DAGs_fromOrder$incidence # List of DAGs sampled under previous beta
    incidence_logscore <- sampled_DAGs_fromOrder$logscore # List of logscores of new sampled DAGs under previous beta
    
    # Update beta matrix using the weights from sampled DAGs
    calculation_beta_values <- calculateBetaScoresArray_hash(incidence_matrices, k = length(incidence_matrices) ,n, base_score = base_score)
    BiDAGscore_propose_list <- calculation_beta_values$target_DAG_score
    beta_values <- calculation_beta_values$allBetaScores
    
    # Update beta matrix using importance sampling
    is_results <- importance_DAG(DAGs = incidence_matrices, score_under_betas = incidence_logscore, target_scores = BiDAGscore_propose_list)
    weights_proposed <- is_results$importance_weights # normalised weights under old beta
    
    # cat("incidence_logscore: ",unlist(incidence_logscore), "\n")
    # cat("is_weights_proposed: ",weights_proposed, "\n")
    # cat("ess: ",is_results$ess_value, "\n")
    #New beta matrix using the normalised weights
    weighted_betas_proposed <- Reduce("+", lapply(1:length(weights_proposed), function(k) beta_values[,,k] * weights_proposed[k]))
    
    compress_DAG[[i+1]] <- is_results$compress_dag

    order_prev <- sampled_DAGs_fromOrder$order
    
    if (length(compress_DAG) > burin_iter) {
      ess_DAGs[i] <- is_results$ess_value
      # print(ess_DAGs[i])
      # print(compress_DAG[[i+1]])
      if (ess_DAGs[i] > 0.4){
        weighted_betas[[i+1]] <- weighted_betas_proposed
        # print(compress_DAG[[i+1]])
        # print(order_prev)
        count_accept <- count_accept + 1
        total_DAG <- total_DAG + compress_DAG[[i+1]]
        current_mat <- total_DAG/count_accept # Average the edges of DAGs after burn in part
      }else{
        current_mat <- current_mat
        weighted_betas[[i+1]] <- weighted_betas[[i]]
      }
    }else{
      weighted_betas[[i+1]] <- weighted_betas_proposed
      ess_DAGs[i] <- 0
      sum_matrix <- Reduce("+", compress_DAG[1:length(compress_DAG)])
      current_mat <- sum_matrix/i
    }
    edge_over_time[,,i] <- current_mat
  }
  
  # Return the results
  return(list(
    edge_prob = edge_over_time[,,-c(1:burin_iter)], 
    essValues = ess_DAGs[-c(1:burin_iter)], 
    acceptCount = count_accept, 
    betas = weighted_betas[-c(1:burin_iter)]
    # ,totalscore_of_DAGs = totalscore_of_DAGs[-c(1:burin_iter)]
  )
  )
}

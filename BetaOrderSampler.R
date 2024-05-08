# n: number of nodes in the DAG.
# iter : the total number of iterations to run.
# order_iter and order_stepsize are parameters for the Order MCMC process.
# moveprobs : a vector of probabilities for moves in the MCMC process.
# base_score : the initial score
# starting_dag : an optional parameter to start with a user-defined DAG. If not provided, it starts with an empty DAG.
# edgesposterior : the matrix of edge probabilities in the target (true) DAG.
# burrnin : burn in iteration, by default 10 percent 

# Output: a list containing the final DAGs, edge differences over each iteration, 
#         ESS values, acceptance counts, and differences from the target DAG.

# iteration = num_iterations
# order_iter = 100
# order_stepsize = 100
# moveprobs = moveprobs
# edgesposterior = edgesposterior
 
BetaOrderSampler <- function(n, iteration, order_iter, order = NULL, 
                             order_stepsize, moveprobs, base_score = 0, 
                             starting_dag = NULL, betas_init = NULL, skeleton = FALSE,
                             edgesposterior, burnin = 0.4) {
  
  elements <- c(1:n) # Define set of elements
  permutations <- permn(elements) # Generate all permutations
  nr_sample <- 10 # Number of generated DAGs under one order
  
  # Initialize order
  if (is.null(order)) {
    order <-  list(seq_len(n))
  }
  # Initialize starting DAG if not provided
  if (is.null(starting_dag)) {
    starting_dag <-list(samplescore(n, matrix(0, nrow = n, ncol = n), order[[1]])$incidence)
  }
  # Initialize beta matrix
  if (is.null(betas_init)) {
    calcultion_betas_init <- calculateBetaScoresArray_hash(starting_dag, k = 1, n, base_score = 0)
    betas_init <- calcultion_betas_init$allBetaScores[,,1]
    # base_score <- calcultion_betas_init$target_DAG_score
    base_score <- 0
    # Idea: subtract the highest order score
    prev_logscore_list <- sapply(1:length(permutations), function(k) {
      sum(orderscore_betas(n, c(1:n), betas_init, permutations[[k]]))
      })
    totalscore_orders_prev <- calculate_final_score(prev_logscore_list, operation = "sum")
  }
  
  # Initialize variables
  weighted_betas <- list(betas_init)
  ess_DAGs <- numeric()
  count_accept <- numeric()
  diff_BiDAGs <- numeric()
  burin_iter <- floor(burnin*iteration)
  iter <- iteration + burin_iter
  # Initialize format for outputs
  compress_DAG <- starting_dag
  total_DAG <- matrix(0, nrow = n, ncol = n)
  edge_over_time <- array(0, dim = c(n, n, iter))
  edge_diff_over_time <- array(0, dim = c(n, n, iter))
  
  init_sampled_DAGs <- lapply(1:nr_sample, function(x) samplescore(n, betas_init, order[[1]]))
  DAG_prev <- lapply(init_sampled_DAGs, function(dag) dag$incidence) 
  # incidence_logscore<-lapply(init_sampled_DAGs, function(dag) dag$logscore) 
  BiDAGscore_prev_list <- calculateBetaScoresArray_hash(DAG_prev, k = length(DAG_prev), n, base_score)$target_DAG_score
  
  # Looping through iterations
  for (i in 1:iter) {
    beta_prev <- weighted_betas[[i]]
    order_prev <- order[[i]]
    
    example <- orderMCMC_betas(n,startorder = order_prev ,iterations = order_iter, betas = beta_prev,
                              stepsave = order_stepsize, moveprobs)  # Sampling orders with OrderMCMC

    permy <- unlist(example[[4]][length(example[[4]])])#Store the last order from the chain
    
    #  Sample nr_sample DAGs using the last sampled order from OrderMCMC and old beta matrix
    sampled_DAGs <- lapply(1:nr_sample, function(x) samplescore(n, beta_prev, permy))
    
    # Extracting incidence matrices and log scores (old beta matrix)
    incidence_matrices <-lapply(sampled_DAGs, function(dag) dag$incidence) 
    incidence_logscore<-lapply(sampled_DAGs, function(dag) dag$logscore) # scores of new DAGs under old beta
    
    # Update beta matrix using the weights from sampled DAGs
    calculation_beta_values <- calculateBetaScoresArray_hash(incidence_matrices, k = length(incidence_matrices) ,n, base_score = base_score)
    BiDAGscore_propose <- calculation_beta_values$target_DAG_score
    beta_values <- calculation_beta_values$allBetaScores
    
    # Update beta matrix using importance sampling
    is_results <- importance_DAG(DAGs = incidence_matrices, score_under_betas = incidence_logscore, target_scores = BiDAGscore_propose)
    weights_proposed <- is_results$importance_weights # normalised weights under old beta
    
    # weighted_betas_proposed <- calculte_beta_matrix(beta_values, weigths = weights_proposed)
    #New beta matrix using the normalised weights
    weighted_betas_proposed <- Reduce("+", lapply(1:length(weights_proposed),
                                                  function(k) beta_values[,,k] * weights_proposed[k]))
    # weighted_betas_proposed <- Reduce("+", lapply(1:length(weights_proposed),
    #                                               function(k) beta_values[,,k] /nr_sample))

    proposed_logscore_list <- sapply(1:length(permutations), 
                                     function(k) sum(orderscore_betas(n, c(1:n), weighted_betas_proposed, permutations[[k]])))
    proposed_totalscore_orders <- calculate_final_score(proposed_logscore_list, operation = "sum")
    
    # Log score of old DAG set under the new beta
    oDAGnBeta_logscore_list <- calculate_DAG_score(DAG_list = DAG_prev ,permy = order_prev, weights = NULL ,
                                              betas = weighted_betas_proposed, base_score = base_score)
    
    # Acceptance ratio
    # Test
    wB <- calculate_final_score(is_results$log_diff, operation = "mean")
    # wB <- mean(is_results$log_diff)
    inv_wA <- calculate_final_score(oDAGnBeta_logscore_list - BiDAGscore_prev_list, operation = "mean")
    # inv_wA <- mean(oDAGnBeta_logscore_list - BiDAGscore_prev_list)
    log_ratio <- wB + inv_wA  + totalscore_orders_prev - proposed_totalscore_orders
    ratio <- exp(log_ratio)
    
    # Metropolis-Hastings acceptance step
    if(runif(1) < ratio){ 
      compress_DAG[[i+1]] <- is_results$compress_dag
      # compress_DAG[[i+1]] <- Reduce("+", lapply(1:length(incidence_matrices), function(k) incidence_matrices[[k]]/nr_sample))
      weighted_betas[[i+1]] <- weighted_betas_proposed
      order[[i+1]] <- permy
      count_accept[i] <- 1 # Accept 
      BiDAGscore_prev_list <- BiDAGscore_propose
      DAG_prev <- incidence_matrices
      totalscore_orders_prev <- proposed_totalscore_orders
      # print(permy)
      # print(log_ratio)
      # cat("wB",wB, "inv_wA", inv_wA, "\n")
    }else{
      compress_DAG[[i+1]] <- compress_DAG[[i]]
      weighted_betas[[i+1]] <- weighted_betas[[i]]
      order[[i+1]] <- order[[i]]
      count_accept[i] <- 0 # Reject
    }
    
    if (length(compress_DAG)-1 > burin_iter) {
      total_DAG <- total_DAG + compress_DAG[[i+1]]
      current_mat <- total_DAG/(i - burin_iter) # Average the edges of DAGs after burn in part
      diff_mat <- CompareDAG(current_mat, edgesposterior)
    }else{
      sum_matrix <- Reduce("+", compress_DAG[1:length(compress_DAG)])
      current_mat <- sum_matrix/i
      diff_mat <- CompareDAG(current_mat, edgesposterior)
    }
    
    edge_over_time[,,i] <- current_mat
    edge_diff_over_time[,,i] <- diff_mat # Store the difference per edge
    diff_BiDAG <-norm(diff_mat,type = "F") # Store the difference of the DAG matrix
    diff_BiDAGs <- c(diff_BiDAGs, diff_BiDAG)
  }
  
  # Return the results
  return(list(# DAGs = DAG[-c(1:burin_iter)], 
              edgeDifferences = edge_diff_over_time[,,-c(1:burin_iter)], 
              edge_prob = edge_over_time[,,-c(1:burin_iter)], 
              essValues = ess_DAGs[-c(1:burin_iter)], 
              acceptCount = count_accept[-c(1:burin_iter)], 
              betas = weighted_betas[[iter+1]],
              diffBiDAGs = diff_BiDAGs[-c(1:burin_iter)]
              # ,totalscore_of_DAGs = totalscore_of_DAGs[-c(1:burin_iter)]
              )
         )
}


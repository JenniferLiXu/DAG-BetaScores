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

BetaPartitionSampler <- function(n, iteration, party_iter, order = NULL, party = NULL,
                             party_stepsize, moveprobs, base_score = 0, 
                             starting_dag = NULL, betas_init = NULL, skeleton = FALSE,
                             edgesposterior, burnin = 0.4 ) {
  
  nr_sample <- 1 # Number of generated DAGs under one order
  # Initialize starting DAG if not provided
  if (is.null(starting_dag)) {
    starting_dag <- list(matrix(0, nrow = n, ncol = n))
  }
  # Initialize beta matrix
  if (is.null(betas_init)) {
    calcultion_betas_init <- calculateSN_ScoresArray_hash(starting_dag, k = 1, n, base_score = 0)
    betas_init <- calcultion_betas_init$allBetaScores[,,1]
    base_score <- 0
  }
  # Initialize order
  if (is.null(order)) {
    order <-  list(seq_len(n))
  }
  # Initialize partition
  if (is.null(party)) {
    party <- list(c(n)) # a starting partition - c(n) gives the empty DAG
  }

  # Initialize variables
  weighted_betas <- list(betas_init)
  ess_DAGs <- numeric()
  count_accept <- numeric()
  diff_BiDAGs <- numeric()
  burin_iter <- floor(burnin*iteration)
  iter <- iteration + burin_iter
  
  compress_DAG <- starting_dag
  total_DAG <- matrix(0, nrow = n, ncol = n)
  edge_over_time <- array(0, dim = c(n, n, iter))
  edge_diff_over_time <- array(0, dim = c(n, n, iter))
  
  init_sampled_DAGs <- lapply(1:nr_sample, function(x) samplescore_partition(n, betas = betas_init, 
                                                                             permy = order[[1]], party = c(n), 
                                                                             posy = parttolist(n,c(n)), 
                                                                             base_score))
  
  DAG_prev <- lapply(init_sampled_DAGs, function(dag) dag$incidence) 
  BiDAGscore_prev <- calculateSN_ScoresArray_hash(DAG_prev, k = length(DAG_prev), n, base_score)$target_DAG_score
  
  # Looping through iterations
  for (i in 1:iter) {
    beta_prev <- weighted_betas[[i]]
    order_prev <- order[[i]]
    party_prev <- party[[i]]
    
    # Sampling orders with PartitionMCMC
    example <- partitionMCMC_betas(n,startpermy = order_prev,startparty = party_prev,iterations = party_iter ,
                             stepsave = party_stepsize ,betas = beta_prev,moveprobs) 
    
    #Store the last order from the chain
    permy <- unlist(example[[4]][length(example[[4]])])
    partition <- unlist(example[[5]][length(example[[5]])])
    partyscore_prev <- unlist(example[[3]][length(example[[3]])])

    #  Sample 30 DAGs using the last sampled order from PartitionMCMC
    sampled_DAGs <- lapply(1:nr_sample, function(x) samplescore_partition(n, betas = beta_prev, 
                                                                   permy = permy, party = partition, 
                                                                   posy = parttolist(n,partition), 
                                                                   base_score))
    
    # Extracting incidence matrices and log scores
    incidence_matrices <-lapply(sampled_DAGs, function(dag) dag$incidence) 
    incidence_logscore<-lapply(sampled_DAGs, function(dag) dag$logscore) 
    
    # Ordres and Partitions used for the set of DAGs
    proposed_orders <- example[[4]][-1]
    proposed_party <- example[[5]][-1]
    
    # Update beta matrix using the weights from sampled DAGs
    calculation_beta_values <- calculateSN_ScoresArray_hash(incidence_matrices, k = length(incidence_matrices) ,n, base_score)
    BiDAGscore_propose <- calculation_beta_values$target_DAG_score
    beta_values <- calculation_beta_values$allBetaScores
    
    # Update beta matrix using importance sampling
    is_results <- importance_DAG(DAGs = incidence_matrices, score_under_betas = incidence_logscore, target_scores = calculation_beta_values$target_DAG_score)
    weights_proposed <- is_results$importance_weights
    
    #New beta matrix using the normalised weights
    weighted_betas_proposed <- Reduce("+", lapply(1:length(weights_proposed), 
                                                  function(k) beta_values[,,k] * weights_proposed[k]))
    
    partyscore_prop <- sum(partitionscore(n,c(1:n), weighted_betas_proposed, order_prev, party_prev, posy = parttolist(n,party_prev)))
    
    # # Log score of old DAG set under the new beta
    oDAGnBeta_logscore <- calculate_DAG_score(DAG_list = DAG_prev, permy = order_prev, weights = NULL,
                                              betas = weighted_betas_proposed, base_score, party = party_prev, posy = parttolist(n,party_prev))
    
    partyscore_prev <- unlist(example[[3]][-1]) # score of new order, previous beta
    # print(orderscore_prev)
    partyscore_prop <- sapply(1:length(proposed_orders), function(k) sum(partitionscore(n, c(1:n), weighted_betas_proposed, 
                                                                                        proposed_orders[[k]],party = proposed_party[[k]], posy = parttolist(n,proposed_party[[k]]))))
    set_new <- unlist(partyscore_prop) - unlist(partyscore_prev)
    new <- calculate_final_score_mean(set_new)
    # Acceptance ratio
    wB <- is_results$log_diff
    inv_wA <- oDAGnBeta_logscore - BiDAGscore_prev
    
    log_ratio <- wB + inv_wA  - new
    ratio <- exp(log_ratio)
    
    # Metropolis-Hastings acceptance step
    if(runif(1) < ratio){ 
      compress_DAG[[i+1]] <- is_results$compress_dag
      weighted_betas[[i+1]] <- weighted_betas_proposed
      order[[i+1]] <- permy
      party[[i+1]] <- partition
      count_accept[i] <- 1 # Accept 
      BiDAGscore_prev <- BiDAGscore_propose
      DAG_prev <- incidence_matrices
    }else{
      compress_DAG[[i+1]] <- compress_DAG[[i]]
      weighted_betas[[i+1]] <- weighted_betas[[i]]
      order[[i+1]] <- order[[i]]
      party[[i+1]] <- party[[i]]
      count_accept[i] <- 0 # Reject
      # cat("ratio:", ratio, "\n")
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
  return(list(
              edgeDifferences = edge_diff_over_time[,,-c(1:burin_iter)], 
              edge_prob = edge_over_time[,,-c(1:burin_iter)], 
              essValues = ess_DAGs[-c(1:burin_iter)], 
              acceptCount = count_accept[-c(1:burin_iter)], 
              betas = weighted_betas[[iter+1]],
              diffBiDAGs = diff_BiDAGs[-c(1:burin_iter)]))
}
 

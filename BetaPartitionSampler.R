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
# 
# iteration = 10
# party_iter = 100
# party_stepsize = 10
# burnin = 0.2
BetaPartitionSampler <- function(n, iteration, party_iter, order = NULL, party = NULL,
                             party_stepsize, moveprobs, base_score = 0, 
                             starting_dag = NULL, betas_init = NULL, skeleton = FALSE,
                             edgesposterior, burnin = 0.2 ) {
  # Initialize starting DAG if not provided
  if (is.null(starting_dag)) {
    starting_dag <- list(matrix(0, nrow = n, ncol = n))
  }
  # Initialize beta matrix
  if (is.null(betas_init)) {
    calcultion_betas_init <- calculateBetaScoresArray_hash(starting_dag, k = 1, n)
    betas_init <- calcultion_betas_init$allBetaScores[,,1]
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
  
  DAG <- starting_dag
  compress_DAG <- starting_dag
  total_DAG <- matrix(0, nrow = n, ncol = n)
  edge_over_time <- array(0, dim = c(n, n, iter))
  edge_diff_over_time <- array(0, dim = c(n, n, iter))
  prev_weight <- 1
  
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
    sampled_DAGs <- lapply(1:30, function(x) samplescore_partition(n, betas = beta_prev, 
                                                                   permy = permy, party = partition, 
                                                                   posy = parttolist(n,partition), 
                                                                   base_score))
    
    # Extracting incidence matrices and log scores
    incidence_matrices <-lapply(sampled_DAGs, function(dag) dag$incidence) 
    incidence_logscore<-lapply(sampled_DAGs, function(dag) dag$logscore) 
    
    # Update beta matrix using the weights from sampled DAGs
    calculation_beta_values <- calculateBetaScoresArray_hash(incidence_matrices, k = length(incidence_matrices) ,n)
    beta_values <- calculation_beta_values$allBetaScores
    
    # Update beta matrix using importance sampling
    is_results <- importance_DAG(DAGs = incidence_matrices, score_under_betas = incidence_logscore, target_scores = calculation_beta_values$target_DAG_score)
    weights_proposed <- is_results$importance_weights
    
    # if(calculation_beta_values$count_DAGcore !=0){print(calculation_beta_values$count_DAGcore)}
    weighted_betas_proposed <- Reduce("+", lapply(1:length(weights_proposed), 
                                                  function(k) beta_values[,,k] * weights_proposed[k]))
   
    # Sample one DAG from our sampled DAGs, using normalised weights as the probability
    represent_sample <- sample(c(1:length(sampled_DAGs)),size = 1, prob = weights_proposed)
    represent_DAG <- incidence_matrices[[represent_sample]]
    represent_weight <- weights_proposed[represent_sample]
    
    partyscore_prop <- sum(partitionscore(n,c(1:n), weighted_betas_proposed, order_prev, party_prev))
    
    # Log score of Proposed DAG set under the previous beta_matrix 
    proposed_logscore <- Reduce("+", lapply(1:length(weights_proposed), function(k) incidence_logscore[[k]] * weights_proposed[k]))
    #proposed_logscore <- calculate_DAG_score(DAG_list = list(represent_DAG),permy = permy, weights = c(1) ,
    #                                         betas =  beta_prev,
    #                                         party = partition, posy = parttolist(n,partition)) 
    # Calculate current log score(DAG from last iteration under the new beta)
    current_logscore <- calculate_DAG_score(DAG_list = list(compress_DAG[[i]]),permy = order[[i]], weights = c(1) ,betas = weighted_betas_proposed)
    #current_logscore <- calculate_DAG_score(DAG_list = list(DAG[[i]]), permy = order_prev, weights = c(1) ,
    #                                        betas = weighted_betas_proposed,
    #                                        party = party_prev, posy = parttolist(n,party_prev))
    
    # Acceptance ratio
    ratio <- exp(proposed_logscore + partyscore_prop - current_logscore - partyscore_prev)
    
    # Metropolis-Hastings acceptance step
    if(runif(1) < ratio){ 
      DAG[[i+1]] <- represent_DAG
      compress_DAG[[i+1]] <- is_results$compress_dag
      weighted_betas[[i+1]] <- weighted_betas_proposed
      order[[i+1]] <- permy
      party[[i+1]] <- partition
      count_accept[i] <- 1 # Accept 
    }else{
      DAG[[i+1]] <- DAG[[i]]
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
  return(list(DAGs = DAG[-c(1:burin_iter)], 
              edgeDifferences = edge_diff_over_time[,,-c(1:burin_iter)], 
              edge_prob = edge_over_time[,,-c(1:burin_iter)], 
              essValues = ess_DAGs[-c(1:burin_iter)], 
              acceptCount = count_accept[-c(1:burin_iter)], 
              betas = weighted_betas[[iter+1]],
              diffBiDAGs = diff_BiDAGs[-c(1:burin_iter)]))
}
 

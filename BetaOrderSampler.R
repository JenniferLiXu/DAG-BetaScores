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
  nr_sample <- 40 # Number of generated DAGs under one order
  
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
  is_weights <- list(rep(1/nr_sample, nr_sample))

  # Looping through iterations
  for (i in 1:iter) {
    ess <- 0
    ess_threshold <- 0.5
    while (ess < ess_threshold) {
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
      ess <- is_results$ess_value
      
      # if (ess < ess_threshold) {
      #   cat("ESS too low (", ess, ") in iteration", i, ", rerunning iteration...\n")
      # }
    }
    weights_proposed <- is_results$importance_weights # normalised weights under old beta
    
    #New beta matrix using the normalised weights
    # weighted_betas_proposed <- Reduce("+", lapply(1:length(weights_proposed),
    #                                               function(k) beta_values[,,k] * weights_proposed[k]))
    weighted_betas_proposed <- Reduce("+", lapply(1:length(weights_proposed), function(k) beta_values[,,k] /nr_sample))

    proposed_logscore_list <- sapply(1:length(permutations), 
                                     function(k) sum(orderscore_betas(n, c(1:n), weighted_betas_proposed, permutations[[k]])))
    proposed_totalscore_orders <- calculate_final_score(proposed_logscore_list, operation = "sum")
    
    # orderscore_prev <- unlist(example[[3]][length(example[[3]])]) #new order, previous beta
    
    # Log score of old DAG set under the new beta
    oDAGnBeta_logscore_list <- calculate_DAG_score(DAG_list = DAG_prev ,permy = order_prev, weights = NULL ,
                                              betas = weighted_betas_proposed, base_score = base_score)
    
    is_results_prev <- importance_DAG(DAGs = DAG_prev, score_under_betas = oDAGnBeta_logscore_list, target_scores = BiDAGscore_prev_list)
    prev_weight <-is_results_prev$importance_weights
    
    modify_weight <- unlist(lapply(1:length(weights_proposed), function(k) is_weights[[i]][k]/weights_proposed[k]))
    normolized_modify_weight <- modify_weight/sum(modify_weight)
    # modify_incidence_logscore <- lapply(1:length(normolized_modify_weight), function(k) incidence_logscore[[k]]*normolized_modify_weight[k]) 
    # modify_incidence_logscore <- lapply(1:length(modify_weight), function(k) incidence_logscore[[k]]*normolized_modify_weight[k]) 
    # modify_res<- importance_DAG(DAGs = incidence_matrices, score_under_betas = incidence_logscore, target_scores = BiDAGscore_propose)
    modify_log_diff <- unlist(lapply(1:length(normolized_modify_weight), function(k) is_results$log_diff[k]*prev_weight[k]))
    # log_diff <- unlist(lapply(1:length(is_weights[[i]]), function(k) (oDAGnBeta_logscore_list - BiDAGscore_prev_list)[k]*is_weights[[i]][k]))
    modify_log_diff_prev <- unlist(lapply(1:length(prev_weight), function(k) is_results_prev$log_diff[k]*weights_proposed[k]))
    
    # Acceptance ratio
    # Test
    # wB <- calculate_final_score(modify_log_diff, operation = "sum")
    # inv_wA <- calculate_final_score(oDAGnBeta_logscore_list - BiDAGscore_prev_list, operation = "mean")

    wB <- sum(modify_log_diff)
    inv_wA <- sum(modify_log_diff_prev)
    # inv_wA <- sum(unlist(lapply(1:length(weights_proposed), function(k) unlist(oDAGnBeta_logscore_list - BiDAGscore_prev_list)[k]*weights_proposed[k])))
    
    # print(inv_wA)
    # print(sum( unlist(lapply(1:length(weights_proposed), function(k) (oDAGnBeta_logscore_list - BiDAGscore_prev_list)[[k]]*weights_proposed[k]))))
    cat("wB +  inv_wA",wB +  inv_wA, "\n")
    log_ratio <- wB +  inv_wA + totalscore_orders_prev - proposed_totalscore_orders
    ratio <- exp(log_ratio)
    
    # Metropolis-Hastings acceptance step
    if(runif(1) < ratio){ 
      compress_DAG[[i+1]] <- is_results$compress_dag
      # compress_DAG[[i+1]] <- Reduce("+", lapply(1:length(incidence_matrices), function(k) incidence_matrices[[k]] /nr_sample))
      weighted_betas[[i+1]] <- weighted_betas_proposed
      is_weights [[i+1]] <- weights_proposed
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
      is_weights [[i+1]] <- is_weights [[i]]
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
              betas = weighted_betas[-c(1:burin_iter)],
              diffBiDAGs = diff_BiDAGs[-c(1:burin_iter)]
              # ,totalscore_of_DAGs = totalscore_of_DAGs[-c(1:burin_iter)]
              )
         )
}


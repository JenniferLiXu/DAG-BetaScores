
SNOrderSampler_OneDAG <- function(n, iteration, order_iter, order = NULL, 
                                  order_stepsize, moveprobs, base_score = 0, 
                                  starting_dag = NULL, betas_init = NULL, skeleton = FALSE,
                                  edgesposterior, burnin = 0.2) {
  
  elements <- c(1:n) # Define set of elements
  # permutations <- permn(elements) # Generate all permutations
  nr_sample <- 1 # Number of generated DAGs under one order
  
  # Initialize order
  if (is.null(order)) {
    order <-  list(seq_len(n))
  }
  # Initialize starting DAG if not provided
  if (is.null(starting_dag)) {
    # starting_dag <- list(matrix(0, nrow = n, ncol = n))
    starting_dag <-list(samplescore(n, betas = matrix(0, nrow = n, ncol = n), order[[1]])$incidence)
  }
  # Initialize beta matrix
  if (is.null(betas_init)) {
    calcultion_betas_init <- calculateBetaScoresArray_hash(starting_dag, k = 1, n, base_score = 0)
    betas_init <- calcultion_betas_init$allBetaScores[,,1]
    base_score <- 0
    # Idea: subtract the highest order score
    # prev_logscore_list <- sapply(1:length(permutations), function(k) sum(orderscore_betas(n, c(1:n), betas_init, permutations[[k]])))
    # totalscore_orders_prev <- calculate_final_score(prev_logscore_list, "sum")
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
  
  init_sampled_DAGs <- lapply(1:nr_sample, function(x) samplescore(n, betas_init, order[[1]]))
  DAG_prev <- lapply(init_sampled_DAGs, function(dag) dag$incidence) 
  BiDAGscore_prev <- calculateBetaScoresArray_hash(DAG_prev, k = length(DAG_prev), n, base_score)$target_DAG_score
  # trace_totalscore <- numeric()
  
  # Estimation for normal constant
  order_set_prev <- order
  DAG_set_prev <- DAG_prev
  logscore_set_prev <- lapply(init_sampled_DAGs, function(dag) dag$logscore) 
  oDAGoBeta_set_logscore <- lapply(1:length(DAG_set_prev), function(k) calculate_DAG_score(DAG_list = DAG_set_prev[k], permy = order_set_prev[[k]], weights = NULL ,
                                                                                           betas = betas_init, base_score = base_score))
  
  # Looping through iterations
  for (i in 1:iter) {
    beta_prev <- weighted_betas[[i]]
    order_prev <- order[[i]]
    
    # Sampling orders with OrderMCMC
    example <- orderMCMC_betas(n,startorder = order_prev ,iterations = order_iter, betas = beta_prev,
                               stepsave = order_stepsize, moveprobs)  
    
    permy <- unlist(example[[4]][length(example[[4]])])#Store the last order from the chain
    
    #  Sample 1 DAG using the new order from OrderMCMC under previous beta matrix
    sampled_DAGs <- lapply(1:nr_sample, function(x) samplescore(n, beta_prev, permy))
    
    # Extracting incidence matrices and log scores (old beta matrix)
    incidence_matrices <-lapply(sampled_DAGs, function(dag) dag$incidence)
    incidence_logscore<-lapply(sampled_DAGs, function(dag) dag$logscore) # scores of new DAGs under old beta
    
    # Ordres used for the set of DAGs
    proposed_orders <- example[[4]][-1]
    
    # Update beta matrix using the weights from sampled DAGs
    calculation_beta_values <- calculateBetaScoresArray_hash(incidence_matrices, k = length(incidence_matrices) ,n, base_score = base_score)
    BiDAGscore_propose <- calculation_beta_values$target_DAG_score
    beta_values <- calculation_beta_values$allBetaScores
    
    # Update beta matrix using importance sampling
    is_results <- importance_DAG(DAGs = incidence_matrices, score_under_betas = incidence_logscore, target_scores = BiDAGscore_propose)
    weights_proposed <- is_results$importance_weights # normalized weights under old beta
    
    #New beta matrix using the normalised weights
    weighted_betas_proposed <- Reduce("+", lapply(1:length(weights_proposed), 
                                                  function(k) beta_values[,,k] * weights_proposed[k]))
    
    # proposed_logscore_list <- sapply(1:length(permutations),
    #                                  function(k) sum(orderscore_betas(n, c(1:n), weighted_betas_proposed, permutations[[k]])))
    # proposed_totalscore_orders <-calculate_final_score(proposed_logscore_list, "sum")
    
    # orderscore_prev <- unlist(example[[3]][length(example[[3]])]) #new order, previous beta
    # orderscore_prop <- sum(orderscore_betas(n,c(1:n), weighted_betas_proposed, order_prev))
    prev_order_logscore_list <- sapply(1:length(order_set_prev), function(k) sum(orderscore_betas(n, c(1:n), beta_prev, order_set_prev[[k]])))
    prop_order_logscore_list <- sapply(1:length(order_set_prev), function(k) sum(orderscore_betas(n, c(1:n), weighted_betas_proposed, order_set_prev[[k]])))
    orderscore_prev <- unlist(example[[3]][-1]) # score of new order, previous beta
    orderscore_prop <- sapply(1:length(proposed_orders), function(k) sum(orderscore_betas(n, c(1:n), weighted_betas_proposed, proposed_orders[[k]])))
    
    set_A_order <- calculate_final_score(c(unlist(prev_order_logscore_list), unlist(orderscore_prev)), "sum")
    set_B_order <- calculate_final_score(c(unlist(prop_order_logscore_list), unlist(orderscore_prop)), "sum")
    alternative_4 <- set_A_order-set_B_order
    # Log score of old DAG set under the new beta
    oDAGnBeta_logscore <- calculate_DAG_score(DAG_list = DAG_prev ,permy = order_prev, weights = NULL ,
                                              betas = weighted_betas_proposed, base_score = base_score)
    # Acceptance ratio
    # Test
    wB <- is_results$log_diff
    inv_wA <- oDAGnBeta_logscore - BiDAGscore_prev
    # print(wB + inv_wA)
    # trace_totalscore[i] <- totalscore_orders_prev - proposed_totalscore_orders
    
    log_ratio <- wB + inv_wA  + alternative_4
    ratio <- exp(log_ratio)
    # cat("alternative_4",alternative_4,"trace_totalscore",trace_totalscore[i], "\n")   
    # Metropolis-Hastings acceptance step
    if(runif(1) < ratio){ 
      compress_DAG[[i+1]] <- is_results$compress_dag
      weighted_betas[[i+1]] <- weighted_betas_proposed
      order[[i+1]] <- permy
      count_accept[i] <- 1 # Accept 
      BiDAGscore_prev <- BiDAGscore_propose
      DAG_prev <- incidence_matrices
      # totalscore_orders_prev <- proposed_totalscore_orders
      order_set_prev <- proposed_orders
      # cat("wB",wB, "inv_wA", inv_wA, "\n")    
      # print(compress_DAG[[i+1]])
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
    }
    
    edge_over_time[,,i] <- current_mat
    
  }
  
  # Return the results
  return(list(DAGs = compress_DAG[-c(1:burin_iter)], 
              edge_prob = edge_over_time[,,-c(1:burin_iter)], 
              essValues = ess_DAGs[-c(1:burin_iter)], 
              acceptCount = count_accept[-c(1:burin_iter)], 
              betas = weighted_betas[-c(1:burin_iter)]
              # trace_totalscore = trace_totalscore
              # ,totalscore_of_DAGs = totalscore_of_DAGs[-c(1:burin_iter)]
  )
  )
}

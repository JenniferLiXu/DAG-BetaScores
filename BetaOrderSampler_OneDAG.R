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
# order_stepsize = 20
# moveprobs = moveprobs
# edgesposterior = edgesposterior

BetaOrderSampler_OneDAG <- function(n, iteration, order_iter, order = NULL, 
                                    order_stepsize, moveprobs, base_score = 0, 
                                    starting_dag = NULL, betas_init = NULL, skeleton = FALSE,
                                    edgesposterior, burnin = 0.4) {
  
  elements <- c(1:n) # Define set of elements
  permutations <- permn(elements) # Generate all permutations
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
    # base_score <- calcultion_betas_init$target_DAG_score
    base_score <- 0
    # Idea: subtract the highest order score
    prev_logscore_list <- sapply(1:length(permutations), function(k) sum(orderscore_betas(n, c(1:n), betas_init, permutations[[k]])))
    totalscore_orders_prev <- calculate_final_score(prev_logscore_list, "sum")
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
  BiDAGscore_prev <- calculateBetaScoresArray_hash(DAG_prev, k = length(DAG_prev), n, base_score)$target_DAG_score
  
  trace_totalscore <- numeric()
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
    
    # incidence_matrices <- example[[1]][length(example[[1]])]
    # incidence_logscore <- example[[2]][length(example[[2]])]
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
    
    proposed_logscore_list <- sapply(1:length(permutations),
                                     function(k) sum(orderscore_betas(n, c(1:n), weighted_betas_proposed, permutations[[k]])))
    proposed_totalscore_orders <-calculate_final_score(proposed_logscore_list, "sum")
    # cat("time for order MCMC:",endtime_order, ",DAGs",endtime_DAGs, ",beta",endtime_beta, "\n")
    
    # orderscore_prev <- unlist(example[[3]][length(example[[3]])]) #new order, previous beta
    # orderscore_prop <- sum(orderscore_betas(n,c(1:n), weighted_betas_proposed, order_prev))
    
    # Log score of old DAG set under the new beta
    oDAGnBeta_logscore <- calculate_DAG_score(DAG_list = DAG_prev ,permy = order_prev, weights = NULL ,
                                              betas = weighted_betas_proposed, base_score = base_score)
    # Acceptance ratio
    # Test
    wB <- is_results$log_diff
    inv_wA <- oDAGnBeta_logscore - BiDAGscore_prev
    # print(wB + inv_wA)
    log_ratio <- wB + inv_wA  + totalscore_orders_prev - proposed_totalscore_orders
    ratio <- exp(log_ratio)
    trace_totalscore[i] <- totalscore_orders_prev - proposed_totalscore_orders
    # cat("trace_totalscore",trace_totalscore[i], "\n")   
    # Metropolis-Hastings acceptance step
    if(runif(1) < ratio){ 
      compress_DAG[[i+1]] <- is_results$compress_dag
      weighted_betas[[i+1]] <- weighted_betas_proposed
      # print(weighted_betas[[i+1]])
      # print(compress_DAG[[i+1]])
      order[[i+1]] <- permy
      count_accept[i] <- 1 # Accept 
      BiDAGscore_prev <- BiDAGscore_propose
      DAG_prev <- incidence_matrices
      totalscore_orders_prev <- proposed_totalscore_orders
      # cat("wB",wB, "inv_wA", inv_wA, "\n")    
      # print(permy)
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
  return(list(DAGs = compress_DAG[-c(1:burin_iter)], 
              edgeDifferences = edge_diff_over_time[,,-c(1:burin_iter)], 
              edge_prob = edge_over_time[,,-c(1:burin_iter)], 
              essValues = ess_DAGs[-c(1:burin_iter)], 
              acceptCount = count_accept[-c(1:burin_iter)], 
              betas = weighted_betas[-c(1:burin_iter)],
              diffBiDAGs = diff_BiDAGs[-c(1:burin_iter)],
              trace_totalscore = trace_totalscore
              # ,totalscore_of_DAGs = totalscore_of_DAGs[-c(1:burin_iter)]
  )
  )
}


BetaOrderSampler_OneDAG_test <- function(n, iteration, order_iter, order = NULL, 
                                    order_stepsize, moveprobs, base_score = 0, 
                                    starting_dag = NULL, betas_init = NULL, skeleton = FALSE,
                                    edgesposterior, burnin = 0.2) {
  
  elements <- c(1:n) # Define set of elements
  permutations <- permn(elements) # Generate all permutations
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
    # base_score <- calcultion_betas_init$target_DAG_score
    base_score <- 0
    # Idea: subtract the highest order score
    prev_logscore_list <- sapply(1:length(permutations), function(k) sum(orderscore_betas(n, c(1:n), betas_init, permutations[[k]])))
    totalscore_orders_prev <- calculate_final_score(prev_logscore_list, "sum")
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
  BiDAGscore_prev <- calculateBetaScoresArray_hash(DAG_prev, k = length(DAG_prev), n, base_score)$target_DAG_score
  
  # Estimation for normal constant
  order_set_prev <- order
  DAG_set_prev <- DAG_prev
  logscore_set_prev <- lapply(init_sampled_DAGs, function(dag) dag$logscore) 
  oDAGoBeta_set_logscore <- lapply(1:length(DAG_set_prev), function(k) calculate_DAG_score(DAG_list = DAG_set_prev[k], permy = order_set_prev[[k]], weights = NULL ,
                                                                                           betas = betas_init, base_score = base_score))
  trace_totalscore <- numeric()
  # Looping through iterations
  for (i in 1:iter) {
    beta_prev <- weighted_betas[[i]]
    order_prev <- order[[i]]
    
    # Sampling orders with OrderMCMC
    example <- orderMCMC_betas(n,startorder = order_prev ,iterations = order_iter, betas = beta_prev,
                               stepsave = order_stepsize, moveprobs)  
    
    permy <- unlist(example[[4]][length(example[[4]])])#Store the last order from the chain
    
    # #  Sample 1 DAG using the new order from OrderMCMC under previous beta matrix
    sampled_DAGs <- lapply(1:nr_sample, function(x) samplescore(n, beta_prev, permy))
    # Extracting incidence matrices and log scores (old beta matrix)
    incidence_matrices <-lapply(sampled_DAGs, function(dag) dag$incidence)
    incidence_logscore<-lapply(sampled_DAGs, function(dag) dag$logscore) # scores of new DAGs under old beta
    
    # Ordres used for the set of DAGs
    proposed_orders <- example[[4]][-1]
    # sampled_DAGset_fromOrder <- DAGs_from_order(order_list = proposed_orders, nr_sample = 1, beta_matrix = beta_prev)
    # DAG_set_prop <- sampled_DAGset_fromOrder$incidence # List of DAGs sampled under previous beta
    # logscore_set_prop <- sampled_DAGset_fromOrder$logscore # List of logscores of new sampled DAGs under previous beta, nDAGoBeta_set_logscore
    
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
    
    proposed_logscore_list <- sapply(1:length(permutations), function(k) sum(orderscore_betas(n, c(1:n), weighted_betas_proposed, permutations[[k]])))
    proposed_totalscore_orders <-calculate_final_score(proposed_logscore_list, "sum")

    # Log score of old DAG set under the new beta
    oDAGnBeta_logscore <- calculate_DAG_score(DAG_list = DAG_prev ,permy = order_prev, weights = NULL ,
                                              betas = weighted_betas_proposed, base_score = base_score)
    
    # oDAGnBeta_set_logscore <- lapply(1:length(DAG_set_prev), function(k) calculate_DAG_score(DAG_list = DAG_set_prev[k], permy = order_set_prev[[k]], weights = NULL ,
    #                                                                                          betas = weighted_betas_proposed, base_score = base_score))
    # 
    # nDAGnBeta_set_logscore <- lapply(1:length(DAG_set_prop), function(k) calculate_DAG_score(DAG_list = DAG_set_prop[k], permy = sampled_DAGset_fromOrder$order[[k]], weights = NULL ,
    #                                                                                          betas = weighted_betas_proposed, base_score = base_score))
    
    prev_order_logscore_list <- sapply(1:length(order_set_prev), function(k) sum(orderscore_betas(n, c(1:n), beta_prev, order_set_prev[[k]])))
    prev_orderscore_sum <- calculate_final_score(prev_order_logscore_list, "sum")

    prop_order_logscore_list <- sapply(1:length(order_set_prev), function(k) sum(orderscore_betas(n, c(1:n), weighted_betas_proposed, order_set_prev[[k]])))
    prop_orderscore_sum <- calculate_final_score(prop_order_logscore_list, "sum")

    orderscore_prev <- unlist(example[[3]][-1]) # score of new order, previous beta
    prev_new_orderscore_sum <- calculate_final_score(orderscore_prev, "sum")
    
    orderscore_prop <- sapply(1:length(proposed_orders), function(k) sum(orderscore_betas(n, c(1:n), weighted_betas_proposed, proposed_orders[[k]])))
    prop_new_orderscore_sum <- calculate_final_score(orderscore_prop, "sum")
    # cat("prev_orderscore_sum",prev_orderscore_sum, "prev_new_orderscore_sum",prev_new_orderscore_sum, "\n") 
    # cat("prev_new_orderscore_sum",prev_new_orderscore_sum, "prop_new_orderscore_sum",prop_new_orderscore_sum, "\n")
    
    alternative_3 <- prev_orderscore_sum - prop_orderscore_sum + prev_new_orderscore_sum - prop_new_orderscore_sum
    
    set_A_order <- calculate_final_score(c(unlist(prev_order_logscore_list), unlist(orderscore_prev)), "sum")
    set_B_order <- calculate_final_score(c(unlist(prop_order_logscore_list), unlist(orderscore_prop)), "sum")
    alternative_4 <- set_A_order-set_B_order
    
    # set_B <- calculate_final_score(unlist(oDAGnBeta_set_logscore), "sum")
    # set_A <- calculate_final_score(unlist(oDAGoBeta_set_logscore), "sum")
    # approx_order_score <- set_A-set_B
    
    # set_B_ON <- calculate_final_score(c(unlist(oDAGnBeta_set_logscore), unlist(nDAGnBeta_set_logscore)), "sum")
    # set_A_ON <- calculate_final_score(c(unlist(logscore_set_prop),unlist(oDAGoBeta_set_logscore)), "sum")
    # alternative <- set_A_ON - set_B_ON
    
    # alternative_2 <- prev_orderscore_sum - prop_orderscore_sum
    # if(i == 1){alternative_2 = 0}
    # Acceptance ratio
    # Test
    wB <- is_results$log_diff
    inv_wA <- oDAGnBeta_logscore - BiDAGscore_prev
    trace_totalscore[i] <- totalscore_orders_prev - proposed_totalscore_orders
    log_ratio <- wB + inv_wA  +  trace_totalscore[i]
    # cat("approx_order_score",approx_order_score, "alternative",alternative, "\n")
    cat("alternative_3",alternative_3,  "alternative_4",alternative_4,"\n")
    cat("trace_totalscore",trace_totalscore[i], "\n")
    ratio <- exp(log_ratio)

    # Metropolis-Hastings acceptance step
    if(runif(1) < ratio){ 
      compress_DAG[[i+1]] <- is_results$compress_dag
      weighted_betas[[i+1]] <- weighted_betas_proposed
      order[[i+1]] <- permy
      count_accept[i] <- 1 # Accept 
      BiDAGscore_prev <- BiDAGscore_propose
      totalscore_orders_prev <- proposed_totalscore_orders
      order_set_prev <- proposed_orders
      print(is_results$compress_dag)
      # DAG_set_prev <- DAG_set_prop
      # oDAGoBeta_set_logscore <- nDAGnBeta_set_logscore
      # cat("wB",wB, "inv_wA", inv_wA, "\n")    
      # cat("approx_order_score",approx_order_score, "alternative",alternative, "alternative_2",alternative_2, "\n") 
      # cat("trace_totalscore",trace_totalscore[i], "\n")   
    }else{
      compress_DAG[[i+1]] <- compress_DAG[[i]]
      weighted_betas[[i+1]] <- weighted_betas[[i]]
      order[[i+1]] <- order[[i]]
      count_accept[i] <- 0 # Reject
    }
    
    if (length(compress_DAG)-1 > burin_iter) {
      total_DAG <- total_DAG + compress_DAG[[i+1]]
      current_mat <- total_DAG/(i - burin_iter) # Average the edges of DAGs after burn in part
      # diff_mat <- CompareDAG(current_mat, edgesposterior)
    }else{
      sum_matrix <- Reduce("+", compress_DAG[1:length(compress_DAG)])
      current_mat <- sum_matrix/i
      # diff_mat <- CompareDAG(current_mat, edgesposterior)
    }
    
    edge_over_time[,,i] <- current_mat
    # edge_diff_over_time[,,i] <- diff_mat # Store the difference per edge
    # diff_BiDAG <-norm(diff_mat,type = "F") # Store the difference of the DAG matrix
    # diff_BiDAGs <- c(diff_BiDAGs, diff_BiDAG)
  }
  
  # Return the results
  return(list(DAGs = compress_DAG[-c(1:burin_iter)], 
              edgeDifferences = edge_diff_over_time[,,-c(1:burin_iter)], 
              edge_prob = edge_over_time[,,-c(1:burin_iter)], 
              essValues = ess_DAGs[-c(1:burin_iter)], 
              acceptCount = count_accept[-c(1:burin_iter)], 
              betas = weighted_betas[-c(1:burin_iter)]
              # diffBiDAGs = diff_BiDAGs[-c(1:burin_iter)]
              # trace_totalscore = trace_totalscore
              # ,totalscore_of_DAGs = totalscore_of_DAGs[-c(1:burin_iter)]
  )
  )
}


###################### DAG sets
# 10 DAGs in set 
# iteration 3000: 0.000659785 0.9957892 
# iteration 5000: 0.0006612062 0.9958414 

# new 10 DAGs in set 
# iteration 3000: 0.001177898 0.9821196
# iteration 5000: 0.0005802811 0.9942173  

# new 100 DAGs
# iteration 2000: 0.00179923 0.9724962 
# iteration 5000:  0.0009047404 0.9874083 
# iteration 10000: 0.0008012927 0.991951 

# 100 DAGs in set
# iteration 3000: 0.004165124 0.9493413 
# iteration 5000: 0.003184019 0.9667746 


# true beta: results_seed123$betas
# [,1]      [,2]      [,3]      [,4]
# [1,]  0.000000 -6.290409 42.525494  1.477611
# [2,]  3.320621  0.000000  6.585970 -4.721817
# [3,] 42.525494 16.197000  0.000000 -3.730036
# [4,]  1.477611 -5.038134 -3.730036  0.000000

###################### Gibbs sampling
num_iterations <- 2000
# weighted version         order_betaOrder123:  0.00448463 0.9778996 
# equally weighted version order_betaOrder300:  0.004134801 0.9740147 

num_iterations <- 3000
# weighted version         order_betaOrder123:  0.007491891 0.9452089 
# equally weighted version order_betaOrder300:  0.003967594 0.9751654 

num_iterations <- 5000
# weighted version         order_betaOrder123:  0.006823928 0.9670896 
# equally weighted version order_betaOrder300:  0.003714876 0.9757753

## With ess > 0.3 (10 sampled orders, each with 10 DAGs)
# iteration 3000:  0.002242364 0.9668648  burnin:0.5  acceptCount_after_burnin:131
# iteration 4000:  0.001860352 0.9717618   burnin:0.5  acceptCount_after_burnin:290

# With ess > 0.4 (10 sampled orders, each with 10 DAGs) ---- modified beta
num_iterations <- 3000
# Method 2 with MSE: 0.005685223  R2: 0.9355907 

num_iterations <- 5000
# Method 2 with MSE: 0.006130747 R2: 0.9439851 

num_iterations <- 8000
# Method 2 with MSE:0.00535461 R2:0.9424123 

###################### End of Gibbs sampling
# > pedges_comp
# [[1]]
# crim          zn      indus        chas
# crim  0.000000000 0.001248439 0.38202247 0.538077403
# zn    0.001248439 0.000000000 0.19475655 0.007490637
# indus 0.617977528 0.805243446 0.00000000 0.063670412
# chas  0.136079900 0.001248439 0.01248439 0.000000000

wB <- sum(is_results$log_diff)/nr_sample
inv_wA <- sum((oDAGnBeta_logscore_list - BiDAGscore_prev_list))/nr_sample
compress_DAG[[i+1]] <- is_results$compress_dag
# [[2]]
# crim          zn      indus        chas
# crim  0.0000000000 0.001320000 0.34600000 0.521212726
# zn    0.0000822884 0.000000000 0.41920000 0.008688032
# indus 0.6540000000 0.580800000 0.00000000 0.042286833
# chas  0.1346322698 0.004206033 0.02186849 0.000000000

compress_DAG[[i+1]] <- Reduce("+", lapply(1:length(incidence_matrices), function(k) incidence_matrices[[k]]/nr_sample))
# [[2]]
# crim      zn   indus    chas
# crim  0.00000 0.00132 0.34600 0.51978
# zn    0.00008 0.00000 0.41920 0.00990
# indus 0.65400 0.58080 0.00000 0.05522
# chas  0.13516 0.00412 0.03588 0.00000

wB <- sum(is_results$log_diff*weights_proposed)
inv_wA <- sum((oDAGnBeta_logscore_list - BiDAGscore_prev_list)*weights_prev)
compress_DAG[[i+1]] <- is_results$compress_dag
# [[2]]
# crim          zn      indus        chas
# crim  0.000000000 0.001879281 0.23280000 0.489124754
# zn    0.001678616 0.000000000 0.05880000 0.008113045
# indus 0.767200000 0.941200000 0.00000000 0.075329203
# chas  0.178180226 0.002812892 0.02384516 0.000000000
results_seed100$betas
# [,1]      [,2]      [,3]       [,4]
# [1,]  0.0000000 -6.290409 42.525494  0.5029098
# [2,] -6.3312354  0.000000 16.197000 -4.0258435
# [3,] 42.0381171 16.197000  0.000000 -1.2932571
# [4,] -0.9591683 -5.038134 -1.293257  0.0000000


wB <- sum(is_results$log_diff*weights_proposed)
inv_wA <- sum((oDAGnBeta_logscore_list - BiDAGscore_prev_list)*weights_prev)
compress_DAG[[i+1]] <- Reduce("+", lapply(1:length(incidence_matrices), function(k) incidence_matrices[[k]]/nr_sample))
# [[2]]
# crim      zn   indus    chas
# crim  0.00000 0.04470 0.23280 0.49600
# zn    0.12220 0.00000 0.05880 0.00930
# indus 0.76720 0.94120 0.00000 0.07516
# chas  0.16414 0.00414 0.04244 0.00000

wB <- sum(is_results$log_diff*weights_proposed)
inv_wA <- sum((oDAGnBeta_logscore_list - BiDAGscore_prev_list)*weights_for_prevD)
compress_DAG[[i+1]] <- is_results$compress_dag
# [[2]]
# crim          zn      indus       chas
# crim  0.000000000 0.001816183 0.18580000 0.47896311
# zn    0.002069068 0.000000000 0.05180000 0.00748412
# indus 0.814200000 0.948200000 0.00000000 0.07932004
# chas  0.176958202 0.006009231 0.02318009 0.00000000
results_seed100$betas
# [,1]      [,2]      [,3]       [,4]
# [1,]  0.0000000 -6.290409 42.525494  0.5029098
# [2,] -6.3312354  0.000000 16.197000 -4.0258435
# [3,] 42.0381171 16.197000  0.000000 -1.2932571
# [4,] -0.9591683 -5.038134 -1.293257  0.0000000


results_seed123$betas
# [,1]      [,2]      [,3]      [,4]
# [1,]  0.0000000  3.320621 32.914464  1.477611
# [2,] -6.4945315  0.000000 16.197000 -3.350996
# [3,] 40.0887149 16.197000  0.000000 -1.293257
# [4,] -0.9591683 -3.350996 -2.980395  0.000000



BetaOrderSampler_OneDAG_test<- function(n, iteration, order_iter, order = NULL, 
                                        order_stepsize, moveprobs, base_score = 0, 
                                        starting_dag = NULL, betas_init = NULL, skeleton = FALSE,
                                        edgesposterior, burnin = 0.2) {
  
  elements <- c(1:n) # Define set of elements
  nr_sample <- 1 # Number of generated DAGs under one order
  # permutations <- permn(elements) # Generate all permutations
  
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
  
  # Sampled DAG
  init_sampled_DAGs <- lapply(1:nr_sample, function(x) samplescore(n, betas_init, order[[1]]))
  DAG_prev <- lapply(init_sampled_DAGs, function(dag) dag$incidence) 
  
  # Estimation for normal constant
  order_set_prev <- order
  DAG_set_prev <- DAG_prev
  logscore_set_prev <- lapply(init_sampled_DAGs, function(dag) dag$logscore) 
  oDAGoBeta_set_logscore <- lapply(1:length(DAG_set_prev), function(k) calculate_DAG_score(DAG_list = DAG_set_prev[k], permy = order_set_prev[[k]], weights = NULL ,
                                                                                           betas = betas_init, base_score = base_score))
  
  BiDAGscore_prev <- calculateBetaScoresArray_hash(DAG_prev, k = length(DAG_prev), n, base_score)$target_DAG_score
  
  trace_totalscore <- numeric()
  # Looping through iterations
  for (i in 1:iter) {
    beta_prev <- weighted_betas[[i]]
    order_prev <- order[[i]]
    # Sampling orders with OrderMCMC
    example <- orderMCMC_betas(n,startorder = order_prev ,iterations = order_iter, betas = beta_prev,
                               stepsave = order_stepsize, moveprobs)  
    
    permy <- unlist(example[[4]][length(example[[4]])])# Store the last order from the chain
    #  Sample 1 DAG using the new order from OrderMCMC under previous beta matrix
    sampled_DAGs <- lapply(1:nr_sample, function(x) samplescore(n, beta_prev, permy))
    # Extracting incidence matrices and log scores (old beta matrix)
    incidence_matrices <-lapply(sampled_DAGs, function(dag) dag$incidence) 
    incidence_logscore<-lapply(sampled_DAGs, function(dag) dag$logscore) # scores of new DAGs under old beta
    
    # Ordres used for the set of DAGs
    proposed_orders <- example[[4]][-1]
    sampled_DAGset_fromOrder <- DAGs_from_order(order_list = proposed_orders, nr_sample = 5, beta_matrix = beta_prev)
    DAG_set_prop <- sampled_DAGset_fromOrder$incidence # List of DAGs sampled under previous beta
    logscore_set_prop <- sampled_DAGset_fromOrder$logscore # List of logscores of new sampled DAGs under previous beta, nDAGoBeta_set_logscore
    
    # Update beta matrix using the weights from sampled DAGs
    calculation_beta_values <- calculateBetaScoresArray_hash(incidence_matrices, k = length(incidence_matrices) ,n, base_score = base_score)
    BiDAGscore_propose <- calculation_beta_values$target_DAG_score
    beta_values <- calculation_beta_values$allBetaScores
    
    # Update beta matrix using importance sampling
    is_results <- importance_DAG(DAGs = incidence_matrices, score_under_betas = incidence_logscore, target_scores = BiDAGscore_propose)
    weights_proposed <- is_results$importance_weights # normalized weights under old beta
    #New beta matrix using the normalised weights
    weighted_betas_proposed <- Reduce("+", lapply(1:length(weights_proposed), function(k) beta_values[,,k] * weights_proposed[k]))
    
    proposed_logscore_list <- sapply(1:length(permutations),
                                     function(k) sum(orderscore_betas(n, c(1:n), weighted_betas_proposed, permutations[[k]])))
    proposed_totalscore_orders <-calculate_final_score(proposed_logscore_list, "sum")
    
    # orderscore_prev <- unlist(example[[3]][length(example[[3]])]) #new order, previous beta
    # orderscore_prop <- sum(orderscore_betas(n,c(1:n), weighted_betas_proposed, order_prev))
    
    # Log score of old DAG set under the new beta
    oDAGnBeta_logscore <- calculate_DAG_score(DAG_list = DAG_prev ,permy = order_prev, weights = NULL ,
                                              betas = weighted_betas_proposed, base_score = base_score)
    
    oDAGnBeta_set_logscore <- lapply(1:length(DAG_set_prev), function(k) calculate_DAG_score(DAG_list = DAG_set_prev[k], permy = order_set_prev[[k]], weights = NULL ,
                                                                                             betas = weighted_betas_proposed, base_score = base_score))
    
    nDAGnBeta_set_logscore <- lapply(1:length(DAG_set_prop), function(k) calculate_DAG_score(DAG_list = DAG_set_prop[k], permy = sampled_DAGset_fromOrder$order[[k]], weights = NULL ,
                                                                                             betas = weighted_betas_proposed, base_score = base_score))
    
    
    prev_order_logscore_list <- sapply(1:length(order_set_prev), function(k) sum(orderscore_betas(n, c(1:n), beta_prev, order_set_prev[[k]])))
    prev_orderscore_sum <- calculate_final_score(prev_order_logscore_list, "sum")
    
    prop_order_logscore_list <- sapply(1:length(order_set_prev), function(k) sum(orderscore_betas(n, c(1:n), weighted_betas_proposed, order_set_prev[[k]])))
    prop_orderscore_sum <- calculate_final_score(prop_order_logscore_list, "sum")
    
    # Acceptance ratio
    # Test
    wB <- is_results$log_diff
    inv_wA <- oDAGnBeta_logscore - BiDAGscore_prev
    
    set_B <- calculate_final_score(unlist(oDAGnBeta_set_logscore), "sum")
    set_A <- calculate_final_score(unlist(oDAGoBeta_set_logscore), "sum")
    # approx_order_score <- calculate_final_score(unlist(oDAGnBeta_set_logscore) - unlist(logscore_set_prev), "sum")
    approx_order_score <- set_A-set_B
    
    set_B_ON <- calculate_final_score(c(unlist(oDAGnBeta_set_logscore), unlist(nDAGnBeta_set_logscore)), "sum")
    set_A_ON <- calculate_final_score(c(unlist(logscore_set_prop),unlist(oDAGoBeta_set_logscore)), "sum")
    
    alternative <- set_A_ON - set_B_ON
    
    alternative_2 <- prev_orderscore_sum - prop_orderscore_sum
    cat("alternative",alternative, "alternative_2",alternative_2, "\n") 
    
    trace_totalscore[i] <- totalscore_orders_prev - proposed_totalscore_orders
    cat("trace_totalscore",trace_totalscore[i], "\n")    
    
    log_ratio <- wB + inv_wA + trace_totalscore[i]
    # print(log_ratio)
    # cat("approx_order_score",approx_order_score, "\n") 
    # print(wB + inv_wA  + totalscore_orders_prev - proposed_totalscore_orders)
    # cat("trace_totalscore",totalscore_orders_prev - proposed_totalscore_orders, "\n")   
    # log_ratio <- wB + inv_wA  + totalscore_orders_prev - proposed_totalscore_orders
    ratio <- exp(log_ratio)
    
    
    # Metropolis-Hastings acceptance step
    if(runif(1) < ratio){ 
      compress_DAG[[i+1]] <- is_results$compress_dag
      weighted_betas[[i+1]] <- weighted_betas_proposed
      order[[i+1]] <- permy
      count_accept[i] <- 1 # Accept 
      BiDAGscore_prev <- BiDAGscore_propose
      DAG_prev <- incidence_matrices
      # totalscore_orders_prev <- proposed_totalscore_orders
      order_set_prev <- sampled_DAGset_fromOrder$order
      DAG_set_prev <- DAG_set_prop
      oDAGoBeta_set_logscore <- nDAGnBeta_set_logscore
      print("accepted")
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
              betas = weighted_betas[[iter+1]],
              diffBiDAGs = diff_BiDAGs[-c(1:burin_iter)],
              trace_totalscore = trace_totalscore
              # ,totalscore_of_DAGs = totalscore_of_DAGs[-c(1:burin_iter)]
  )
  )
}


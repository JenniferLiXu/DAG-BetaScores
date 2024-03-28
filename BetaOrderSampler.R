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

# time for order MCMC: 0.009 0 0.009 
# time for Sample 30 DAGs: 0.007 0 0.008 
# time for Beta matrix: 0.007 0 0.007 

BetaOrderSampler <- function(n, iteration, order_iter, order = NULL, 
                             order_stepsize, moveprobs, base_score = 0, 
                             starting_dag = NULL, betas_init = NULL, skeleton = FALSE,
                             edgesposterior, burnin = 0.2 ) {
  nr_sample <- 30
  # Initialize starting DAG if not provided
  if (is.null(starting_dag)) {
    starting_dag <- list(matrix(0, nrow = n, ncol = n))
  }
  # Initialize beta matrix
  if (is.null(betas_init)) {
    calcultion_betas_init <- calculateBetaScoresArray_hash(starting_dag, k = 1, n)
    betas_init <- calcultion_betas_init$allBetaScores[,,1]
    totalscore_of_DAGs <- numeric()
    totalscore_of_DAGs[1] <- nr_sample*calcultion_betas_init$target_DAG_score
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
  burin_iter <- floor(burnin*iteration)
  iter <- iteration + burin_iter
  
  # DAG <- starting_dag
  init_sampled_DAGs <- lapply(1:nr_sample, function(x) samplescore(n, betas_init, order[[1]], base_score))
  prev_DAG <- lapply(init_sampled_DAGs, function(dag) dag$incidence) 
  totalscore_of_DAGs_prev <- calculateBetaScoresArray_hash(prev_DAG, k = length(prev_DAG), n)$target_DAG_score
  compress_DAG <- starting_dag
  total_DAG <- matrix(0, nrow = n, ncol = n)
  edge_over_time <- array(0, dim = c(n, n, iter))
  edge_diff_over_time <- array(0, dim = c(n, n, iter))
  prev_weight <- 1
  
  
  # Looping through iterations
  for (i in 1:iter) {
    beta_prev <- weighted_betas[[i]]
    order_prev <- order[[i]]
    
    starttime_order<-proc.time() # for timing the problem
    # Sampling orders with OrderMCMC
    example <- orderMCMC_betas(n,startorder = order_prev ,iterations = order_iter, betas = beta_prev,
                               stepsave = order_stepsize, moveprobs) # run the Order MCMC code
    
    endtime_order<-proc.time()
    endtime_order<-endtime_order-starttime_order
    
    #Store the last order from the chain
    permy <- unlist(example[[4]][length(example[[4]])])
    
    starttime_DAGs<-proc.time() # for timing the problem
    #  Sample nr_sample DAGs using the last sampled order from OrderMCMC and old beta matrix
    sampled_DAGs <- lapply(1:nr_sample, function(x) samplescore(n, beta_prev, permy, base_score))
    
    endtime_DAGs<-proc.time()
    endtime_DAGs<-endtime_DAGs-starttime_DAGs
    
    # Extracting incidence matrices and log scores (old beta matrix)
    incidence_matrices <-lapply(sampled_DAGs, function(dag) dag$incidence) 
    incidence_logscore<-lapply(sampled_DAGs, function(dag) dag$logscore) # scores of new DAGs under old beta
    
    starttime_beta<-proc.time() # for timing the problem
    # Update beta matrix using the weights from sampled DAGs
    calculation_beta_values <- calculateBetaScoresArray_hash(incidence_matrices, k = length(incidence_matrices) ,n)
    totalscore_DAGs <- sum(calculation_beta_values$target_DAG_score)
    beta_values <- calculation_beta_values$allBetaScores
    
    endtime_beta<-proc.time()
    endtime_beta<-endtime_beta-starttime_beta
    
    # Update beta matrix using importance sampling
    is_results <- importance_DAG(DAGs = incidence_matrices, score_under_betas = incidence_logscore, target_scores = calculation_beta_values$target_DAG_score)
    weights_proposed <- is_results$importance_weights # normalised weights under old beta
    
    #New beta matrix using the normalised weights
    weighted_betas_proposed <- Reduce("+", lapply(1:length(weights_proposed), 
                                                  function(k) beta_values[,,k] * weights_proposed[k]))
    
    # cat("time for order MCMC:",endtime_order, ",DAGs",endtime_DAGs, ",beta",endtime_beta, "\n")
    
    orderscore_prev <- unlist(example[[3]][length(example[[3]])]) #new order, previous beta
    orderscore_prop <- sum(orderscore_betas(n,c(1:n), weighted_betas_proposed, order_prev))
    
    # Log score of new DAG set under the old beta
    # nDAGoBeta_logscore <- Reduce("+", lapply(1:length(weights_proposed), function(k) incidence_logscore[[k]] * weights_proposed[k]))
    nDAGoBeta_logscore <- Reduce("+", lapply(1:length(weights_proposed), function(k) incidence_logscore[[k]]))
    # Log score of old DAG set under the new beta
    oDAGnBeta_logscore <- calculate_DAG_score(DAG_list = prev_DAG ,permy = order_prev, weights = NULL ,betas = weighted_betas_proposed)
    
    # Acceptance ratio
    # Test
    wA <- is_results$log_diff/nr_sample
    inv_wB <- sum(oDAGnBeta_logscore - totalscore_of_DAGs_prev)/nr_sample
    
    # wA <-  (totalscore_DAGs - nDAGoBeta_logscore)/nr_sample
    # inv_wB <- (sum(oDAGnBeta_logscore) - sum(totalscore_of_DAGs_prev))/nr_sample
    ratio <- exp(wA + inv_wB)
    # cat("ratio",ratio,"totalscore_DAGs",totalscore_DAGs, "totalscore_of_DAGs_prev", totalscore_of_DAGs_prev, "\n")
    # Test
    #ratio <- exp(totalscore_DAGs - totalscore_of_DAGs[i] + orderscore_prop -orderscore_prev + current_logscore -proposed_logscore)

    # Metropolis-Hastings acceptance step
    if(runif(1) < ratio){ 
      compress_DAG[[i+1]] <- is_results$compress_dag
      weighted_betas[[i+1]] <- weighted_betas_proposed
      order[[i+1]] <- permy
      count_accept[i] <- 1 # Accept 
      totalscore_of_DAGs_prev <- totalscore_DAGs
      #cat("ratio:", ratio, "\n")
      prev_DAG <- incidence_matrices
      cat("ratio",ratio,"wA",wA, "inv_wB", inv_wB, "\n")
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
              diffBiDAGs = diff_BiDAGs[-c(1:burin_iter)],
              totalscore_of_DAGs = totalscore_of_DAGs[-c(1:burin_iter)]
              )
         )
}

# Incorporation of the Auxiliary Variable Method in Bayesian Network Learning
aux_OneDAG <- function(n, iteration, order_iter, order = NULL, 
                       order_stepsize, moveprobs, base_score = 0, 
                       starting_dag = NULL, betas_init = NULL, skeleton = FALSE,
                       edgesposterior, burnin = 0.4) {
  
  elements <- c(1:n) # Define set of elements
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
    calcultion_betas_init <- calculateSN_ScoresArray_hash(starting_dag, k = 1, n, base_score = 0)
    betas_init <- calcultion_betas_init$allBetaScores[,,1]
    base_score <- 0
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
  sum_matrix <- matrix(0, nrow = n, ncol = n)
  edge_over_time <- array(0, dim = c(n, n, iter))
  
  init_sampled_DAGs <- lapply(1:nr_sample, function(x) samplescore(n, betas_init, order[[1]]))
  sampled_prev_matrix <- lapply(init_sampled_DAGs, function(dag) dag$incidence) 
  sampled_prev_logscore <- lapply(init_sampled_DAGs, function(dag) dag$logscore) 
  sampled_prev_DAG_BiDAG <- calculateSN_ScoresArray_hash(sampled_prev_matrix, k = length(sampled_prev_matrix), n, base_score)$target_DAG_score
  
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
    
    # Update beta matrix using the weights from sampled DAGs
    calculation_beta_values <- calculateSN_ScoresArray_hash(incidence_matrices, k = length(incidence_matrices) ,n, base_score = base_score)
    BiDAGscore_propose <- calculation_beta_values$target_DAG_score
    beta_values <- calculation_beta_values$allBetaScores
    
    # Update beta matrix using importance sampling
    is_results <- importance_DAG(DAGs = incidence_matrices, score_under_betas = incidence_logscore, target_scores = BiDAGscore_propose)
    weights_proposed <- is_results$importance_weights # normalized weights under old beta
    
    #New beta matrix using the normalised weights
    weighted_betas_proposed <- Reduce("+", lapply(1:length(weights_proposed), 
                                                  function(k) beta_values[,,k] * weights_proposed[k]))
    
    
    sampled_new_DAG <- lapply(1:nr_sample, function(x) samplescore(n, weighted_betas_proposed, permy))
    sampled_new_matrix <- lapply(sampled_new_DAG, function(dag) dag$incidence)
    sampled_new_logscore <- lapply(sampled_new_DAG, function(dag) dag$logscore)
    
    sampled_new_DAG_BiDAG <- calculateSN_ScoresArray_hash(sampled_new_matrix, k = length(sampled_new_matrix) ,n, base_score = base_score)$target_DAG_score
    
    
    # # Log score of old DAG set under the new beta
    nDAGoBeta_logscore <- general_calculate_DAG_score(DAG_list = incidence_matrices, weights = NULL ,
                                                      betas = weighted_betas[[i]], base_score = base_score)
    
    nDAGnBeta_logscore <- general_calculate_DAG_score(DAG_list = incidence_matrices, weights = NULL ,
                                                      betas = weighted_betas_proposed, base_score = base_score)
    
    
    BiDAG_ratio <- unlist(sampled_new_DAG_BiDAG) - unlist(sampled_prev_DAG_BiDAG)
    yB_ratio <- unlist(nDAGnBeta_logscore) - unlist(incidence_logscore) 
    sampled_ratio <- unlist(sampled_prev_logscore) - unlist(sampled_new_logscore)
    
    # cat("BiDAG_ratio",BiDAG_ratio,"yB_ratio", yB_ratio, "sampled_ratio", sampled_ratio, "\n")
    log_ratio <- BiDAG_ratio + yB_ratio + sampled_ratio
    ratio <- exp(log_ratio)
    
    # Metropolis-Hastings acceptance step
    if(runif(1) < ratio){ 
      compress_DAG[[i+1]] <- sampled_new_matrix[[1]]
      weighted_betas[[i+1]] <- weighted_betas_proposed
      order[[i+1]] <- permy
      count_accept[i] <- 1 # Accept 
      BiDAGscore_prev <- BiDAGscore_propose
      
      sampled_prev_DAG_BiDAG <- sampled_new_DAG_BiDAG
      sampled_prev_logscore <- sampled_new_logscore
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
      sum_matrix <- sum_matrix + compress_DAG[[i+1]]
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
              # trace_totalscore = trace_totalscore,
              # ,totalscore_of_DAGs = totalscore_of_DAGs[-c(1:burin_iter)]
  )
  )
}
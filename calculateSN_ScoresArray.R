#This Document contains the following functions:
#calculate_DAG_score
#calculateBetaScoresArray_hash
#importance_DAG

# Input: a set of DAGs as input
# Output : logscores of the DAGs in the set with given beta
calculate_DAG_score <- function(DAG_list, permy, weights ,betas, base_score, party = NULL, posy = NULL){
  sampledscore <- numeric()
  for (dagIndex in seq_along(DAG_list)) {
    incidence <- DAG_list[[dagIndex]]
    sampledscore[dagIndex] <- base_score
    
    for (child in 1:n){
      #Store the position of each element in the permutation (i.e. an inverse permutation)
      positions <- order(permy)
      
      if(!is.null(party) && !is.null(posy)){  #partitionMCMC
        partyelement<-posy[positions[child]]
        if(partyelement==length(party)){# no parents are allowed
        }else{
          parentnodes<-permy[which(posy > partyelement)]
          for (parent in parentnodes){
            sampledscore[dagIndex] <- sampledscore[dagIndex] + betas[parent, child]*incidence[parent, child]                      
          }
        }
      }else{                                    #orderMCMC
        if(positions[child]==n){ # no parents are allowed
        } else {
          #Nodes that are allowed to be its parents
          parentnodes <- permy[c((positions[child]+1):n)]
          for (parent in parentnodes){
            sampledscore[dagIndex] <- sampledscore[dagIndex] + betas[parent, child]*incidence[parent, child]                      
          }
        }
      }
    }
  }
  if(is.null(weights)){
    logscore <- sampledscore
  }else{
    logscore <- sampledscore*weights
  }
  return(logscore)
}

general_calculate_DAG_score <- function(DAG_list, weights ,betas, base_score){
  sampledscore <- numeric()
  for (dagIndex in seq_along(DAG_list)) {
    incidence <- DAG_list[[dagIndex]]
    sampledscore[dagIndex] <- base_score
    for (child in 1:n){
      for (parent in 1:n) {
        if (parent != child){
          sampledscore[dagIndex] <- sampledscore[dagIndex] + betas[parent, child]*incidence[parent, child]  
        }
      }
    }
  }
  if(is.null(weights)){
    logscore <- sampledscore
  }else{
    logscore <- sampledscore*weights
  }
  return(logscore)
}

calculate_final_score_mean <- function(log_scores){
  max_logscore <- max(log_scores) # Find the maximum log score to adjust the scores for numerical stability
  adjusted_log_scores <- log_scores - max_logscore
  result <- mean(exp(adjusted_log_scores))
  # Take the log of the result and add back the max log score
  final_score <- log(result) + max_logscore
  return(final_score)
}

# Function that computes the final score based on log scores.
# The 'operation' parameter determines whether to sum or average the scores.
calculate_final_score <- function(log_scores, operation) {
  # operation <- match.arg(c("sum", "mean")) # Ensure 'operation' is one of the allowed values
  max_logscore <- max(log_scores) # Find the maximum log score to adjust the scores for numerical stability
  adjusted_log_scores <- log_scores - max_logscore
  
  # Depending on the operation, sum or average the exponentiated adjusted scores
  if (operation == "sum") {
    result <- sum(exp(adjusted_log_scores))
  } else if (operation == "mean"){
    result <- mean(exp(adjusted_log_scores))
  } 
  # Take the log of the result and add back the max log score
  final_score <- log(result) + max_logscore
  
  return(final_score)
}
# 
# calculte_beta_matrix <- function(beta_matrix, weigths){
#   matrix <- matrix(0, nrow = n, ncol = n)
#   for (child in 1:n){
#     for (parent in 1:n){
#       matrix[parent, child] <- calculate_final_score(beta_matrix[parent, child, ], weigths, operation = "weigths" )
#     }
#   }
#   return(matrix)
# }


#importance_DAG 
importance_DAG <- function(DAGs, score_under_betas, target_scores){
  differ_score <- numeric()
  for (i in 1:length(DAGs)){
    differ_score[i] <- target_scores[i] - score_under_betas[[i]]
  }
  # To avoid numerical issues, subtract the max score before exponentiating
  # print(differ_score)
  max_score <- max(differ_score)
  adjusted_scores <- differ_score - max_score
  importance_weights <- exp(adjusted_scores)

  # Normalize the exponentiated scores to get the importance weights
  normalised_weights <- importance_weights / sum(importance_weights)
  
  ess_value <- (1 / sum(normalised_weights^2))/length(DAGs)
  
  compress_dag <- Reduce("+", lapply(1:length(normalised_weights), 
                                     function(k) DAGs[[k]] * normalised_weights[k]))
  
  return(list(unnorm_weights = importance_weights, importance_weights = normalised_weights, log_diff = differ_score, ess_value = ess_value, compress_dag = compress_dag))
}

DAGs_from_order <- function(order_list, nr_sample, beta_matrix){
  all_DAGs <- list()
  all_DAGs_logscore <- list()
  order_of_DAGs <- list()
  for(order in order_list){
    sampled_DAGs = lapply(1:nr_sample, function(x) samplescore(n, beta_matrix, order))
    DAGs_for_order <- lapply(sampled_DAGs, function(dag) dag$incidence) 
    logscores_DAGs <- lapply(sampled_DAGs, function(dag) dag$logscore) 
    # Append these DAGs to the main list
    all_DAGs <- c(all_DAGs, DAGs_for_order)
    all_DAGs_logscore <- c(all_DAGs_logscore,logscores_DAGs)
    order_of_DAGs <- c(order_of_DAGs, rep(list(order), nr_sample))
  }
  
  return(list(incidence = all_DAGs, logscore = all_DAGs_logscore, order = order_of_DAGs))
}

# This file contains the function for importance sampling

#importance_DAG_prev uses the the previous beta_matrix to calculate the DAG score
importance_DAG <- function(DAGs, score_under_betas){
  differ_score <- numeric()
  target_scores <- numeric()
  
  for (i in 1:length(DAGs)){
    target_score <- BiDAG::DAGscore(scoreParam, DAGs[[i]])
    beta_score <- score_under_betas[[i]]
    differ_score[i] <- target_score - beta_score
  }
  # To avoid numerical issues, subtract the max score before exponentiating
  max_score <- max(differ_score)
  exp_scores <- exp(differ_score - max_score)
  
  # Normalize the exponentiated scores to get the importance weights
  importance_weights <- exp_scores / sum(exp_scores)

  ess_value <- 1 / sum(importance_weights^2)
  
  compress_dag <- Reduce("+", lapply(1:length(importance_weights), 
                                     function(k) DAGs[[k]] * importance_weights[k]))

  return(list(importance_weights = importance_weights, ess_value = ess_value, compress_dag = compress_dag))
}


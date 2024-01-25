# This file contains the function for importance sampling

#importance_DAG_prev uses the the previous beta_matrix to calculate the DAG score
importance_DAG_prev <- function(DAGs, s_betas){
  differ_score <- numeric()
  target_scores <- numeric()
  for (i in 1:length(DAGs)){
    target_score <- BiDAG::DAGscore(scoreParam, DAGs[[i]])
    beta_score <- s_betas[[i]]
    #cat("beta_score",i, ":" , beta_score, "\n")
    differ_score[i] <- target_score - beta_score
  }
  # To avoid numerical issues, subtract the max score before exponentiating
  max_score <- max(differ_score)
  exp_scores <- exp(differ_score - max_score)
  
  # Normalize the exponentiated scores to get the importance weights
  importance_weights <- exp_scores / sum(exp_scores)
  #cat("importance_weights", ":" , importance_weights, "\n")
  ess_value <- 1 / sum(importance_weights^2)
  #cat("ess_value", ":" , ess_value, "\n")
  return(list(importance_weights = importance_weights, ess_value = ess_value))
}

#Importance Weights for the sampled DAGs
importance_DAG <- function(DAGs, betas){
  differ_score <- numeric()
  
  #betas = beta_values
  for (i in 1:length(DAGs)){
    target_score <- BiDAG::DAGscore(scoreParam, DAGs[[i]])
    #cat("target_score",i, ":" , target_score, "\n")
    beta_score <- DAGscore_under_betas(DAGs[[i]], betas[,,i])
    #cat("beta_score",i, ":" , beta_score, "\n")
    differ_score[i] <- target_score - beta_score
  }
  # To avoid numerical issues, subtract the max score before exponentiating
  max_score <- max(differ_score)
  exp_scores <- exp(differ_score - max_score)
  
  # Normalize the exponentiated scores to get the importance weights
  importance_weights <- exp_scores / sum(exp_scores)
  #cat("importance_weights", ":" , importance_weights, "\n")
  ess_value <- 1 / sum(importance_weights^2)
  
  #cat("ess_value", ":" , ess_value, "\n")
  
  return(list(importance_weights = importance_weights, ess_value = ess_value))
}


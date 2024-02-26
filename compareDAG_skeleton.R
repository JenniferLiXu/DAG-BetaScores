# Compare two DAGs in skeleton space
combineProbabilities <- function(dag_probs) {
  # Assuming dag_probs is a square matrix of directed edge probabilities
  undirected_probs <- dag_probs + t(dag_probs)
  #undirected_probs <- dag_probs | t(dag_probs)
  diag(undirected_probs) <- 0  # Remove self-loops if present
  return(undirected_probs)
}

CompareDAG <-function(sampled_dag, target_dag, skeleton = FALSE){
  num_nodes <- nrow(target_dag)
  edge_differences <- array(0, dim = c(num_nodes, num_nodes))
  
  if(skeleton == TRUE){
    sampled <- combineProbabilities(sampled_dag)
    target <- combineProbabilities(target_dag)
    
    # Compute the absolute differences
    for (row in 1:(num_nodes - 1)) {
      for (col in (row + 1):num_nodes) {
        edge_differences[row, col] <- abs(sampled[row, col] - target[row, col])
      }
    }
  }
  else{
    sampled <- sampled_dag
    target <- target_dag
    # Compute the absolute differences
    edge_differences <- abs(target_dag-sampled_dag)
  }
  
  return(edge_differences)
}




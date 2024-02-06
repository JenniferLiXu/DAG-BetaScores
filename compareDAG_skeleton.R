# Compare two DAGs in skeleton space
combineProbabilities <- function(dag_probs) {
  # Assuming dag_probs is a square matrix of directed edge probabilities
  undirected_probs <- dag_probs + t(dag_probs)
  diag(undirected_probs) <- 0  # Remove self-loops if present
  return(undirected_probs)
}

CompareDAG_skeleton <-function(sampled_dag_skeleton, target_dag){
  target_dag_skeleton <- combineProbabilities(target_dag)
  
  num_nodes <- nrow(target_dag)
  edge_differences <- array(0, dim = c(num_nodes, num_nodes))
  
  # Compute the absolute differences
  for (row in 1:(num_nodes - 1)) {
    for (col in (row + 1):num_nodes) {
      edge_differences[row, col] <- abs(sampled_dag_skeleton[row, col] - target_dag_skeleton[row, col])
    }
  }
  
  return(edge_differences)
}



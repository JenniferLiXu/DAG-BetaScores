#This function calculate beta for each possible parent node (e.g., node 2) 
#in relation to the child node (e.g., node 1) as follows

#Beta for a parent node is calculated by taking the score of the child node 
#when it has that parent node and other existing parents, 
#minus the score of the child node when it only has the other existing parents (excluding the specific parent node in question).

#This function will return a 3D array where allBetaScores[i, j, k] 
#represents the beta score of parent node i on child node j in the k-th sampled DAG.

calculateBetaScoresArray <- function(sampledDAGs, k, n) {
  # The array dimensions are [number of parents, number of children, number of sampled DAGs]
  allBetaScores <- array(NA, dim = c(n, n, k))
  
  for (dagIndex in seq_along(sampledDAGs)) {
    incidence <- sampledDAGs[[dagIndex]]
    
    for (childNode in 1:n) {
      #childNode <- 3
      # Identify the parent nodes for this child in the current DAG
      parentNodes <- which(incidence[, childNode] == 1)
      
      # Calculate the score for the child with all its current parents
      scoreWithAllParents <- BiDAG:::DAGcorescore(childNode, parentnodes = parentNodes,n, scoreParam)

      for (parentNode in 1:n) {
        #parentNode <- 1
        if (parentNode != childNode) {
          # Check if parentNode is actually a parent in this DAG
          if (parentNode %in% parentNodes) {
            # Remove parentNode and calculate the score
            WithoutparentNodes <- parentNodes[!parentNodes == parentNode]
            scoreWithoutParent <- BiDAG:::DAGcorescore(childNode, parentnodes = WithoutparentNodes, n, scoreParam)
            
            # Calculate beta score
            allBetaScores[parentNode, childNode, dagIndex] <- scoreWithAllParents - scoreWithoutParent
            
          } else {
            WithparentNodes <-c(parentNode, parentNodes)
            scoreWithParentNode <- BiDAG:::DAGcorescore(childNode, parentnodes = WithparentNodes, n, scoreParam)
    
            
            # Calculate the beta score
            allBetaScores[parentNode, childNode, dagIndex] <- scoreWithParentNode - scoreWithAllParents
            
            # Print the scoreWithParentNode
            # print(paste("Child node",childNode,
            #             " has parents:",paste(parentNodes, collapse = ", "),
            #             " with Added node:",parentNode,
            #             " scoreWithParentNode is", scoreWithParentNode,
            #             " scoreWithAllParents is", scoreWithAllParents))
            
          }
        } else {
          # Beta score is not applicable for self (node cannot be its own parent)
          allBetaScores[parentNode, childNode, dagIndex] <- NA
        }
      }
    }
  }
  
  return(allBetaScores)
}





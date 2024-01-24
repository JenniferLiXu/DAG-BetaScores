#This function calculate beta for each possible parent node (e.g., node 2) 
#in relation to the child node (e.g., node 1) as follows

#Beta for a parent node j : 
#beta[j,i] = (score of i with parents j and x, y, z) - (score of i with parents x, y, z)

#This function will return a 3D array where allBetaScores[j, i, k] 
#represents the beta score of parent node j on child node i in the k-th sampled DAG.

calculateBetaScoresArray <- function(sampledDAGs, k, n) {
  # The array dimensions are [number of parents, number of children, number of sampled DAGs]
  allBetaScores <- array(NA, dim = c(n, n, k))

  for (dagIndex in seq_along(sampledDAGs)) {
    incidence <- sampledDAGs[[dagIndex]]
    for (childNode in 1:n) {
      # Identify the parent nodes for this child in the current DAG
      parentNodes <- which(incidence[, childNode] == 1)

      # Calculate the score for the child with all its current parents
      scoreWithAllParents <- BiDAG:::DAGcorescore(childNode, parentnodes = parentNodes, scoreParam$n, scoreParam)

      for (parentNode in 1:n) {
        if (parentNode != childNode) {
          # Check if parentNode is actually a parent in this DAG
          if (parentNode %in% parentNodes) {
            # Remove parentNode and calculate the score
            withoutParentNodes <- parentNodes[!parentNodes == parentNode]
            scoreWithoutNode <- BiDAG:::DAGcorescore(childNode, parentnodes = withoutParentNodes, scoreParam$n, scoreParam)

            # Calculate beta score
            allBetaScores[parentNode, childNode, dagIndex] <- scoreWithAllParents - scoreWithoutNode

          } else {
            withParentNodes <-c(parentNodes, parentNode)
            scoreWithNode <- BiDAG:::DAGcorescore(childNode, parentnodes = withParentNodes, scoreParam$n, scoreParam)

            # Calculate the beta score
            allBetaScores[parentNode, childNode, dagIndex] <- scoreWithNode - scoreWithAllParents

          }
        } else {
          # Beta score is not applicable for self (node cannot be its own parent)
          allBetaScores[parentNode, childNode, dagIndex] <- 0
        }
      }
    }
  }
  return(allBetaScores)
}

# calculateBetaScoresArray <- function(sampledDAGs, k, n) {
#   # The array dimensions are [number of parents, number of children, number of sampled DAGs]
#   allBetaScores <- array(NA, dim = c(n, n, k))
# 
#   for (dagIndex in seq_along(sampledDAGs)) {
#     incidence <- sampledDAGs[[dagIndex]]
#     for (childNode in 1:n){
#       # Identify the parent nodes for this child in the current DAG
#       parentNodes <- which(incidence[, childNode] == 1)
# 
#       for (parentNode in 1:n) {
#         if (length(parentNodes) == 0){
#           # Calculate the score for the child node if they only have one parent
#           allBetaScores[parentNode, childNode, dagIndex] <- BiDAG:::DAGcorescore(childNode, parentnodes = parentNode, scoreParam$n, scoreParam)
#         }else{
#           # Calculate the score for the child with all its current parents
#           scoreWithAllParents <- BiDAG:::DAGcorescore(childNode, parentnodes = parentNodes, scoreParam$n, scoreParam)
#           if (parentNode != childNode){
#             # Check if parentNode is actually a parent in this DAG
#             if (parentNode %in% parentNodes){
#               # Remove parentNode and calculate the score
#               withoutParentNodes <- parentNodes[!parentNodes == parentNode]
#               scoreWithoutNode <- BiDAG:::DAGcorescore(childNode, parentnodes = withoutParentNodes, scoreParam$n, scoreParam)
# 
#               # Calculate beta score
#               allBetaScores[parentNode, childNode, dagIndex] <- scoreWithAllParents - scoreWithoutNode
#             }else{
#               withParentNodes <-c(parentNodes, parentNode)
#               scoreWithNode <- BiDAG:::DAGcorescore(childNode, parentnodes = withParentNodes, scoreParam$n, scoreParam)
# 
#               # Calculate the beta score
#               allBetaScores[parentNode, childNode, dagIndex] <- scoreWithNode - scoreWithAllParents
#             }
#           }else{
#             # Beta score is not applicable for self (node cannot be its own parent)
#             allBetaScores[parentNode, childNode, dagIndex] <- 0
#           }
#         }
# 
#       }
#     }
#   }
#   return(allBetaScores)
# }


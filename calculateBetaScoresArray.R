#This Document contains the following functions:

#calculate_DAG_score
#calculateBetaScoresArray_hash
#importance_DAG

#This function calculate beta for each possible parent node (e.g., node 2) 
#in relation to the child node (e.g., node 1) as follows

#Beta for a parent node j : 
#beta[j,i] = (score of i with parents j and x, y, z) - (score of i with parents x, y, z)

# Input: a set of DAGs as input
# Output : total logscores, by summing  weights times ( logscore of the DAGs under the beta_matrix )
calculate_DAG_score <- function(DAG_list, permy, weights ,betas, party = NULL, posy = NULL){
  sampledscore <- numeric()
  for (dagIndex in seq_along(DAG_list)) {
    incidence <- DAG_list[[dagIndex]]
    sampledscore[dagIndex] <- base_score
    
    for (child in 1:n){
      #Store the position of each element in the permutation (i.e. an inverse permutation)
      positions <- order(permy)
      
      if(!is.null(party) && !is.null(posy)){.  #partitionMCMC
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

####### Hash Tables ########
# Initialize two environments as hash tables
parentSetBetaCache <- new.env(hash = TRUE, parent = emptyenv())
individualParentBetaCache <- new.env(hash = TRUE, parent = emptyenv()) # Column of Beta values for each individual parents on the child node

# Generate a unique key for each parent set and child node combination
generateKey <- function(childNode, parentSet) {
  # Sort parentSet to ensure consistency, then concatenate with childNode
  paste0(childNode, "-", toString(sort(parentSet)))
}


# Check if a key exists and return the stored value if it does
getStoredValue <- function(hashTable, key) {
  if(exists(key, envir=hashTable)) {
    return(hashTable[[key]])
  }else{
    return(NULL)
  }
}

# Store a value in a cache
storeValue <- function(env, key, value) {
  env[[key]] <- value
}

# Compute the column of beta values for a child node given its current parent set
computeAndCacheBetaValues <- function(childNode, parentSet, n) {
  count_DAGcore <- 0
  # Generate key for the current parent set
  
  parentSetKey <- generateKey(childNode, parentSet)
  
  # Attempt to retrieve pre-computed beta values for the individual parents
  betaValues <- getStoredValue(individualParentBetaCache, parentSetKey)
  
  if (!is.null(betaValues)) {
    return(list(betaValues = betaValues[-(n+1)], count_DAGcore = count_DAGcore, DAG_parentset = betaValues[n+1])) # Cached beta values found, return them
  } else {
    # Beta values not found in cache, compute them
    # Ensure the score for the current parent set is cached
    parentSetScore <- getStoredValue(parentSetBetaCache, parentSetKey)
    betaValues <- rep(0, times = (n + 1)) # Store the beta values(colunm) + the score of the current parent set!
    if (is.null(parentSetScore)){
      # Compute and cache the score for the current parent set
      parentSetScore <- BiDAG:::DAGcorescore(childNode, parentnodes = parentSet, scoreParam$n, scoreParam)# Compute score for childNode with parentSet
      count_DAGcore <- count_DAGcore+1
      #print(parentSetKey)
      storeValue(parentSetBetaCache, parentSetKey, parentSetScore)
      betaValues[n+1] <- parentSetScore
    }

    # Compute beta values for individual parents by adding/removing them
    
    for (parentNode in 1:n) {
      if (parentNode != childNode) {
        # Check if parentNode is actually a parent in this DAG
        if (parentNode %in% parentSet) {
          # Remove parentNode and calculate the score
          reducedParentSet <- parentSet[!parentSet == parentNode]
          reducedKey <- generateKey(childNode, reducedParentSet)
          reducedScore <- getStoredValue(parentSetBetaCache, reducedKey)
          
          if (is.null(reducedScore)) {
            reducedScore <- BiDAG:::DAGcorescore(childNode, parentnodes = reducedParentSet, scoreParam$n, scoreParam)# Compute score for childNode with reducedParentSet
            count_DAGcore <- count_DAGcore+1
            #print(reducedKey)
            storeValue(parentSetBetaCache, reducedKey, reducedScore)
          }
          
          # Calculate beta value as the difference in scores
          betaValue <- parentSetScore - reducedScore
          betaValues[parentNode] <- betaValue
          
        } else {
          withParentNodes <-c(parentSet, parentNode)
          withParentNodesKey <- generateKey(childNode, withParentNodes)
          withParentNodesScore <- getStoredValue(parentSetBetaCache, withParentNodesKey)
          
          if (is.null(withParentNodesScore)) {
            withParentNodesScore <- BiDAG:::DAGcorescore(childNode, parentnodes = withParentNodes, scoreParam$n, scoreParam)# Compute score for childNode with reducedParentSet
            count_DAGcore <- count_DAGcore+1
            #print(withParentNodesKey)
            storeValue(parentSetBetaCache, withParentNodesKey, withParentNodesScore)
          }
          
          # Calculate beta value as the difference in scores
          betaValue <- withParentNodesScore - parentSetScore
          betaValues[parentNode] <- betaValue
        }
        
      } else {
        # Beta score is not applicable for self (node cannot be its own parent)
        betaValues[parentNode] <- 0
      }
    }
    # Cache the computed beta values for future reference
    storeValue(individualParentBetaCache, parentSetKey, betaValues)
    
    return(list(betaValues = betaValues[-(n+1)], count_DAGcore = count_DAGcore, DAG_parentset = betaValues[(n+1)])) 
  }
}
 
#This function will return a 3D array where allBetaScores[j, i, k] 
#represents the beta score of parent node j on child node i in the k-th sampled DAG.
calculateBetaScoresArray_hash <- function(sampledDAGs, k, n){
  count <-  0 
  # The array dimensions are [number of parents, number of children, number of sampled DAGs]
  allBetaScores <- array(NA, dim = c(n, n, k))
  incidence_scores <- numeric()
  for (dagIndex in seq_along(sampledDAGs)) {
    incidence <- sampledDAGs[[dagIndex]]
    incidence_score <- 0
    for (childNode in 1:n){
      # Identify the parent nodes for this child in the current DAG
      parentNodes <- which(incidence[, childNode] == 1)
      calculation_column_beta <- computeAndCacheBetaValues(childNode = childNode, parentSet = parentNodes, n)
      column_beta <- calculation_column_beta$betaValues
      count <- count + calculation_column_beta$count_DAGcore
      incidence_score <- incidence_score + calculation_column_beta$DAG_parentset
      
      allBetaScores[ , childNode, dagIndex] <- array(column_beta, dim = c(length(column_beta), 1, 1))
    }
    incidence_scores[dagIndex] <- incidence_score
  }
  return(list(allBetaScores = allBetaScores, count_DAGcore = count, target_DAG_score = incidence_scores))
}

# DAGs = incidence_matrices
# score_under_betas = incidence_logscore
# target_scores = calculation_beta_values$target_DAG_score


#importance_DAG 
importance_DAG <- function(DAGs, score_under_betas, target_scores){
  differ_score <- numeric()
  
  starttime_DAG <-proc.time()
  for (i in 1:length(DAGs)){
    differ_score[i] <- target_scores[i] - score_under_betas[[i]]
  }
  endtime_DAG <- proc.time()
  endtime_DAG <- endtime_DAG - starttime_DAG 
  
  # To avoid numerical issues, subtract the max score before exponentiating
  max_score <- max(differ_score)
  importance_weights <- exp(differ_score - max_score)
  log_diff <- sum(differ_score)
  # importance_weights <- exp(differ_score)
  # cat("importance_weights",importance_weights, "\n")
  
  # Normalize the exponentiated scores to get the importance weights
  normalised_weights <- importance_weights / sum(importance_weights)
  ess_value <- 1 / sum(importance_weights^2)
  
  compress_dag <- Reduce("+", lapply(1:length(normalised_weights), 
                                     function(k) DAGs[[k]] * normalised_weights[k]))
  
  return(list(importance_weights = normalised_weights, log_diff = log_diff, ess_value = ess_value, compress_dag = compress_dag))
}


# Example usage
# Assuming you're about to calculate values for child node 2 with parents 1 and 3
# childNode <- 2
# parentSet <- c(1,3)
# key <- generateKey(childNode, parentSet)
# 
# computeAndCacheBetaValues(childNode, parentSet, n)
# calculateBetaScoresArray_hash(results_seed123$DAGs, length(results_seed123$DAGs),n)


# ####### Another version of using hash table ########
# # Compute the column of beta values for a child node given its current parent set
# computeAndCacheBetaValues_short <- function(childNode, parentSet, n) {
#   count_DAGcore <- 0
#   # Generate key for the current parent set
#   parentSetKey <- generateKey(childNode, parentSet)
# 
#   # Beta values not found in cache, compute them
#   # Ensure the score for the current parent set is cached
#   parentSetScore <- getStoredValue(parentSetBetaCache, parentSetKey)
#   if (is.null(parentSetScore)){
#     # Compute and cache the score for the current parent set
#     parentSetScore <- BiDAG:::DAGcorescore(childNode, parentnodes = parentSet, scoreParam$n, scoreParam)# Compute score for childNode with parentSet
#     count_DAGcore <- count_DAGcore+1
#     storeValue(parentSetBetaCache, parentSetKey, parentSetScore)
#     #print(parentSetKey)
#     
#   }
#   
#   # Compute beta values for individual parents by adding/removing them
#   betaValues <- rep(0, times = n)
#   for (parentNode in 1:n) {
#     if (parentNode != childNode) {
#       # Check if parentNode is actually a parent in this DAG
#       if (parentNode %in% parentSet) {
#         # Remove parentNode and calculate the score
#         reducedParentSet <- parentSet[!parentSet == parentNode]
#         reducedKey <- generateKey(childNode, reducedParentSet)
#         reducedScore <- getStoredValue(parentSetBetaCache, reducedKey)
#         
#         if (is.null(reducedScore)) {
#           reducedScore <- BiDAG:::DAGcorescore(childNode, parentnodes = reducedParentSet, scoreParam$n, scoreParam)# Compute score for childNode with reducedParentSet
#           count_DAGcore <- count_DAGcore+1
#           storeValue(parentSetBetaCache, reducedKey, reducedScore)
#           #print(reducedKey)
#         }
#         
#         # Calculate beta value as the difference in scores
#         betaValue <- parentSetScore - reducedScore
#         betaValues[parentNode] <- betaValue
#         
#       } else {
#         withParentNodes <-c(parentSet, parentNode)
#         withParentNodesKey <- generateKey(childNode, withParentNodes)
#         withParentNodesScore <- getStoredValue(parentSetBetaCache, withParentNodesKey)
#         
#         if (is.null(withParentNodesScore)) {
#           withParentNodesScore <- BiDAG:::DAGcorescore(childNode, parentnodes = withParentNodes, scoreParam$n, scoreParam)# Compute score for childNode with reducedParentSet
#           count_DAGcore <- count_DAGcore+1
#           storeValue(parentSetBetaCache, withParentNodesKey, withParentNodesScore)
#           #print(withParentNodesKey)
#         }
#         
#         # Calculate beta value as the difference in scores
#         betaValue <- withParentNodesScore -  parentSetScore
#         betaValues[parentNode] <- betaValue
#       }
#       
#     } else {
#       # Beta score is not applicable for self (node cannot be its own parent)
#       betaValues[parentNode] <- 0
#     }
#   }
#   # Cache the computed beta values for future reference
#   storeValue(individualParentBetaCache, parentSetKey, betaValues)
#   
#   return(list(betaValues = betaValues, count_DAGcore=count_DAGcore))
# }
# 
# calculateBetaScoresArray_hash_short <- function(sampledDAGs, k, n){
#   count <-  0
#   # The array dimensions are [number of parents, number of children, number of sampled DAGs]
#   allBetaScores <- array(NA, dim = c(n, n, k))
# 
#   for (dagIndex in seq_along(sampledDAGs)) {
#     incidence <- sampledDAGs[[dagIndex]]
#     for (childNode in 1:n){
#       # Identify the parent nodes for this child in the current DAG
#       parentNodes <- which(incidence[, childNode] == 1)
#       
#       calculation_column_beta <- computeAndCacheBetaValues_short(childNode = childNode, parentSet = parentNodes, n)
#       column_beta <- calculation_column_beta$betaValues
#       count <- count + calculation_column_beta$count_DAGcore
#       
#       allBetaScores[ , childNode, dagIndex] <- array(column_beta, dim = c(length(column_beta), 1, 1))
#     }
#   }
#   return(list(allBetaScores = allBetaScores, count_DAGcore = count))
# }





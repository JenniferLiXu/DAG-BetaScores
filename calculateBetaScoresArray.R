#This function calculate beta for each possible parent node (e.g., node 2) 
#in relation to the child node (e.g., node 1) as follows

#Beta for a parent node j : 
#beta[j,i] = (score of i with parents j and x, y, z) - (score of i with parents x, y, z)

#This function will return a 3D array where allBetaScores[j, i, k] 
#represents the beta score of parent node j on child node i in the k-th sampled DAG.

calculateBetaScoresArray <- function(sampledDAGs, k, n) {
  # The array dimensions are [number of parents, number of children, number of sampled DAGs]
  allBetaScores <- array(NA, dim = c(n, n, k))
  positions <- order(permy) 

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
  
  logscore <- sum(sampledscore*weights)
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
  if(exists("key", envir=hashTable)) {
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
  # Generate key for the current parent set
  parentSetKey <- generateKey(childNode, parentSet)
  
  # Attempt to retrieve pre-computed beta values for the individual parents
  betaValues <- getStoredValue(individualParentBetaCache, parentSetKey)
  
  if (!is.null(betaValues)) {
    return(betaValues) # Cached beta values found, return them
  } else {
    # Beta values not found in cache, compute them
    # Ensure the score for the current parent set is cached
    parentSetScore <- getStoredValue(parentSetBetaCache, parentSetKey)
    if (is.null(parentSetScore)){
      # Compute and cache the score for the current parent set
      parentSetScore <- BiDAG:::DAGcorescore(childNode, parentnodes = parentSet, scoreParam$n, scoreParam)# Compute score for childNode with parentSet
      storeValue(parentSetBetaCache, parentSetKey, parentSetScore)
    }

    # Compute beta values for individual parents by adding/removing them
    betaValues <- rep(0, times = n)
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
    
    return(betaValues)
  }
}

calculateBetaScoresArray_hash <- function(sampledDAGs, k, n){
  # The array dimensions are [number of parents, number of children, number of sampled DAGs]
  allBetaScores <- array(NA, dim = c(n, n, k))

  for (dagIndex in seq_along(sampledDAGs)) {
    incidence <- sampledDAGs[[dagIndex]]
    for (childNode in 1:n){
      # Identify the parent nodes for this child in the current DAG
      parentNodes <- which(incidence[, childNode] == 1)
      column_beta <- computeAndCacheBetaValues(childNode = childNode, parentSet = parentNodes, n)
      
      allBetaScores[ , childNode, dagIndex] <- array(column_beta, dim = c(length(column_beta), 1, 1))
    }
  }
  return(allBetaScores)
}

# Example usage
# Assuming you're about to calculate values for child node 2 with parents 1 and 3
# childNode <- 2
# parentSet <- c(1,3)
# key <- generateKey(childNode, parentSet)
# 
# computeAndCacheBetaValues(childNode, parentSet, n)
# calculateBetaScoresArray_hash(results_seed123$DAGs, length(results_seed123$DAGs),n)


####### Another version of using hash table ########
# Compute the column of beta values for a child node given its current parent set
computeAndCacheBetaValues_short <- function(childNode, parentSet, n) {
  # Generate key for the current parent set
  parentSetKey <- generateKey(childNode, parentSet)

  # Beta values not found in cache, compute them
  # Ensure the score for the current parent set is cached
  parentSetScore <- getStoredValue(parentSetBetaCache, parentSetKey)
  if (is.null(parentSetScore)){
    # Compute and cache the score for the current parent set
    parentSetScore <- BiDAG:::DAGcorescore(childNode, parentnodes = parentSet, scoreParam$n, scoreParam)# Compute score for childNode with parentSet
    storeValue(parentSetBetaCache, parentSetKey, parentSetScore)
  }
  
  # Compute beta values for individual parents by adding/removing them
  betaValues <- rep(0, times = n)
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
          storeValue(parentSetBetaCache, withParentNodesKey, withParentNodesScore)
        }
        
        # Calculate beta value as the difference in scores
        betaValue <- withParentNodesScore -  parentSetScore
        betaValues[parentNode] <- betaValue
      }
      
    } else {
      # Beta score is not applicable for self (node cannot be its own parent)
      betaValues[parentNode] <- 0
    }
  }
  # Cache the computed beta values for future reference
  storeValue(individualParentBetaCache, parentSetKey, betaValues)
  
  return(betaValues)
}

calculateBetaScoresArray_hash_short <- function(sampledDAGs, k, n){
  # The array dimensions are [number of parents, number of children, number of sampled DAGs]
  allBetaScores <- array(NA, dim = c(n, n, k))
  
  for (dagIndex in seq_along(sampledDAGs)) {
    incidence <- sampledDAGs[[dagIndex]]
    for (childNode in 1:n){
      # Identify the parent nodes for this child in the current DAG
      parentNodes <- which(incidence[, childNode] == 1)
      column_beta <- computeAndCacheBetaValues_short(childNode = childNode, parentSet = parentNodes, n)
      
      allBetaScores[ , childNode, dagIndex] <- array(column_beta, dim = c(length(column_beta), 1, 1))
    }
  }
  return(allBetaScores)
}

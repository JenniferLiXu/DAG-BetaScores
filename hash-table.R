####### Hash Tables ########
# Initialize two environments as hash tables
parentSetSNCache <- new.env(hash = TRUE, parent = emptyenv()) 
individualParentSNCache <- new.env(hash = TRUE, parent = emptyenv()) # Column of Beta values for each individual parents on the child node

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


# Compute the column of single node values for a child node given its current parent set
computeAndCacheSN_Values <- function(childNode, parentSet, n) {
  count_DAGcore <- 0
  # Generate key for the current parent set
  parentSetKey <- generateKey(childNode, parentSet)
  
  # Attempt to retrieve pre-computed single node values for the individual parents
  SN_Values <- getStoredValue(individualParentSNCache, parentSetKey)
  if (!is.null(SN_Values)) {
    return(list(SN_Values = SN_Values[-(n+1)], count_DAGcore = count_DAGcore, DAG_parentset = SN_Values[n+1])) # Cached single node values found, return them
  } else {
    # single node values not found in cache, compute them
    parentSetScore <- getStoredValue(parentSetSNCache, parentSetKey)
    SN_Values <- numeric() # Store the beta values(colunm) + the score of the current parent set!
    if (is.null(parentSetScore)){
      # Compute and cache the score for the current parent set
      parentSetScore <- BiDAG:::DAGcorescore(childNode, parentnodes = parentSet, scoreParam$n, scoreParam)# Compute score for childNode with parentSet
      count_DAGcore <- count_DAGcore+1
      storeValue(parentSetSNCache, parentSetKey, parentSetScore)
    }
    SN_Values[n+1] <- parentSetScore # the score of the current parent set
    # Compute beta values for individual parents by adding/removing them
    for (parentNode in 1:n) {
      if (parentNode != childNode) {
        # Check if parentNode is actually a parent in this DAG
        if (parentNode %in% parentSet) {
          # Remove parentNode and calculate the score
          reducedParentSet <- parentSet[!parentSet == parentNode]
          reducedKey <- generateKey(childNode, reducedParentSet)
          reducedScore <- getStoredValue(parentSetSNCache, reducedKey)
          if (is.null(reducedScore)) {
            reducedScore <- BiDAG:::DAGcorescore(childNode, parentnodes = reducedParentSet, scoreParam$n, scoreParam)# Compute score for childNode with reducedParentSet
            count_DAGcore <- count_DAGcore+1
            storeValue(parentSetSNCache, reducedKey, reducedScore)
          }
          # Calculate beta value as the difference in scores
          betaValue <- parentSetScore - reducedScore
          SN_Values[parentNode] <- betaValue
          
        } else {
          withParentNodes <-c(parentSet, parentNode)
          withParentNodesKey <- generateKey(childNode, withParentNodes)
          withParentNodesScore <- getStoredValue(parentSetSNCache, withParentNodesKey)
          if (is.null(withParentNodesScore)) {
            withParentNodesScore <- BiDAG:::DAGcorescore(childNode, parentnodes = withParentNodes, scoreParam$n, scoreParam)# Compute score for childNode with reducedParentSet
            count_DAGcore <- count_DAGcore+1
            storeValue(parentSetSNCache, withParentNodesKey, withParentNodesScore)
          }
          # Calculate beta value as the difference in scores
          betaValue <- withParentNodesScore - parentSetScore
          SN_Values[parentNode] <- betaValue
        }
      } else {
        # Beta score is not applicable for self (node cannot be its own parent)
        SN_Values[parentNode] <- 0
      }
    }
    # Cache the computed beta values for future reference
    storeValue(individualParentSNCache, parentSetKey, SN_Values)
    
    return(list(SN_Values = SN_Values[-(n+1)], count_DAGcore = count_DAGcore, DAG_parentset = SN_Values[(n+1)])) 
  }
}

#This function will return a 3D array where allBetaScores[j, i, k] 
#represents the beta score of parent node j on child node i in the k-th sampled DAG.
calculateSN_ScoresArray_hash <- function(sampledDAGs, k, n, base_score){
  count <-  0 
  # The array dimensions are [number of parents, number of children, number of sampled DAGs]
  allBetaScores <- array(NA, dim = c(n, n, k))
  incidence_scores <- numeric()
  for (dagIndex in seq_along(sampledDAGs)) {
    incidence <- sampledDAGs[[dagIndex]]
    incidence_score <- base_score
    for (childNode in 1:n){
      # Identify the parent nodes for this child in the current DAG
      parentNodes <- which(incidence[, childNode] == 1)
      calculation_column_beta <- computeAndCacheSN_Values(childNode = childNode, parentSet = parentNodes, n)
      column_beta <- calculation_column_beta$SN_Values
      count <- count + calculation_column_beta$count_DAGcore
      # cat("DAG_parentset:",calculation_column_beta$DAG_parentset, "childNode", childNode ,",parentNodes",parentNodes, "\n")
      incidence_score <- incidence_score + calculation_column_beta$DAG_parentset
      allBetaScores[ , childNode, dagIndex] <- array(column_beta, dim = c(length(column_beta), 1, 1))
      # print(allBetaScores[ , childNode, dagIndex])
    }
    incidence_scores[dagIndex] <- incidence_score
  }
  return(list(allBetaScores = allBetaScores, count_DAGcore = count, target_DAG_score = incidence_scores))
}


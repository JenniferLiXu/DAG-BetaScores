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

#Beta for a parent node j : 
#beta[j,i] = (score of i with parents j, x, y, z) - (score of i with parents x, y, z)

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
    parentSetScore <- getStoredValue(parentSetBetaCache, parentSetKey)
    betaValues <- numeric() # Store the beta values(colunm) + the score of the current parent set!
    if (is.null(parentSetScore)){
      # Compute and cache the score for the current parent set
      parentSetScore <- BiDAG:::DAGcorescore(childNode, parentnodes = parentSet, scoreParam$n, scoreParam)# Compute score for childNode with parentSet
      count_DAGcore <- count_DAGcore+1
      storeValue(parentSetBetaCache, parentSetKey, parentSetScore)
    }
    betaValues[n+1] <- parentSetScore # the score of the current parent set
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
calculateBetaScoresArray_hash <- function(sampledDAGs, k, n, base_score){
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
      calculation_column_beta <- computeAndCacheBetaValues(childNode = childNode, parentSet = parentNodes, n)
      column_beta <- calculation_column_beta$betaValues
      count <- count + calculation_column_beta$count_DAGcore
      # cat("DAG_parentset:",calculation_column_beta$DAG_parentset, "childNode", childNode ,",parentNodes",parentNodes, "\n")
      incidence_score <- incidence_score + calculation_column_beta$DAG_parentset
      allBetaScores[ , childNode, dagIndex] <- array(column_beta, dim = c(length(column_beta), 1, 1))
    }
    incidence_scores[dagIndex] <- incidence_score
  }
  return(list(allBetaScores = allBetaScores, count_DAGcore = count, target_DAG_score = incidence_scores))
}

#importance_DAG 
importance_DAG <- function(DAGs, score_under_betas, target_scores){
  differ_score <- numeric()
  for (i in 1:length(DAGs)){
    differ_score[i] <- target_scores[i] - score_under_betas[[i]]
  }
  # To avoid numerical issues, subtract the max score before exponentiating
  max_score <- max(differ_score)
  importance_weights <- exp(differ_score - max_score)

  # Normalize the exponentiated scores to get the importance weights
  normalised_weights <- importance_weights / sum(importance_weights)
  ess_value <- 1 / sum(importance_weights^2)
  
  compress_dag <- Reduce("+", lapply(1:length(normalised_weights), 
                                     function(k) DAGs[[k]] * normalised_weights[k]))
  
  return(list(importance_weights = normalised_weights, log_diff = differ_score, ess_value = ess_value, compress_dag = compress_dag))
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
    order_of_DAGs <- order_of_DAGs <- c(order_of_DAGs, rep(list(order), nr_sample))
  }
  
  return(list(incidence = all_DAGs, logscore = all_DAGs_logscore, order = order_of_DAGs))
}

# In each order, we sample a set of DAGs, take thier weigted average and their corresponding weighted beta matrix
# We take the average of all beta matrices and the average of compressed DAGs from each order
DAGs_from_order_ver3 <- function(order_list, nr_sample, beta_matrix){
  base_score = 0
  all_DAGs <- list()
  all_DAGs_logscore <- list()
  order_of_DAGs <- list()
  betas_final <- matrix(0, nrow = n, ncol = n)
  total_DAG <- matrix(0, nrow = n, ncol = n)
  for(order in order_list){
    sampled_DAGs = lapply(1:nr_sample, function(x) samplescore(n, beta_matrix, order))
    DAGs_for_order <- lapply(sampled_DAGs, function(dag) dag$incidence) 
    logscores_DAGs <- lapply(sampled_DAGs, function(dag) dag$logscore) 
    
    calculation_beta_values <- calculateBetaScoresArray_hash(DAGs_for_order, k = length(DAGs_for_order) ,n, base_score = base_score)
    BiDAGscore_propose_list <- calculation_beta_values$target_DAG_score
    BiDAGscore_propose <-  calculate_final_score(BiDAGscore_propose_list, operation = "mean")
    
    beta_values <- calculation_beta_values$allBetaScores
    is_results <- importance_DAG(DAGs = DAGs_for_order, score_under_betas = logscores_DAGs, target_scores = BiDAGscore_propose_list)
    
    DAG_ave <- is_results$compress_dag
    weights_proposed <- is_results$importance_weights # normalised weights under old beta
    
    #New beta matrix using the normalised weights
    weighted_betas_proposed <- Reduce("+", lapply(1:length(weights_proposed), function(k) beta_values[,,k] * weights_proposed[k]))
    
    betas_final = betas_final + weighted_betas_proposed
    total_DAG = total_DAG + DAG_ave
    # Append these DAGs to the main list
    all_DAGs <- c(all_DAGs, DAGs_for_order)
    all_DAGs_logscore <- c(all_DAGs_logscore,logscores_DAGs)
    order_of_DAGs <- order_of_DAGs <- c(order_of_DAGs, rep(list(order), nr_sample))
  }
  betas_final_ave <- betas_final/length(order_list)
  total_DAG_ave <- total_DAG/length(order_list)
  
  return(list(incidence = all_DAGs, logscore = all_DAGs_logscore, order = order_of_DAGs, 
              betas_final_ave = betas_final_ave, total_DAG_ave = total_DAG_ave))
}

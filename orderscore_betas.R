# This function scores an order given the beta values

#Input: an ordering of nodes (permutation), the beta values, which DAG, the dataset we are analyzing, etc
#Output: returns a score for that ordering

# For each node in the ordering, calculate its score based on the beta values of its parents.
# Aggregate individual scores to get a total score for the ordering.


#beta value is a 3D array where betas[i, j, k] 
#represents the beta score of parent node i on child node j in the k-th sampled DAG.
#betas is a 2D array [i, j] 

orderscore_betas <- function(n,scorenodes, betas ,permy) {
  
  nodescores<-array(NA, dim = n)
  
  #Store the position of each element in the permutation (i.e. an inverse permutation)
  positions <- order(permy) 
  
  for (i in scorenodes){
    if(positions[i]==n){ # no parents are allowed
      nodescores[i]<- 1 # there is only one score
    } else {
      if(positions[i]>1){
        #Nodes that are allowed to be its parents
        parentnodes<-permy[c((positions[i]+1):n)]

        # For each node i, calculate the scores
        #nodescores[i] <- sum(log(1 + exp(betas[i, parentnodes])))
        node_scores <- numeric(length(parentnodes))  # Initialize node scores
        for (j in 1:length(parentnodes)) {
          a <- 0  # This is the '0' in your 'log(1 + exp(...))' expression
          b <- betas[i, parentnodes[j]]
          
          c <- max(a, b)
          d <- min(a, b)
          
          # Numerically stable computation
          node_scores[j] <- c + log1p(exp(d - c))
        }
        nodescores[i] <- sum(node_scores)
        
      } else{ # all parents are allowed
        nodescores[i] <- sum(log(1 + exp(betas[i, -i])))
      }
    }
  }
  
  
  scores <- nodescores
  #scores$orderscores<-sum(nodescores)
  return(scores)
  
}


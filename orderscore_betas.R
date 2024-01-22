# This function scores an order given the beta values

#Input: an ordering of nodes (permutation), the beta values
#Output: returns a score for that ordering

# For each node in the ordering, calculate its score based on the beta values of its parents.
# Aggregate individual scores to get a total score for the ordering.

#betas is a 2D array [parent, child] 
orderscore_betas <- function(n,scorenodes, betas, permy) {
  nodescores<- rep(0,n)
  
  #Store the position of each element in the permutation (i.e. an inverse permutation)
  positions <- order(permy) 
  
  for (i in scorenodes){
    if(positions[i]==n){ # no parents are allowed
      nodescores[i]<- 0 # there is only one score
    } else {
      #Nodes that are allowed to be its parents
      parentnodes<-permy[c((positions[i]+1):n)]
      # For each node i, calculate the scores
      node_scores <- numeric(length(parentnodes))  # Initialize node scores
      for (j in 1:length(parentnodes)) {
        # Want : z = log (exp a + exp b)
        a <- 0  # This is the '0' in the 'log(1 + exp(...))' expression
        b <- betas[parentnodes[j],i]
        #cat("beta :", parentnodes[j],i, "is:",betas[parentnodes[j],i], "\n")
        # For numerically stable computation: z = c + log(1 + exp (d-c))
        c <- max(a, b)
        d <- min(a, b)
        node_scores[j] <- c + log1p(exp(d - c))
      }
      nodescores[i] <- sum(node_scores)
    }
  }
  scores <- nodescores
  #scores$orderscores<-sum(nodescores)
  return(scores)
  
}

# currentpermy <- c(1, 6, 7, 5, 8, 2, 4, 3 )
# proposdedpermy <- c(7 ,6, 1, 5, 8, 2, 4, 3 )
# (currentpermy_score <- orderscore_betas (n,scorenodes = c(1,6,7), betas, currentpermy))
# (proposedpermy_score <- orderscore_betas (n,scorenodes = c(1,6,7) , betas, proposdedpermy))
#   
# beta : 6 1 is: -2.828407 
# beta : 7 1 is: -3.088243 
# beta : 8 1 is: -3.092051 
# beta : 7 6 is: 928.7908 
# 
# 
# beta : 8 1 is: -3.092051 
# beta : 1 6 is: -2.828407 
# beta : 6 7 is: 928.7908 
# beta : 1 7 is: -3.088243 

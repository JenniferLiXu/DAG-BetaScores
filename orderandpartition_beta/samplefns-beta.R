
# This function samples an element from a vector properly

propersample <- function(x){if(length(x)==1) x else sample(x,1)} 

# This function samples a single DAG score according to a list of possible betas

samplescore<-function(n, betas, permy){
  
  incidence<-matrix(numeric(n*n),nrow=n) # store the adjacency matrix
  sampledscore<- 0
  #Store the position of each element in the permutation (i.e. an inverse permutation)
  positions <- order(permy) 
  
  max_val <- 700  # Set a threshold to avoid overflow in exp
  
  for (child in 1:n){

    if(positions[child]==n){ # no parents are allowed
    } else {
      #Nodes that are allowed to be its parents
      parentnodes <- permy[c((positions[child]+1):n)]
      
      for (parent in parentnodes){
        if (is.na(betas[parent, child])){
          cat("For parent:", parent, child, " we have:",betas[parent, child], "\n")
          print(betas)
          print(permy)
        }

        incidence[parent, child] <- sample(c(0,1), 1, prob = c(1,exp(min(betas[parent, child], max_val))))
        sampledscore <- sampledscore + betas[parent, child]*incidence[parent, child]                      
      }
    }
  }
  DAG<-list()
  DAG$incidence<-incidence
  DAG$logscore<-sampledscore
  return(DAG)
}


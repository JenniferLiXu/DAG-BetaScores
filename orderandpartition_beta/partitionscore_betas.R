
# This function gives the sum of scores compatible with a given partition
# n
# scorenodes <- c(1:n)
# betas <- weighted_betas_proposed
# permy <- order_prev
# party <- partition

partitionscore<-function(n,scorenodes,betas,permy,party,posy = NULL){
  if(is.null(posy)){
    posy <- parttolist(n,party)
  }
  nodescores<- rep(0,n)
  m<-length(party)
  
  for (i in scorenodes){
    position<-which(permy==i)
    partyelement<-posy[position]
    if(partyelement==m){# no parents are allowed
      nodescores[i]<- 0 # there is only one score
    } else {
      parentnodes<-permy[which(posy > partyelement)] 

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

# This function takes in a partition and lists which partition 
# element each node belongs to

parttolist<-function(n,party){
  posy<-rep(0,n)
  kbot<-1
  m<-length(party)
  for (j in 1:m){
    ktop<-kbot+party[j]-1
    posy[kbot:ktop]<-j
    kbot<-ktop+1
  }
  return(posy)
}

# This function calculates permutations between neighbouring partition elements
# given the partition

parttopermneighbourposs<-function(n,party){
  m<-length(party)
  possibs<-rep(0,m-1)
  if(m>1){
    for (i in 1:(m-1)){
      possibs[i]<-party[i]*party[i+1]
    }	
  }
  return(possibs)
}

# This function calculates permutations excluding nodes in the same partition element
# given the partition

parttopermdiffelemposs<-function(n,party){
  m<-length(party)
  possibs<-rep(0,m-1)
  if(m>1){
    remainder<-n
    for (i in 1:(m-1)){
      remainder<-remainder-party[i]
      possibs[i]<-party[i]*remainder
    }	
  }
  return(possibs)
}

# This function returns the number of ways a partition can be joined or split
# note that adding an if statement for elemsize=1 seems to make the function slower

partysteps<-function(n,party){
  partypossibs<-rep(0,n)
  kbot<-1
  m<-length(party)
  for (j in 1:m){
    elemsize<-party[j]
    ktop<-kbot+elemsize-1
    partypossibs[kbot:ktop]<-choose(elemsize,1:elemsize)
    kbot<-ktop+1
  }
  return(partypossibs[1:(n-1)])
}

# This function returns the number of ways a node can be joined to another partition element

partyjoin<-function(n,party,posy){
  m<-length(party)
  joinpossibs<-rep(0,n)
  for(k in 1:n){
    joinpossibs[k]<-m-1
    nodeselement<-posy[k]
    if(party[nodeselement]==1){ #nodes in a partition element of size 1
      if(nodeselement<m){
        if(party[nodeselement+1]==1){ #and if the next partition element is also size 1
          joinpossibs[k]<-m-2 #we only allow them to jump to the left to count the swap only once
        }
      }
    }
  }
  return(joinpossibs)
}

# This function returns the number of ways a node can be move to a new partition element

partyhole<-function(n,party,posy){
  m<-length(party)
  holepossibs<-rep(0,n)
  for(k in 1:n){
    nodeselement<-posy[k]
    if(party[nodeselement]==1){ #nodes in a partition element of size 1 cannot move to the neighbouring holes
      holepossibs[k]<-m-1
      if(nodeselement<m){
        if(party[nodeselement+1]==1){ #and if the next partition element is also size 1
          holepossibs[k]<-m-2 #we only allow them to jump to the left to count the swap only once
        }
      }
    } else if(party[nodeselement]==2){ #nodes in a partition element of size 2 cannot move to the hole on the left
      holepossibs[k]<-m #since this would count the same splitting twice 
    } else {
      holepossibs[k]<-m+1
    }
  }
  return(holepossibs)
}

# this function takes in an adjacency matrix and returns the partition and permutation

DAGtopartition <- function(n,incidence){
  party<-c() # to store the partition
  permy <- numeric(n) # to store the permutation
  m <- n # counter
  while(m>0){
    topnodes<-which(colSums(incidence)==0) # find the outpoints
    incidence[topnodes,]<-0 # remove their edges
    incidence[cbind(topnodes,topnodes)]<-1 # add a one to their columns so they are no longer counted
    l<-length(topnodes) # the number of outpoints
    m<-m-l
    permy[(m+1):(m+l)]<-topnodes
    party<-c(l,party)
  }
  
  partypermy<-list()
  partypermy$party<-party
  partypermy$permy<-permy
  return(partypermy)
}


## Code for the manuscript ``Partition MCMC for inference on acyclic digraphs''
#
## Jack Kuipers and Giusi Moffa
## ETH
#
## Last modified: April 19, 2015
#
## Disclaimer: The code in this archive is not guaranteed to be optimised or free of bugs.
##        Please report any issues to the authors 
##        (jack.kuipers@bsse.ethz.ch and giusi.moffa@gmail.com).

# score a complete incidence matrix

DAGscore <- function(incidence, n){
	P_local <- numeric(n)   
   
	for (j in 1:n)  {
		parentnodes <- which(incidence[,j]==1)
		P_local[j]<-DAGcorescore(j,parentnodes,n)
	}

return(sum(P_local))
}

# score just certain nodes

DAGnodescore <- function(incidence, n, rescorenodes){
	P_local <- numeric(n)   
   
	for (j in rescorenodes)  {
		parentnodes <- which(incidence[,j]==1)
    P_local[j]<-DAGcorescore(j,parentnodes,n)
	}
		
return(P_local)
}

# Now we take in a matrix whose rows are the parent sets

TableDAGscore <- function(parentrows, j, n){
	nrows<-nrow(parentrows)
	P_local <- numeric(nrows)   
   
	for (i in 1:nrows)  {
		parentnodes <- parentrows[i,which(parentrows[i,]>0)]
		P_local[i]<-DAGcorescore(j,parentnodes,n)
	}

return(P_local)
}

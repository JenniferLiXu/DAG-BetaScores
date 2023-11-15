#Creating a small DAG model in R to 
#test the proposed method for obtaining beta values

#install.packages("bnlearn")
install.packages("BiDAG")
install.packages("gRain")

#Create a Small DAG
#library(bnlearn)
# Manually define a small DAG
#dag <- model2network("[A][B|A][C|A][D|B:C]")

seedset<-1 # set a seed?
seednumber<-101 # an example used

# Choose size of DAGs to consider
n<-4

# Choose maximum number of parents
maxparents<-3 # Maximum number of parents to allow


# load the necessary functions

source('./edgerevandstructure/structurefns.R')
source('./edgerevandstructure/structureMCMC.R')
#source('./edgerevandstructure/newedgerevfns.R')
#source('./edgerevandstructure/newedgerevmove.R')

source('./orderandpartition/samplefns.R')
source('./scoring/combinations.R')
source('./scoring/scorefns.R')
source('./scoring/scoretables.R') 

# load a simple score proportional to the number of edges in the DAG
source('./scoring/numedgescore.R')

#Use the function to calculate beta values
source('./calculateBetaScoresArray.R')

# Fill up a matrix with possible parents
parenttable<-listpossibleparents(maxparents,c(1:n))
tablelength<-nrow(parenttable[[1]]) # size of the table

scoretable<-scorepossibleparents(parenttable,n) 

# Define other parameters
iterations <- 20
moveprobs<-c(1)
startDAG<-matrix(numeric(n*n),nrow=n) # starting DAG is empty say
revallowed<-1 # allow standard edge reversals
stepsave<-1 #stepsave<-iterations/1000 #how often to save the result

if(seedset>0){
  set.seed(seednumber) # choose one?
}


example<-structureMCMC(n,startDAG,iterations,stepsave,maxparents,parenttable,scoretable,revallowed,moveprobs)
#Store all the DAGs
sampledDAGs<-example[[1]]
calculateBetaScoresArray(sampledDAGs, n, parenttable, scoretable)



#Use an existing MCMC method to sample new DAGs and update β values.
# This step will involve resampling DAGs, updating beta values, and using MCMC methods


# In this file, we are going to give a comparision between our method and order MCMC
.libPaths("/cluster/home/xuwli/R/x86_64-pc-linux-gnu-library/4.3") 
library(BiDAG) 
library(pcalg)
library(igraph)
library(bnlearn)
library(ggplot2)
library(combinat)

# load the necessary functions
#source('./edgerevandstructure/newedgerevfns.R')
#source('./edgerevandstructure/newedgerevmove.R')
source('./orderandpartition_beta/orderMCMC_betas.R')
source('./orderandpartition_beta/orderscore_betas.R')
source('./orderandpartition_beta/partitionMCMC_betas.R')
source('./orderandpartition_beta/partitionscore_betas.R')
source('./orderandpartition_beta/partitionmoves_betas.R')
source('./orderandpartition_beta/samplefns-beta.R')
source('./orderandpartition_beta/samplefns_party-beta.R')

source('./calculateBetaScoresArray.R')
source('./compareDAG_skeleton.R')
source('./BetaOrderSampler.R') #our method
source('./BetaOrderSampler_OneDAG.R') 
source('./BetaOrderSampler_gibbs.R') 







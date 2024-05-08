###################### DAG sets
# 10 DAGs in set 
# iteration 3000: 0.000659785 0.9957892 
# iteration 5000: 0.0006612062 0.9958414 

# 100 DAGs in set
# iteration 3000: 0.004165124 0.9493413 
# iteration 5000: 0.003184019 0.9667746 

###################### Gibbs sampling
num_iterations <- 2000
# weighted version         order_betaOrder123:  0.00448463 0.9778996 
# equally weighted version order_betaOrder300:  0.004134801 0.9740147 

num_iterations <- 3000
# weighted version         order_betaOrder123:  0.007491891 0.9452089 
# equally weighted version order_betaOrder300:  0.003967594 0.9751654 

num_iterations <- 5000
# weighted version         order_betaOrder123:  0.006823928 0.9670896 
# equally weighted version order_betaOrder300:  0.003714876 0.9757753

## With ess > 0.3 (10 sampled orders, each with 10 DAGs)
# iteration 3000:  0.002242364 0.9668648  burnin:0.5  acceptCount_after_burnin:131
# iteration 4000:  0.001860352 0.9717618   burnin:0.5  acceptCount_after_burnin:290

# With ess > 0.3 (10 sampled orders, each with 10 DAGs) ---- modified beta
# iteration 3000: 0.004157752 0.9466148 
# iteration 5000: 0.004156788 0.9488965 

###################### End of Gibbs sampling
# > pedges_comp
# [[1]]
# crim          zn      indus        chas
# crim  0.000000000 0.001248439 0.38202247 0.538077403
# zn    0.001248439 0.000000000 0.19475655 0.007490637
# indus 0.617977528 0.805243446 0.00000000 0.063670412
# chas  0.136079900 0.001248439 0.01248439 0.000000000

wB <- sum(is_results$log_diff)/nr_sample
inv_wA <- sum((oDAGnBeta_logscore_list - BiDAGscore_prev_list))/nr_sample
compress_DAG[[i+1]] <- is_results$compress_dag
# [[2]]
# crim          zn      indus        chas
# crim  0.0000000000 0.001320000 0.34600000 0.521212726
# zn    0.0000822884 0.000000000 0.41920000 0.008688032
# indus 0.6540000000 0.580800000 0.00000000 0.042286833
# chas  0.1346322698 0.004206033 0.02186849 0.000000000

compress_DAG[[i+1]] <- Reduce("+", lapply(1:length(incidence_matrices), function(k) incidence_matrices[[k]]/nr_sample))
# [[2]]
# crim      zn   indus    chas
# crim  0.00000 0.00132 0.34600 0.51978
# zn    0.00008 0.00000 0.41920 0.00990
# indus 0.65400 0.58080 0.00000 0.05522
# chas  0.13516 0.00412 0.03588 0.00000

wB <- sum(is_results$log_diff*weights_proposed)
inv_wA <- sum((oDAGnBeta_logscore_list - BiDAGscore_prev_list)*weights_prev)
compress_DAG[[i+1]] <- is_results$compress_dag
# [[2]]
# crim          zn      indus        chas
# crim  0.000000000 0.001879281 0.23280000 0.489124754
# zn    0.001678616 0.000000000 0.05880000 0.008113045
# indus 0.767200000 0.941200000 0.00000000 0.075329203
# chas  0.178180226 0.002812892 0.02384516 0.000000000
results_seed100$betas
# [,1]      [,2]      [,3]       [,4]
# [1,]  0.0000000 -6.290409 42.525494  0.5029098
# [2,] -6.3312354  0.000000 16.197000 -4.0258435
# [3,] 42.0381171 16.197000  0.000000 -1.2932571
# [4,] -0.9591683 -5.038134 -1.293257  0.0000000


wB <- sum(is_results$log_diff*weights_proposed)
inv_wA <- sum((oDAGnBeta_logscore_list - BiDAGscore_prev_list)*weights_prev)
compress_DAG[[i+1]] <- Reduce("+", lapply(1:length(incidence_matrices), function(k) incidence_matrices[[k]]/nr_sample))
# [[2]]
# crim      zn   indus    chas
# crim  0.00000 0.04470 0.23280 0.49600
# zn    0.12220 0.00000 0.05880 0.00930
# indus 0.76720 0.94120 0.00000 0.07516
# chas  0.16414 0.00414 0.04244 0.00000

wB <- sum(is_results$log_diff*weights_proposed)
inv_wA <- sum((oDAGnBeta_logscore_list - BiDAGscore_prev_list)*weights_for_prevD)
compress_DAG[[i+1]] <- is_results$compress_dag
# [[2]]
# crim          zn      indus       chas
# crim  0.000000000 0.001816183 0.18580000 0.47896311
# zn    0.002069068 0.000000000 0.05180000 0.00748412
# indus 0.814200000 0.948200000 0.00000000 0.07932004
# chas  0.176958202 0.006009231 0.02318009 0.00000000
results_seed100$betas
# [,1]      [,2]      [,3]       [,4]
# [1,]  0.0000000 -6.290409 42.525494  0.5029098
# [2,] -6.3312354  0.000000 16.197000 -4.0258435
# [3,] 42.0381171 16.197000  0.000000 -1.2932571
# [4,] -0.9591683 -5.038134 -1.293257  0.0000000


results_seed123$betas
# [,1]      [,2]      [,3]      [,4]
# [1,]  0.0000000  3.320621 32.914464  1.477611
# [2,] -6.4945315  0.000000 16.197000 -3.350996
# [3,] 40.0887149 16.197000  0.000000 -1.293257
# [4,] -0.9591683 -3.350996 -2.980395  0.000000
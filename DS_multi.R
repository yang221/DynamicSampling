rm(list = ls())

library(doMPI)
library(MASS)
library(GA)

# source all function files in the folder "functions"
function_files = list.files(path = 'functions_2', pattern="*.R")
for(f in function_files){
  source(paste('functions_2/', f, sep=''))
}

#Conduct a case study for dynamic sampling design
load("sampleData/tool_wear_data_36_45.RData")


## Presettings
step1 = 6 # Level 1 stepsize
step2 = 3 # Level 2 stepsize

grid_mat <- grid_seg(tool_wear_data$no.row, tool_wear_data$no.col, step2) # Grid segmentation
no_grid <- max(max(grid_mat))

lambda <- 1.3 # Tunning parameter in the cost function
# Higher lambda --> better precision & higher cost
# Lower lambda --> poorer precision & lower cost

c1 <- 1 # Unit cost for Level 1
c2 <- 5 # Unit cost for Level 2


W <- tool_wear_data$W

sigma1sq <- 1
sigma2sq <- 0.1

# confidence_level = 0.05
confidence_level <- 0.05

D1 <- level1_design(tool_wear_data$no.row, tool_wear_data$no.col, step1) # Level 1 design is fixed
# Generate M1 first
no_sampling <- tool_wear_data$no.sampling
M1 <- array(list(), no_sampling)
for (i in 1:no_sampling){
  M1[[i]] = measurement(D1, tool_wear_data$locs.index, unlist(tool_wear_data$Z[i]), 1, sigma1sq)
}

## Design
design_results = array(list(), no_sampling)

Phi = array(list(), no_sampling)
A = array(list(), no_sampling)
B = array(list(), no_sampling)
Atrue = array(list(), no_sampling)
Btrue = array(list(), no_sampling)
C = array(list(), no_sampling)
D = array(list(), no_sampling)
m1 = array(list(), no_sampling)
m2 = array(list(), no_sampling)
coords1 = array(list(), no_sampling)
coords2 = array(list(), no_sampling)
precision = array(list(), no_sampling)
Sigma_e = array(list(), no_sampling)
cov_eta = array(list(), no_sampling)
Zhat = array(list(), no_sampling)
Zhat_mean = array(NaN, no_sampling)
h_decision = array(NaN, no_sampling)
D1_final = array(list(), no_sampling)
D2 = array(list(), no_sampling)
# M1 = array(list(), no_sampling)
M2 = array(list(), no_sampling)
precision_true = array(list(), no_sampling)
GA_result = array(list(), no_sampling)

decision1 = array(NaN, no_sampling)
decision1[1] = 1
decision2 = array(NaN, no_sampling)
decision2[1] = 1

plot_on = 1

l = 1
n = tool_wear_data$n
MU = array(0, n)
cov_mat = matrix(NaN, nrow = n, ncol = n)
for(i in 1:n){
  for(j in 1:n){
    x1 = tool_wear_data$coords.all[i,]
    x2 = tool_wear_data$coords.all[j,]
    cov_mat[i, j] <- gp_exp_kernel(x1, x2, l)
  }
}
phi_hat = array(NaN, no_sampling)
phi_hat[1] = tool_wear_data$phi[1]

A0 = 0.05*matrix(1, nrow = n, ncol = n) + 0.95*diag(1, n, n)

Zhat[[1]] = phi_hat[1]*tool_wear_data$Z0

for(i in 1:no_sampling){
  Phi[[i]] = simplify2array(tool_wear_data$Phi[[i]])[,,]
}
Phi_hat = Phi

# create a cluster
np <- mpi.universe.size() - 1
cl <- startMPIcluster(np)
registerDoMPI(cl)

for(i in 1:no_sampling){
  current_time = i
  
  cov_eta[[i]] = cov_mat
  
  # Test whether phi_t changes, phi_t needs to be estimated first in
  # order to estimate Bt, based on Level 1 measurements
  if (i > 1){  # When i == 1, phi_hat(1) = phi(1), no need to estimate
    phi_prev = phi_hat[current_time-1]
    #     [Yprev, Ynow, mutural_number1, mutural_number2] = ...
    #         Y_calculation(M1, M2, current_time); % Old version
    Y_temp_result = Y_calculation_level1(M1, current_time)
    Yprev = Y_temp_result[1]
    Ynow = Y_temp_result[2]
    mutural_number1 = nrow(M1[[1]]) # Fixed number of Level 1 measurement points
    decision = phi_hypothesis_level1(Yprev, Ynow, mutural_number1, phi_prev, sigma1sq, confidence_level)
    decision1[i] = decision
    phi_hat[i] = phi_update(decision, current_time, Yprev, Ynow, phi_prev)
  } 
  
  if (i == 1){ # First design involves A0
    B[[1]] = Phi[[1]] %*% A0 %*% t((Phi[[1]])) + cov_mat
    Btrue[[1]] = Phi[[1]] %*% A0 %*% t((Phi[[1]])) + cov_mat
  } else{ # Update B{i} using A{i-1} and phi_hat{i}
    Phi_hat[[i]] = diag(1, n, n)*phi_hat[i] # update phi_hat
    B[[i]] = Phi_hat[[i]] %*% A[[i-1]] %*% t(Phi_hat[[i]]) + cov_eta[[i]]
    Btrue[[i]] = Phi[[i]] %*% Atrue[[i-1]] %*% t((Phi[[i]])) + cov_eta[[i]] # true B matrix, used to calculate Vt
  }
  Bnow = B[[i]]
  Bnow_true = Btrue[[i]]
  #save.image(file= paste("presetting/presetting-",i,".RData", sep=''))

  print(paste('Dynamic sampling stage',current_time,'starts.'))
  # call multi_design function which in turn call ga_mpi function to do GA in parallel
  multi_design_result <- multi_design(plot_on, current_time, tool_wear_data$Z, tool_wear_data$locs.index, no_grid, grid_mat, D1, c1, c2, sigma1sq, sigma2sq, Bnow, Bnow_true, W, lambda, tool_wear_data$coords.all)

  design_results[[i]] = multi_design_result
  
  m1[[i]] = multi_design_result$m1
  m2[[i]] = multi_design_result$m2
  coords1[[i]] = multi_design_result$coords1
  coords2[[i]] = multi_design_result$coords2
  D[[i]] = multi_design_result$D
  C[[i]] = multi_design_result$C
  Sigma_e[[i]] = multi_design_result$Sigma_e
  precision[[i]] = multi_design_result$precision
  D1_final[[i]] = multi_design_result$D1
  D2[[i]] = multi_design_result$D2
  #     M1[[i]] = multi_design_result$M1
  M2[[i]] = multi_design_result$M2
  precision_true[[i]] = multi_design_result$precision_true
  GA_result[[i]] = multi_design_result$GA

  print(paste('Number of iterations: ',multi_design_result$GA@iter,'.'))
  print(paste('Value of optimal cost: ',multi_design_result$GA@fitnessValue,'.')) 
  
  # Estimated A
  A[[i]] = B[[i]] - B[[i]] %*% t(D[[i]]) %*% ginv(D[[i]] %*% B[[i]] %*% t(D[[i]]) + Sigma_e[[i]]) %*% D[[i]] %*% B[[i]]
  # True A
  Atrue[[i]] = Btrue[[i]] - Btrue[[i]] %*% t(D[[i]]) %*% ginv(D[[i]] %*% Btrue[[i]] %*% t(D[[i]]) + Sigma_e[[i]]) %*% D[[i]] %*% Btrue[[i]]
  
  if(i == 1){
    Zhat[[1]] = phi_hat[1]*tool_wear_data$Z0 + unlist(tool_wear_data$eta[[1]])
  }else{
    # Test whether phi_t changes, AGAIN, using both Levels 1 & 2 data and update its value
    phi_prev = phi_hat[current_time-1]
    Y_calculation_result = Y_calculation(M1, M2, current_time)
    Yprev = Y_calculation_result[1]
    Ynow = Y_calculation_result[2]
    mutural_number1 = Y_calculation_result[3]
    mutural_number2 = Y_calculation_result[4]
    
    decision = phi_hypothesis(Yprev, Ynow, mutural_number1,mutural_number2, phi_prev, sigma1sq, sigma2sq, confidence_level)
    decision2[i] = decision
    phi_hat[i] = phi_update(decision, current_time, Yprev, Ynow, phi_prev)
    
    # Predict Z at time i (height)
    #         Zhat{i} = phi_hat(i)*Zhat{i-1} + eta0;
    Zhat[[i]] = phi_hat[i]*Zhat[[i-1]] + unlist(tool_wear_data$eta[[i]])
  }
  
  save.image(file = paste("results/dynamic_sampling_tool_wear_result_temp.RData", sep=''))
  
  print(paste('Dynamic sampling stage',current_time,'has been completed. '))
}

save.image(file = "results/dynamic_sampling_tool_wear_result.RData")

closeCluster(cl)
mpi.finalize()
q()

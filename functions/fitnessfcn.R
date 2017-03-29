# Calculate the objective function for GA optimization
# Output is the value of the objective function
# Inputs include:
  #   d: the design vector (n_g by 1), where 0 means no measurement, 1 means
#   Level 2 measurement for the grid
#   c1 and c2 are the unit costs for Levels 1 and 2, respectively

fitnessfcn <- function(d, grid_mat, D1, c1, c2, sigma1sq, sigma2sq, Bnow, W, lambda){
  
  # First convert d to D2
  D2 = d_to_D2(d, grid_mat)
  
  # For certain locations where both D1 and D2 measure, change Sigma_e to
  # Level 2 repeatability, but the cost remains
  
  m1 = sum(D1)
  m2 = sum(D2)
  
  cost = m1*c1 + m2*c2
  if(m2>0){
    D1_locs = matrix_locs(D1)
    D2_locs = matrix_locs(D2)
    overlap_locs = which(D1_locs %in% D2_locs)
    if(length(overlap_locs) > 0){ # Remove locations from D1 where D2 also measures
      D1 = D1[-overlap_locs,]
    }
  }
  
  m1_new = sum(sum(D1))
  
  # Assign values for Sigma_e (measurement error variance matrix)
  Sigma_e = matrix(0, nrow = m1_new+m2, ncol = m1_new+m2)
  if(m1_new>0){
    Sigma_e[1:m1_new,1:m1_new] = sigma1sq*diag(m1_new)
    Sigma_e[-1:-m1_new, -1:-m1_new] = sigma2sq*diag(m2)
  }else{
    Sigma_e=sigma2sq*diag(m2)
  }
  
  
  D = rbind(D1, D2)
  D_locs = matrix_locs(D)
  
  A = Bnow - Bnow[,D_locs] %*% ginv(Bnow[D_locs, D_locs] + Sigma_e) %*% Bnow[D_locs,]
  
  precision = sum(diag(A %*% W))
  
  obj_value = cost + lambda*precision
  return(-obj_value)
}

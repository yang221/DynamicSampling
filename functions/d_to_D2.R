# Coding from d to K2
# Construct the design for Level 2 from a 0-1 vector, d
# grid_mat is a matrix whose elements are the grid incides
# Each element in d is corresponding to a square region
# grid_mat is a matrix, the element values represent which grid it belongs to

d_to_D2 <- function(d, grid_mat){
  n = length(grid_mat) # Total number of locations
  
  grid_vec = c(t(grid_mat))
  
  sel_grid = which(d==1)
  # locations where need Level 2 measurements
  level2_index = which(grid_vec %in% sel_grid == 1)
  m2 = length(level2_index)
  if(m2 > 0){
    D2 = matrix(0, nrow = m2, ncol = n)
    for(i in 1:m2){
      D2[i, level2_index[i]] = 1 
    }
    return(D2) 
  }
  else{
    return(c())
  }
}
# Generate the measurement data from the design matrix, D
# M is a m by 4 matrix:
  # The first two columns of M are the coordinates, the third column is the
# height, and the fouth column is the Level (1 or 2)
# Calculate the coordinates from a design matrix D
measurement <- function(D, locs_index, Z_true, level, sigmasq){
  if(length(D)>0){
    m <- nrow(D)
    measurement_error <- rnorm(n = m, mean = 0, sd = sqrt(sigmasq))
    M <- matrix(nrow = m, ncol = 4)
    M[,4] <- level
    for(i in 1:m){
      M[i, 1:2] <- locs_index[[which(D[i,]==1)]][[1]]
      M[i, 3] <- Z_true[which(D[i,]==1)] + measurement_error[i]
    }
    return(M)
  }else{
    return(c())
  }
}
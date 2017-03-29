# Generate K1 matrix with a prespecified stepsize

level1_design <- function(no_row, no_col, step){
  n <- no_row*no_col
  m1 <- ceiling(no_row/step)*ceiling(no_col/step)
  K1 <- matrix(0, nrow = m1, ncol = n)
  temp <- grid_seg(no_row, no_col, 1)
  flag <- 0
  
  for (i in 1:no_row){
    for (j in 1:no_col){
      if (i%%step == 1 && j%%step == 1){
        flag <- flag + 1
        K1[flag, temp[i,j]] <- 1
      }
    }
  }
  return(K1)
}

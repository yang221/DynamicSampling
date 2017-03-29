# Segment a surface into equal square grids
# Assume the surface takes a rectangular shape
# The inputs include:
#   no_row and no_col are the size of the surface;
#   step is the stepsize (lateral resolution) for Level 2
# The output is a matrix with equal size, and the elements are the grid number

# Please note that row is correponding to y, column is corresponding to x

grid_seg <- function(no_row, no_col, step){
  no_grid_col <- ceiling(no_col/step)
  grid_mat <- matrix(, nrow = no_row, ncol = no_col)
  for(y in 1:no_row){
    for(x in 1:no_col){
      grid_mat[y,x] <- floor((y-1)/step)*no_grid_col + ceiling(x/step)
    }
  }
  return(grid_mat)
}

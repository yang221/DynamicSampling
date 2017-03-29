# Calculate the coordinates from a design matrix D
# coords is a m by 2 matrix, each row of which corresponds to one measurement point

sampling_locs <- function(D, locs_index){
  m = nrow(D)
  if(length(D)>0){
    coords = matrix(NA, nrow = m, ncol = 2)
    for(i in 1:m){
      coords[i,] = unlist(locs_index[which(D[i,] == 1)])
    }
    return(coords)
  }else{
    return(c())
  }

}
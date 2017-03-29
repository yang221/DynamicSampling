matrix_locs <- function(D){
  len = nrow(D)
  result = array(0,len)
  for(i in 1:len){
    result[i] = which(D[i,]==1)
  }
  return(result)
}

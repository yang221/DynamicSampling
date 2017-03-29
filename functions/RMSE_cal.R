RMSE_cal <- function(Z,Zhat,n){
  no_sampling = length(Z)
  RMSE = array(NaN, no_sampling)
  for (i in 1:no_sampling){
    error = unlist(Z[[i]]) - Zhat[[i]]
    RMSE[i] = norm(error, type = '2')/sqrt(n)
  }
  return(RMSE)
}

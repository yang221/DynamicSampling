Y_calculation_level1 <- function(M1, current_time){
  data_prev = M1[[current_time-1]]
  data_now = M1[[current_time]]
  Yprev = mean(data_prev[, 3])
  Ynow = mean(data_now[, 3])
  return(c(Yprev, Ynow))
}


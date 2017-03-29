# Update the value of phi
phi_update <- function(decision, current_time, Yprev, Ynow, phi_prev){
  h = Ynow/Yprev
  if (decision == 0){ # shift occurs
    phi_hat = h
  } else{ # no change
    # update using historical data
    phi_hat = ((current_time-1)*phi_prev + h)/current_time
  } 
  return(phi_hat)
}
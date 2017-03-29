Y_calculation <- function(M1, M2, current_time){
  coords1_prev = M1[[current_time-1]][,1:2]
  coords2_prev = M2[[current_time-1]][,1:2]
  coords1_now = M1[[current_time]][,1:2]
  coords2_now = M2[[current_time]][,1:2]
  
  mutural1_prev = ismemberbyrow(coords1_prev, coords1_now)
  mutural1_now = ismemberbyrow(coords1_now, coords1_prev)
  mutural2_prev = ismemberbyrow(coords2_prev, coords2_now)
  mutural2_now = ismemberbyrow(coords2_now, coords2_prev)
  
  M1_mutural_prev = M1[[current_time-1]][mutural1_prev,]
  M1_mutural_now = M1[[current_time]][mutural1_now,]
  M2_mutural_prev = M2[[current_time-1]][mutural2_prev,]
  M2_mutural_now = M2[[current_time]][mutural2_now,]
  M_mutural_prev = rbind(M1_mutural_prev, M2_mutural_prev)
  M_mutural_now = rbind(M1_mutural_now, M2_mutural_now)
  
  
  mutural_number1 = length(mutural1_prev)
  mutural_number2 = length(mutural2_prev)
  Yprev = mean(M_mutural_prev[, 3])
  Ynow = mean(M_mutural_now[, 3])
  
  return(c(Yprev, Ynow, mutural_number1, mutural_number2))
}


level2_loc_freq <- function(coords2, coords_all){
  no_sampling = length(coords2)
  no_locs = sqrt(nrow(coords_all))
  freq_mat_level2 = matrix(0, no_locs, no_locs)
  
  for(i in 1:no_sampling){
    temp_coords2 = coords2[[i]]
    for(j in 1:nrow(temp_coords2)){
      temp_coord = temp_coords2[j,]
      freq_mat_level2[temp_coord[2],temp_coord[1]] = freq_mat_level2[temp_coord[2],temp_coord[1]] + 1
    }
  }
  return(freq_mat_level2)
}
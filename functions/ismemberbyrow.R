ismemberbyrow <- function(D1, D2){
  overlap_locs = c()
  if(length(D1)>0&&length(D2)>0){
    for(i in 1:nrow(D1)){ # Remove locations from D1 where D2 also measures
      for(j in 1:nrow(D2)){
        if(isTRUE(all.equal(D1[i,],D2[j,]))){
          overlap_locs = c(overlap_locs,i)
        }
      }
    }
  }
  return(overlap_locs)
}

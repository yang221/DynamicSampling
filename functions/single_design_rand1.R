single_design_rand1 <- function(plot_on, suffix, Z, locs_index, no_grid, grid_mat, D1, D2, sigma1sq, sigma2sq, Bnow, Bnow_true, W, coords_all){
  # Perform random design strategy
  # Level 1 locations are fixed, randomly select a certain number of Level 2 rectangles

  m1 = nrow(D1)
  m2 = nrow(D2)
  coords1 = sampling_locs(D1, locs_index)
  
  grid_mat_available = grid_mat
  # for i = 1:m1
  #     temp_removal_coord = coords1(i,:);
  #     removal_grid_index = grid_mat(temp_removal_coord(2),temp_removal_coord(1));
  #     grid_mat_available(grid_mat_available == removal_grid_index) = NaN;
  # end
  grid_available_index = unique(c(grid_mat_available))
  grid_available_index = grid_available_index[which(is.nan(grid_available_index) == 0)]
  m2_grid = round(m2/grid_size)
  level2_grid = sort(sample(grid_available_index,m2_grid))
  no_grid = max(max(grid_mat))
  d = array(0, no_grid)
  d[level2_grid] = 1
  D2 = d_to_D2(d, grid_mat) # D2 for random sampling
  coords2 = sampling_locs(D2, locs_index)
  
  # Deal with partial grids
  m2_rand = nrow(D2)
  Sigma_e_rand = matrix(0, nrow = m1+m2_rand, ncol = m1+m2_rand)
  Sigma_e_rand[1:m1,1:m1] = sigma1sq*diag(m1)
  Sigma_e_rand[-1:-m1, -1:-m1] = sigma2sq*diag(m2_rand)
  D = rbind(D1, D2)
  
  A = Bnow - Bnow %*% t(D) %*% ginv(D %*% Bnow %*% t(D) + Sigma_e_rand) %*% D %*% Bnow
  A_true = Bnow_true - Bnow_true %*% t(D) %*% ginv(D %*% Bnow_true %*% t(D) + Sigma_e_rand) %*% D %*% Bnow_true
  precision = sum(diag(A %*% W))
  precision_true = sum(diag(A_true %*% W)) # true precision calculated based on true Phi
  
  if(plot_on == 1){
    #png(filename = paste("figures/design_time_",suffix,"_dynamic.png", sep=''), height = 480, width = 620, res = 100)
    pdf(file = paste("figures/design_time_",suffix,"_random.pdf", sep=''),
        height = 5.2, width = 7)
    par(xpd = TRUE, mar = c(4.1, 4.1, 3.1, 8.1))
    plot(coords_all[,1], coords_all[,2], pch = 20, cex = 0.5, 
         xlab="x-axis", ylab="y-axis", xlim = c(-0.5, ncol(grid_mat)+1.25), asp = 1)
    points(coords1[,1], coords1[,2], pch =4, col = 2)
    points(coords2[,1], coords2[,2], pch = 1 , col = 4)
    text(x= -0.5, y = seq(nrow(grid_mat)), labels = seq(1, length(grid_mat)-ncol(grid_mat)+1, ncol(grid_mat)), cex = 0.6)
    text(x= ncol(grid_mat) + 1.25, y = seq(nrow(grid_mat)), labels = seq(ncol(grid_mat), length(grid_mat), ncol(grid_mat)), cex = 0.6)
    legend('right', inset = -0.3, legend = c('All locations', 'Level1', 'Level 2'), 
           pch = c(20, 4, 1), col = c(1, 2, 4), cex = 0.8, bty = 'n')
    dev.off()
  }
  
  single_design_result <- list(m1 = nrow(D1), m2 = nrow(D2), coords1 = sampling_locs(D1, locs_index), 
                               coords2 = sampling_locs(D2, locs_index), D = D, C = C, Sigma_e_rand = Sigma_e_rand,
                               precision = precision, D1 = D1, D2 = D2, 
                               M1 = measurement(D1, locs_index, unlist(Z[current_time]), 1, sigma1sq),
                               M2 = measurement(D2, locs_index, unlist(Z[current_time]), 2, sigma2sq),
                               precision_true = precision_true)
  return(single_design_result)
}



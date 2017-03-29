single_design <- function(plot_on, suffix, Z, locs_index, no_grid, grid_mat, D1, c1, c2, sigma1sq, sigma2sq, Bnow, Bnow_true, W, lambda, coords_all){
  # Perform one-time design at time t
  
  # Given a design of Level 1, choose the best Level 2 design
  GA <- ga('binary',fitnessfcn,grid_mat, D1, c1, c2, sigma1sq, sigma2sq, Bnow, W, lambda, 
           nBits = no_grid, popSize = 100, maxiter = 2000, run = 30, parallel = TRUE)
  d = GA@solution
  C = -GA@fitnessValue
  D2 = d_to_D2(d, grid_mat)
  
  # Disable to ensure there are always mutually measured locations
  # overlap_locs = ismember(D1,D2,'rows');
  # D1(overlap_locs, :) = []; % Remove locations from D1 where D2 also measures
  
  m1_new = nrow(D1)
  m2 = nrow(D2)
  
  if(length(D2)<=0){
    m2 = 0
  }
  
  coords1 = sampling_locs(D1, locs_index)
  coords2 = sampling_locs(D2, locs_index)
  
  Sigma_e = matrix(0, nrow = m1_new+m2, ncol = m1_new+m2)
  Sigma_e[1:m1_new,1:m1_new] = sigma1sq*diag(m1_new)
  Sigma_e[-1:-m1_new, -1:-m1_new] = sigma2sq*diag(m2)
  
  D = rbind(D1, D2)
  A = Bnow - Bnow %*% t(D) %*% ginv(D %*% Bnow %*% t(D) + Sigma_e) %*% D %*% Bnow
  A_true = Bnow_true - Bnow_true %*% t(D) %*% ginv(D %*% Bnow_true %*% t(D) + Sigma_e) %*% D %*% Bnow_true
  precision = sum(diag(A %*% W))
  precision_true = sum(diag(A_true %*% W)) # true precision calculated based on true Phi
  
  if(plot_on == 1){
    #png(filename = paste("figures/design_time_",suffix,"_dynamic.png", sep=''), height = 480, width = 620, res = 100)
    pdf(file = paste("figures/design_time_",suffix,"_dynamic.pdf", sep=''),
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
                               coords2 = sampling_locs(D2, locs_index), D = D, C = C, Sigma_e = Sigma_e,
                               precision = precision, D1 = D1, D2 = D2, 
                               M1 = measurement(D1, locs_index, unlist(Z[current_time]), 1, sigma1sq),
                               M2 = measurement(D2, locs_index, unlist(Z[current_time]), 2, sigma2sq),
                               precision_true = precision_true, GA = GA)
  return(single_design_result)
}



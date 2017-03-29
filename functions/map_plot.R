map_plot <- function(mat, coords_all, clim, fig_name){
  pdf(file = paste("figures/", fig_name, ".pdf", sep=''),
      height = 5.2, width = 6)
  image.plot(1:ncol(mat), 1:nrow(mat), mat, col = clim,  xlim = c(-1.25, ncol(mat)+2), 
        asp = 1, bty = 'n', xlab = 'x', ylab = 'y', xpd = TRUE)
  points(coords_all[,1], coords_all[,2], pch = 20, cex = 0.25, col = 'cadetblue')
  text(x= -0.5, y = seq(nrow(grid_mat)), labels = seq(1, length(grid_mat)-ncol(grid_mat)+1, ncol(grid_mat)), cex = 0.6)
  text(x= ncol(grid_mat) + 1.5, y = seq(nrow(grid_mat)), labels = seq(ncol(grid_mat), length(grid_mat), ncol(grid_mat)), cex = 0.6)
  dev.off()
}
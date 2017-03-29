# Calculate the covariance between two locations, x1 and x2, which are 2-d coordinates
# based on Squared Exponential Covariance Function
# l is the length-scale parameter
gp_exp_kernel <- function(x1, x2, l){
  r = norm(x1 - x2, type = '2')
  cov_gp = exp(-r^2/2/l^2)
  return(cov_gp)
}

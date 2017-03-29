# Test whether the value of phi changes
phi_hypothesis <- function(Yprev, Ynow, mutural_number1, mutural_number2, 
               phi_prev, sigma1sq, sigma2sq, alpha){
  Ynow_hat = phi_prev*Yprev
  sigmasq_prev_now = sigma1sq/mutural_number1 + sigma2sq/mutural_number2
  phi_statistic = abs((Ynow-Ynow_hat)/sqrt((1+phi_prev^2)*sigmasq_prev_now))
  decision = phi_statistic < qnorm(1-alpha) # 0 = reject H0 --> shift occurs
}


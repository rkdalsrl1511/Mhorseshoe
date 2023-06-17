
# rejection sampler function
rejection_sampler <- function(Epsilon, a, b){

  # p-dimensions
  p <- length(Epsilon)
  eta <- rep(1, p)

  # Calculate gradient lambda (0 < a < 1 < b)
  lambda2 <- lambda_calculation(Epsilon, a, 1)
  lambda3 <- lambda_calculation(Epsilon, 1, b)

  # Sampling Probability Calculation
  v <- v_calculation(Epsilon, lambda2, lambda3, a, b)
  total_v <- apply(v, MARGIN = 1, sum)
  prob_v <- v/total_v # p*4 matrix

  # sampling large epsilon eta(ordinary rejection sampler)
  large_eps_idx <- which(Epsilon > 1)
  eta <- sample_large_eps_eta(eta, large_eps_idx, prob_v, Epsilon)

  # sampling small epsilon eta(Johndraw's rejection sampler)
  small_eps_idx <- which(Epsilon <= 1)
  eta <- sample_small_eps_eta(eta, small_eps_idx, prob_v, Epsilon, lambda2, lambda3, a, b)

  return(eta)

}

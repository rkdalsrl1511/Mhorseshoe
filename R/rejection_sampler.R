
# rejection sampler function
rejection_sampler <- function(Epsilon, a = 1/5, b = 10){

  # p-dimensions
  p <- length(Epsilon)
  Epsilon <- ifelse(Epsilon == 0, 10^(-4),
                    ifelse(Epsilon >= 10^12, 10^12, Epsilon))
  eta <- rep(-1, p)
  rejected_index <- 1:p

  # Calculate gradient lambda (0 < a < 1 < b)
  lambda2 <- lambda_calculation(Epsilon, a, 1)
  lambda3 <- lambda_calculation(Epsilon, 1, b)

  # Sampling Probability Calculation
  v <- v_calculation(Epsilon, lambda2, lambda3, a, b)
  total_v <- apply(v, MARGIN = 1, sum)
  prob_v <- v/total_v # p*4 matrix

  # Sampling Local Shrinkage Parameters
  eta <- sample_eta(eta, rejected_index, prob_v, Epsilon, lambda2, lambda3, a, b)

  # return
  eta

}

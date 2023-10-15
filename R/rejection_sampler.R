## rejection sampler to sample local shrinkage parameters
#' @keywords internal
rejection_sampler <- function(Epsilon, a, b) {
  p <- length(Epsilon)
  eta <- rep(1, p)
  lambda2 <- (f_e_x(Epsilon, 1/Epsilon) - f_e_x(Epsilon, a/Epsilon)) /
    (1/Epsilon - a/Epsilon)
  lambda3 <- (f_e_x(Epsilon, b/Epsilon) - f_e_x(Epsilon, 1/Epsilon)) /
    (b/Epsilon - 1/Epsilon)
  A <- f_e_x(Epsilon, a/Epsilon)
  B <- f_e_x(Epsilon, b/Epsilon)
  I <- f_e_x(Epsilon, 1/Epsilon)
  v1 <- log(1+a/Epsilon)
  v2 <- (1/lambda2) * exp(-A) * (1 - exp(-(I-A)))
  v3 <- (1/lambda3) * exp(-I) * (1 - exp(-(B-I)))
  v4 <- (1/Epsilon) * exp(-B)
  v <- cbind(v1, v2, v3, v4)
  total_v <- apply(v, MARGIN = 1, sum)
  prob_v <- v/total_v
  large_eps_idx <- which(Epsilon > 1)
  eta <- sample_large_eps_eta(eta, large_eps_idx, prob_v, Epsilon)
  small_eps_idx <- which(Epsilon <= 1)
  eta <- sample_small_eps_eta(eta, small_eps_idx, prob_v, Epsilon, lambda2,
                              lambda3, a, b)
  return(eta)
}

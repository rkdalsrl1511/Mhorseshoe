## Eta sampling function for small epsilon in rejection sampler
#' @keywords internal
sample_small_eps_eta <- function(eta, idx, prob_v, Epsilon, lambda2, lambda3,
                                 a, b) {
  p <- length(idx)
  while(p != 0){
    u <- stats::runif(p, min = 0, max = 1)
    selected_interval <- ifelse(u < prob_v[idx, 1], 1,
                                ifelse(u < prob_v[idx, 1] + prob_v[idx, 2], 2,
                                       ifelse(u < prob_v[idx, 1] +
                                                prob_v[idx, 2] +
                                                prob_v[idx, 3], 3, 4)))
    u_z <- stats::runif(p, min = 0, max = 1)
    z <- sample_z(selected_interval, Epsilon[idx], lambda2[idx], lambda3[idx],
                  u_z, a, b)
    f_l <- f_L_z(selected_interval, Epsilon[idx], z, a, b, lambda2[idx],
                 lambda3[idx])
    f <- f_e_x(Epsilon[idx], z)
    u <- stats::runif(p, min = 0, max = 1)
    eta[idx] <- ifelse(u < exp(-(f- f_l)), z, -1)
    idx <- which(eta == -1)
    p <- length(idx)
  }
  return(eta)
}

## Eta sampling function for large epsilon in rejection sampler
#' @keywords internal
sample_large_eps_eta <- function(eta, idx, prob_v, Epsilon) {
  p <- length(idx)
  while(p != 0){
    v <- stats::runif(p, min = 0, max = 1)
    z <- - log(1-v) / Epsilon[idx]
    u <- stats::runif(p, min = 0, max = 1)
    eta[idx] <- ifelse(u < 1/(1+z), z, -1)
    idx <- which(eta == -1)
    p <- length(idx)
  }
  return(eta)
}

## Internal function in rejection sampler
#' @keywords internal
f_e_x <- function(Epsilon, x) {
  result <- Epsilon * x + log(1+x)
  return(result)
}

## Internal function in rejection sampler
#' @keywords internal
f_L_z <- function(selected_interval, Epsilon, z, a, b, lambda2, lambda3) {
  result <- ifelse(selected_interval == 1, log(1+z),
                   ifelse(selected_interval == 2, f_e_x(Epsilon, a/Epsilon) +
                            lambda2 * (z - a/Epsilon),
                          ifelse(selected_interval == 3,
                                 f_e_x(Epsilon, 1/Epsilon) + lambda3 *
                                   (z - b/Epsilon),
                                 f_e_x(Epsilon, b/Epsilon) + Epsilon *
                                   (z - b/Epsilon))))
  return(result)
}

## Internal function in rejection sampler
#' @keywords internal
sample_z <- function(selected_interval, Epsilon, lambda2, lambda3, u, a, b) {
  result <- ifelse(selected_interval == 1, (1+a/Epsilon)^u - 1,
                   ifelse(selected_interval == 2,
                          a/Epsilon -
                            log(1 - u + u *
                                  exp(-(f_e_x(Epsilon, 1/Epsilon) -
                                          f_e_x(Epsilon, a/Epsilon)))) /
                            lambda2,
                          ifelse(selected_interval == 3, 1/Epsilon -
                                   log(1 - u +  u *
                                         exp(-(f_e_x(Epsilon, b/Epsilon) -
                                                 f_e_x(Epsilon, 1/Epsilon)))) /
                                   lambda3, b/Epsilon - log(1 - u)/Epsilon)))
  return(result)
}

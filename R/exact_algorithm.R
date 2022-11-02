#' Run exact algorithm
#' @title Exact algorithm
#' @param W predictor matrix
#' @param z response variance
#' @param iteration Number of iterations
#' @param a,b rejection sampler parameter
#' @param s xi's hyperparameter
#' @param fast_algorithm When set to a true value, a calculation is performed to reduce time complexity.
#' @param time_check Output run time
#' @param iteration_check Output iteration
#' @param beta the initial value of beta
#' @param xi the initial value of Global shrinkage parameter
#' @param Sigma the inial value of Sigma
#' @param w Sigma's hyperparameter
#' @return beta posterior samples
#' @importFrom MASS mvrnorm
#' @importFrom invgamma rinvgamma
#' @export
exact_algorithm <- function(W, z, iteration = 1000, a = 1/5, b = 10, s = 0.1,
                            time_check = FALSE, iteration_check = TRUE,
                            beta = matrix(data = rep(1, times = p), nrow = 1),
                            xi = 1, sigma = 1, w = 1) {

  # time_check
  if (time_check == TRUE) {

    start_time <- Sys.time()

  }

  # data size
  N <- nrow(W)
  p <- ncol(W)

  # sampling start
  i <- 1
  while(i <= iteration) {

    # 1. eta sampling
    eta <- rejection_sampler((beta[i, ]^2)*xi/(2 * sigma), a, b)

    # D matrix, M matrix
    D <- diag(1/eta)
    WDW <- W %*% D %*% t(W)
    inverse_M <- solve(diag(N) + WDW/xi)


    # 2. xi sampling
    log_xi <- rnorm(1, mean = log(xi), sd = s)
    new_xi <- exp(log_xi)
    new_inverse_M <- solve(diag(N) + WDW/new_xi)

    # acceptance probability
    acceptance_probability <- probability_a(N,
                                            xi, new_xi,
                                            inverse_M, new_inverse_M,
                                            z, w)
    u <- runif(n = 1, min = 0, max = 1)

    # new xi accept/reject process
    if (u < acceptance_probability) {

      xi <- new_xi
      inverse_M <- new_inverse_M

    }

    # 3. sigma sampling
    sigma <- rinvgamma(1,
                       shape = (w+N)/2,
                       rate = (w + t(z) %*% inverse_M %*% z)/2)


    # 4. beta sampling : using u, f, and v
    u <- mvrnorm(n = 1, mu = rep(0, p), Sigma = D/xi)
    f <- mvrnorm(n = 1, mu = rep(0, N), Sigma = diag(N))
    v <- W %*% u + f
    v_star <- inverse_M %*% (z /sqrt(sigma) - v)

    new_beta <- sqrt(sigma) * (u + D %*% t(W) %*% v_star / xi)

    # save the sampled new_beta
    beta <- rbind(beta, t(new_beta))

    if (iteration_check == TRUE) {

      cat("iteration : ", i, "\n")

    }

    i <- i + 1

  }

  if (time_check == TRUE) {

    cat("총 소요시간 : ", Sys.time() - start_time, "\n")

  }

  # return
  beta

}

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
exact_algorithm <- function(W, z, iteration = 1000, a = 1/5, b = 10,
                            s = 0.1, xi = 1, sigma = 1, w = 1,
                            max_Epsilon = 10^(8), step_check = FALSE) {


  # data size
  N <- nrow(W)
  p <- ncol(W)

  # parameters
  beta <- matrix(0, nrow = iteration+1, ncol = p)
  global_shrinkage_parameters <- rep(0, iteration)
  local_shrinkage_parameters <- matrix(0, nrow = iteration, ncol = p)
  sigma_parameters <- rep(0, iteration)

  # 임시 테스트용 -------------------------------------------------------------
  l0 <- rep(0, p)
  l1 <- rep(1, N)
  l2 <- rep(1, p)

  if (p > N) {
    lambda_star <- sqrt(1/xi) * 1
    U <- as.numeric(lambda_star^2) * t(W)
    u <- stats::rnorm(l2, l0, lambda_star)
    v <- W %*% u + stats::rnorm(N)
    v_star <- solve((W %*% U + diag(N)), ((z/sqrt(sigma)) - v))
    beta[1, ] <- sqrt(sigma) * (u + U %*% v_star)
  }
  else {
    Q_star <- t(W) %*% W
    lambda_star <- sqrt(1/xi) * 1
    L <- chol((1/sigma) * (Q_star + diag(1/as.numeric(lambda_star^2), p, p)))
    v <- solve(t(L), t(t(z) %*% W)/sigma)
    mu <- solve(L, v)
    u <- solve(L, stats::rnorm(p))
    beta[1, ] <- mu + u
  }
  # 임시 테스트용 -------------------------------------------------------------

  if (step_check == TRUE) {

    step_checks <- data.frame(matrix(rep(0, 5), nrow = 1))
    colnames(step_checks) <- c("step1", "step2", "step3", "total_time")

  }

  # sampling start
  for(i in 1:iteration) {

    if(step_check == TRUE)
      iteration_start_time <- Sys.time()

    # 1. eta sampling
    eta <- rejection_sampler((beta[i, ]^2)*xi/(2 * sigma), a, b, max_Epsilon)

    # step1 끝낸 시간
    if(step_check == TRUE)
      step1_time <- Sys.time()

    # 2. xi sampling
    if(s != 0) {

      log_xi <- rnorm(1, mean = log(xi), sd = sqrt(s))
      new_xi <- exp(log_xi)

    }

    dw <- (1/eta) * t(W)
    WDW <- W %*% dw
    M <- diag(N) + WDW/xi
    m <- solve(M, z)
    zmz <- t(z) %*% m

    if (s != 0) {

      new_M <- diag(N) + WDW/new_xi
      k <- sqrt(det(solve(new_M, M)))
      new_m <- solve(new_M, z)
      new_zmz <- t(z) %*% new_m

      acceptance_probability <- probability_a(N, xi, new_xi, k, zmz, new_zmz, w)
      u <- runif(n = 1, min = 0, max = 1)

      # new xi accept/reject process
      if (u < acceptance_probability) {

        xi <- new_xi
        zmz <- new_zmz
        M <- new_M

      }

    }

    if(step_check == TRUE)
      step2_time <- Sys.time()

    # 3. sigma sampling
    sigma <- rinvgamma(1,
                       shape = (w+N)/2,
                       rate = (w + zmz)/2)

    # D matrix
    diagonal <- eta * xi
    diagonal_delta <- 1/diagonal

    # 4. beta sampling : using u, f, and v
    u <- rnorm(n = p, mean = 0, sd = sqrt(1/diagonal))
    f <- rnorm(n = N, mean = 0, sd = 1)
    v <- W %*% u + f
    U <- diagonal_delta * t(W)
    v_star <- solve(M, (z / sqrt(sigma) - v))
    new_beta <- sqrt(sigma) * (u + U %*% v_star)

    if(step_check == TRUE)
      step3_time <- Sys.time()

    # save the sampled value
    beta[i+1, ] <- new_beta
    local_shrinkage_parameters[i, ] <- eta
    global_shrinkage_parameters[i] <- xi
    sigma_parameters[i] <- sigma

    if ((i %% 50) == 0) {

      cat("iteration : ", i, "global : ", xi, "\n")

    }

    if (step_check == TRUE) {

      S <- length(active_set_column_index)
      step1 <- step1_time - iteration_start_time
      step2 <- step2_time - step1_time
      step3 <- step3_time - step2_time
      total <- step3_time - iteration_start_time

      step_checks[i, ] <- c(step1, step2, step3, total)

    }

  }


  # return
  if(step_check == TRUE) {

    return(list(beta = beta[-1, ],
                global_shrinkage_parameter = global_shrinkage_parameters,
                local_shrinkage_parameter = local_shrinkage_parameters,
                sigma_parameter = sigma_parameters,
                spand_time = step_checks))

  } else {

    return(list(beta = beta[-1, ],
                global_shrinkage_parameter = global_shrinkage_parameters,
                local_shrinkage_parameter = local_shrinkage_parameters,
                sigma_parameter = sigma_parameters))

  }

}

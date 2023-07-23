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
exact_algorithm <- function(W, z, iteration = 5000, a = 1/5, b = 10,
                            s = 0.8, xi = 1, sigma = 1, w = 1,
                            step_check = FALSE) {


  # data size
  N <- nrow(W)
  p <- ncol(W)

  # initial values
  eta <- rep(1, p)
  Q <- t(W) %*% W

  # parameters
  beta <- matrix(0, nrow = iteration, ncol = p)
  global_shrinkage_parameters <- rep(0, iteration)
  local_shrinkage_parameters <- matrix(0, nrow = iteration, ncol = p)
  sigma_parameters <- rep(0, iteration)

  if (step_check == TRUE) {

    step_checks <- data.frame(matrix(rep(0, 5), nrow = 1))
    colnames(step_checks) <- c("step1", "step2", "step3", "step4", "total_time")

  }

  # sampling start
  for(i in 1:iteration) {

    if(step_check == TRUE)
      iteration_start_time <- Sys.time()

    # 1. xi sampling
    log_xi <- rnorm(1, mean = log(xi), sd = sqrt(s))
    new_xi <- exp(log_xi)

    dw <- (1/eta) * t(W)
    WDW <- W %*% dw
    M <- diag(N) + WDW/xi
    m <- solve(M, z)
    zmz <- t(z) %*% m

    if (s != 0) {

      new_M <- diag(N) + WDW/new_xi
      new_m <- solve(new_M, z)
      new_zmz <- t(z) %*% new_m

      k <- sqrt(prod(diag(chol(M))^2 / diag(chol(new_M))^2))
      acceptance_probability <- probability_a(N, xi, new_xi, k, zmz, new_zmz, w)
      u <- runif(n = 1, min = 0, max = 1)

      # new xi accept/reject process
      if (u < acceptance_probability) {

        xi <- new_xi
        zmz <- new_zmz
        M <- new_M

      }

    }

    # step1 끝낸 시간
    if(step_check == TRUE)
      step1_time <- Sys.time()

    # 2. sigma sampling
    sigma <- rinvgamma(1,
                       shape = (w+N)/2,
                       rate = (w + zmz)/2)

    if(step_check == TRUE)
      step2_time <- Sys.time()

    # 3. beta sampling : using u, f, and v
    diagonal <- eta * xi
    diagonal_delta <- 1/diagonal

    u <- rnorm(n = p, mean = 0, sd = sqrt(1/diagonal))
    f <- rnorm(n = N, mean = 0, sd = 1)
    v <- W %*% u + f
    U <- diagonal_delta * t(W)
    v_star <- solve(M, (z / sqrt(sigma) - v))
    new_beta <- sqrt(sigma) * (u + U %*% v_star)

    # save the sampled value
    beta[i, ] <- new_beta
    local_shrinkage_parameters[i, ] <- eta
    global_shrinkage_parameters[i] <- xi
    sigma_parameters[i] <- sigma

    if(step_check == TRUE)
      step3_time <- Sys.time()

    # 1. eta sampling
    eta <- rejection_sampler((beta[i, ]^2)*xi/(2 * sigma), a, b)
    eta <- ifelse(eta <= .Machine$double.eps, .Machine$double.eps, eta)

    if(step_check == TRUE)
      step4_time <- Sys.time()

    if ((i %% 50) == 0) {

      cat("iteration : ", i, "global : ", xi, "\n")

    }

    if (step_check == TRUE) {

      S <- length(active_set_column_index)
      step1 <- step1_time - iteration_start_time
      step2 <- step2_time - step1_time
      step3 <- step3_time - step2_time
      step4 <- step4_time - step3_time
      total <- step4_time - iteration_start_time

      step_checks[i, ] <- c(step1, step2, step3, step4, total)

    }

  }


  # return
  if(step_check == TRUE) {

    return(list(beta = beta,
                global_shrinkage_parameter = global_shrinkage_parameters,
                local_shrinkage_parameter = local_shrinkage_parameters,
                sigma_parameter = sigma_parameters,
                spand_time = step_checks))

  } else {

    return(list(beta = beta,
                global_shrinkage_parameter = global_shrinkage_parameters,
                local_shrinkage_parameter = local_shrinkage_parameters,
                sigma_parameter = sigma_parameters))

  }

}

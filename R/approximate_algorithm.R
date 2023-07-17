#' Run approximate algorithm
#' @title Approximate algorithm
#' @param W design matrix
#' @param z response variance
#' @param iteration Number of iterations
#' @param a,b rejection sampler parameter
#' @param s xi's hyperparameter
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
approximate_horseshoe <- function(W, z, iteration = 1000, a = 1/5, b = 10,
                                  s = 0.8, xi = 1, sigma = 1, w = 1,
                                  tolerance=1e-20, step_check = FALSE) {

  # data size
  N <- nrow(W)
  p <- ncol(W)

  # threshold
  if (p >= N)
    threshold <- p
  else
    threshold <- sqrt(N*p)

  # initial values
  eta <- rep(1, p)
  Q <- t(W) %*% W

  # outputs
  beta <- matrix(0, nrow = iteration, ncol = p)
  local_shrinkage_parameters <- matrix(0, nrow = iteration, ncol = p)
  global_shrinkage_parameter <- rep(0, iteration)
  sigma_parameters <- rep(0, iteration)

  if (step_check == TRUE) {

    step_checks <- data.frame(matrix(rep(0, 6), nrow = 1))
    colnames(step_checks) <- c("total_active_column",
                               "step1", "step2", "step3", "step4", "total_time")

  }

  # sampling start
  for(i in 1:iteration) {

    if(step_check == TRUE)
      iteration_start_time <- Sys.time()

    # 1. xi sampling
    log_xi <- rnorm(1, mean = log(xi), sd = s)
    new_xi <- exp(log_xi)
    max_xi <- max(xi, new_xi)

    # active W matrix
    active_index <- which((eta * max_xi < threshold))
    S <- length(active_index)
    Q_s <- Q[active_index, active_index, drop = FALSE]
    W_s <- W[, active_index, drop = FALSE]

    # s < N case
    if (S < N) {

      wz <- t(W_s) %*% z
      z_square <- t(z) %*% z
      Q_star <- xi * diag(eta[active_index], nrow = S) + Q_s
      m <- solve(Q_star, wz)
      zmz <- z_square - t(z) %*% W_s %*% m

      if (s != 0) {

        new_Q_star <- new_xi * diag(eta[active_index], nrow = S) + Q_s
        new_m <- solve(new_Q_star, wz)
        new_zmz <- z_square - t(z) %*% W_s %*% new_m
        # new xi accept/reject process
        k <- sqrt(prod(diag(chol(Q_star))^2 / diag(chol(new_Q_star))^2  * new_xi / xi))
        acceptance_probability <- probability_a(N, xi, new_xi, k, zmz, new_zmz, w)
        u <- runif(n = 1, min = 0, max = 1)
        if (u < acceptance_probability) {

          xi <- new_xi
          zmz <- new_zmz
          Q_star <- new_Q_star

        }

      }

      # s >= N인 경우 inverse_M 계산
    } else {

      dw <- (1/eta[active_index]) * t(W_s)
      WDW <- W_s %*% dw
      M <- diag(N) + WDW/xi
      m <- solve(M, z)
      zmz <- t(z) %*% m

      if(s != 0) {

        new_M <- diag(N) + WDW/new_xi
        new_m <- solve(new_M, z)
        new_zmz <- t(z) %*% new_m
        # new xi accept/reject process
        k <- sqrt(prod(diag(chol(M))^2 / diag(chol(new_M))^2))
        acceptance_probability <- probability_a(N, xi, new_xi, k, zmz, new_zmz, w)
        u <- runif(n = 1, min = 0, max = 1)
        if (u < acceptance_probability) {

          xi <- new_xi
          zmz <- new_zmz
          M <- new_M

        }

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
    diagonal_delta[-active_index] <- 0

    u <- rnorm(n = p, mean = 0, sd = sqrt(1/diagonal))
    f <- rnorm(n = N, mean = 0, sd = 1)
    v <- W %*% u + f
    U <- diagonal_delta * t(W)

    if (S < N) {

      zv <- z / sqrt(sigma) - v
      wzv <- t(W_s) %*% zv
      m <- solve(Q_star, wzv, tol=tolerance)
      m_star <- zv - W_s %*% m
      new_beta <- sqrt(sigma) * (u + U %*% m_star)

    } else {

      v_star <- solve(M, (z / sqrt(sigma) - v), tol=tolerance)
      new_beta <- sqrt(sigma) * (u + U %*% v_star)

    }

    # save the sampled value
    beta[i, ] <- new_beta
    local_shrinkage_parameters[i, ] <- eta
    global_shrinkage_parameter[i] <- xi
    sigma_parameters[i] <- sigma

    if(step_check == TRUE)
      step3_time <- Sys.time()

    # 4. eta sampling
    eta <- rejection_sampler((beta[i, ]^2)*xi/(2 * sigma), a, b)
    eta <- ifelse(eta <= .Machine$double.eps, .Machine$double.eps, eta)

    if(step_check == TRUE)
      step4_time <- Sys.time()

    if ((i %% 50) == 0) {

      cat("iteration : ", i,
          "active : ", length(active_index),
          "global : ", xi, "\n")

    }

    if (step_check == TRUE) {

      S <- length(active_index)
      step1 <- step1_time - iteration_start_time
      step2 <- step2_time - step1_time
      step3 <- step3_time - step2_time
      step4 <- step4_time - step3_time
      total <- step4_time - iteration_start_time

      step_checks[i, ] <- c(S, step1, step2, step3, step4, total)

    }

  }

  # return
  if(step_check == TRUE) {

    return(list(beta = beta,
                global_shrinkage_parameter = global_shrinkage_parameter,
                local_shrinkage_parameter = local_shrinkage_parameters,
                sigma_parameter = sigma_parameters,
                spend_time = step_checks))

  } else {

    return(list(beta = beta,
                global_shrinkage_parameter = global_shrinkage_parameter,
                local_shrinkage_parameter = local_shrinkage_parameters,
                sigma_parameter = sigma_parameters))

  }

}

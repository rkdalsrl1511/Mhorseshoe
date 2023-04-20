#' Run approximate algorithm
#' @title Approximate algorithm
#' @param W predictor matrix
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
approximate_algorithm <- function(W, z, iteration = 1000, a = 1/5, b = 10,
                                  s = 0.01, xi = 1, sigma = 1, w = 0,
                                  step_check = FALSE) {

  # data size
  N <- nrow(W)
  p <- ncol(W)
  max_xi <- 1

  # threshold
  if (p >= N) {

    threshold <- p

  }else {

    threshold <- sqrt(N*p)

  }

  # parameters
  beta <- matrix(0.01, nrow = iteration+1, ncol = p)
  global_shrinkage_parameters <- rep(0, iteration)
  local_shrinkage_parameters <- matrix(0, nrow = iteration, ncol = p)
  sigma_parameters <- rep(0, iteration)

  if (step_check == TRUE) {

    step_checks <- data.frame(matrix(rep(0, 5), nrow = 1))
    colnames(step_checks) <- c("total_active_column",
                               "step1", "step2", "step3", "total_time")

  }

  # sampling start
  for(i in 1:iteration) {

    if(step_check == TRUE)
      iteration_start_time <- Sys.time()

    # 1. eta sampling
    eta <- rejection_sampler((beta[i, ]^2)*xi/(2 * sigma), a, b)

    # active W matrix
    active_set_column_index <- which((eta * max_xi < threshold))
    S <- length(active_set_column_index)
    W_s <- W[, active_set_column_index, drop = FALSE]

    # step1 끝낸 시간
    if(step_check == TRUE)
      step1_time <- Sys.time()

    # 2. xi sampling
    log_xi <- rnorm(1, mean = log(xi), sd = sqrt(s))
    new_xi <- exp(log_xi)
    max_xi <- max(xi, new_xi)

    # s < N인 경우 inverse_M 계산
    if (S < N) {

      Q <- t(W_s) %*% W_s
      Q_star <- xi * diag(eta[active_set_column_index], nrow = S) + Q
      new_Q_star <- new_xi * diag(eta[active_set_column_index], nrow = S) + Q

      k <- sqrt(det(solve(new_Q_star, Q_star) * new_xi / xi))

      wz <- t(W_s) %*% z
      m <- solve(Q_star, wz)
      new_m <- solve(new_Q_star, wz)
      z_square <- t(z) %*% z
      zmz <- z_square - t(z) %*% W_s %*% m
      new_zmz <- z_square - t(z) %*% W_s %*% new_m

      acceptance_probability <- probability_a(N, xi, new_xi, k, zmz, new_zmz, w)
      u <- runif(n = 1, min = 0, max = 1)

      # new xi accept/reject process
      if (u < acceptance_probability) {

        xi <- new_xi
        zmz <- new_zmz
        Q_star <- new_Q_star

      }

      # s >= N인 경우 inverse_M 계산
    } else {

      # M matrix
      dw <- (1/eta[active_set_column_index]) * t(W_s)
      WDW <- W_s %*% dw
      M <- diag(N) + WDW/xi
      new_M <- diag(N) + WDW/new_xi

      k <- sqrt(det(solve(new_M, M)))

      m <- solve(M, z)
      new_m <- solve(new_M, z)
      zmz <- t(z) %*% m
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
    diagonal_delta[-active_set_column_index] <- 0

    # 4. beta sampling : using u, f, and v
    u <- rnorm(n = p, mean = 0, sd = sqrt(1/diagonal))
    f <- rnorm(n = N, mean = 0, sd = 1)
    v <- W %*% u + f
    U <- diagonal_delta * t(W)

    if (S < N) {

      zv <- z / sqrt(sigma) - v
      wzv <- t(W_s) %*% zv
      m <- solve(Q_star, wzv)
      m_star <- zv - W_s %*% m
      new_beta <- sqrt(sigma) * (u + U %*% m_star)

    } else {

      v_star <- solve(M, (z / sqrt(sigma) - v))
      new_beta <- sqrt(sigma) * (u + U %*% v_star)

    }

    if(step_check == TRUE)
      step3_time <- Sys.time()

    # save the sampled value
    beta[i+1, ] <- new_beta
    local_shrinkage_parameters[i, ] <- eta
    global_shrinkage_parameters[i] <- xi
    sigma_parameters[i] <- sigma

    if ((i %% 50) == 0) {

      cat("iteration : ", i,
          "active : ", length(active_set_column_index),
          "global : ", xi, "\n")

    }

    if (step_check == TRUE) {

      S <- length(active_set_column_index)
      step1 <- step1_time - iteration_start_time
      step2 <- step2_time - step1_time
      step3 <- step3_time - step2_time
      total <- step3_time - iteration_start_time

      step_checks[i, ] <- c(S, step1, step2, step3, total)

    }

  }

  # return
  if(step_check == TRUE) {

    return(list(beta = beta[-1, ],
                global_shrinkage_parameter = global_shrinkage_parameters,
                local_shrinkage_parameter = local_shrinkage_parameters,
                sigma_parameter = sigma_parameters,
                spend_time = step_checks))

  } else {

    return(list(beta = beta[-1, ],
                global_shrinkage_parameter = global_shrinkage_parameters,
                local_shrinkage_parameter = local_shrinkage_parameters,
                sigma_parameter = sigma_parameters))

  }

}

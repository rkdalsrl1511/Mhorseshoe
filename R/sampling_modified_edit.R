# Run modified approximate algorithm with dependence prior
#' @importFrom invgamma rinvgamma
#' @export
sampling_modified_edit <- function(W, z, xi = 1, sigma = 1, iteration = 5000,
                                   a = 1/5, b = 10, max_Epsilon = 10^(8), s = 1, w = 1,
                                   m_eff = NULL, t = 50, alpha0 = -0.5, alpha1 = -7*10^(-4),
                                   step_check = FALSE) {

  # data size
  N <- nrow(W)
  p <- ncol(W)

  # initial values
  S <- p
  active_set_column_index <- 1:p
  if (is.null(m_eff))
    m_eff <- p
  eta <- rep(1, p)
  Q <- t(W) %*% W
  s2.vec <- diag(Q)
  if (p < N)
    Q_star <- xi * diag(eta) + Q
  else {
    dw <- (1/eta) * t(W)
    WDW <- W %*% dw
    M <- diag(N) + WDW/xi
  }

  # parameters
  beta <- matrix(0, nrow = iteration, ncol = p)
  local_shrinkage_parameters <- matrix(0, nrow = iteration, ncol = p)
  global_shrinkage_parameter <- rep(0, iteration)
  sigma_parameters <- rep(0, iteration)
  active_sets <- rep(0, iteration)
  meffs <- rep(0, iteration)

  if (step_check == TRUE) {

    step_checks <- data.frame(matrix(rep(0, 6), nrow = 1))
    colnames(step_checks) <- c("number of active set columns",
                               "step1", "step2", "step3", "step4", "total_time")

  }

  # sampling start
  for(i in 1:iteration) {

    if (step_check == TRUE)
      iteration_start_time <- Sys.time()

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
      wzv <- t(W[, active_set_column_index, drop = FALSE]) %*% zv
      m <- solve(Q_star, wzv)
      m_star <- zv - W[, active_set_column_index, drop = FALSE] %*% m
      new_beta <- sqrt(sigma) * (u + U %*% m_star)

    } else {

      v_star <- solve(M, (z / sqrt(sigma) - v))
      new_beta <- sqrt(sigma) * (u + U %*% v_star)

    }

    # save the sampled value
    beta[i, ] <- new_beta
    local_shrinkage_parameters[i, ] <- eta
    global_shrinkage_parameter[i] <- xi
    sigma_parameters[i] <- sigma
    active_sets[i] <- S
    meffs[i] <- m_eff

    if(step_check == TRUE)
      step1_time <- Sys.time()

    # 1. eta sampling
    eta <- rejection_sampler((beta[i, ]^2)*xi/(2 * sigma), a, b, max_Epsilon)

    # meff 계산 및 threshold 수정
    if (i %% t == 0) {

      u_i <- runif(1,0,1)
      p_i <- exp(alpha0 + alpha1 * i)

      if (u_i < p_i)
        m_eff <- sum(1/((eta*xi)/s2.vec + 1))

    }

    threshold <- sort(eta)[ceiling(m_eff)]
    active_set_column_index <- which(eta <= threshold)
    S <- length(active_set_column_index)

    if(step_check == TRUE)
      step2_time <- Sys.time()

    # 2. xi sampling
    log_xi <- rnorm(1, mean = log(xi), sd = sqrt(s))
    new_xi <- exp(log_xi)

    # s < N인 경우 inverse_M 계산
    if (S < N) {

      Q_star <- xi * diag(eta[active_set_column_index], nrow = S) + Q[active_set_column_index, active_set_column_index]
      new_Q_star <- new_xi * diag(eta[active_set_column_index], nrow = S) + Q[active_set_column_index, active_set_column_index]
      k <- sqrt(det(solve(new_Q_star, Q_star) * new_xi / xi))
      wz <- t(W[, active_set_column_index, drop = FALSE]) %*% z
      m <- solve(Q_star, wz)
      new_m <- solve(new_Q_star, wz)
      z_square <- t(z) %*% z
      zmz <- z_square - t(z) %*% W[, active_set_column_index, drop = FALSE] %*% m
      new_zmz <- z_square - t(z) %*% W[, active_set_column_index, drop = FALSE] %*% new_m

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
      dw <- (1/eta[active_set_column_index]) * t(W[, active_set_column_index, drop = FALSE])
      WDW <- W[, active_set_column_index, drop = FALSE] %*% dw
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
      step3_time <- Sys.time()

    # 3. sigma sampling
    sigma <- rinvgamma(1,
                       shape = (w+N)/2,
                       rate = (w + zmz)/2)

    if(step_check == TRUE)
      step4_time <- Sys.time()

    if ((i %% 50) == 0) {

      cat("iteration : ", i,
          "active : ", active_sets[i],
          "global : ", global_shrinkage_parameter[i], "\n")

    }

    # 상세 시간 체크 옵션
    if (step_check == TRUE) {

      step1 <- step1_time - iteration_start_time
      step2 <- step2_time - step1_time
      step3 <- step3_time - step2_time
      step4 <- step4_time - step3_time
      total_time <- step4_time - iteration_start_time

      step_checks[i, ] <- c(active_sets[i], step1, step2, step3, step4, total_time)

    }

  }

  # return
  if(step_check == TRUE){

    return(list(beta = beta[-1, ],
                local_shrinkage_parameter = local_shrinkage_parameters,
                global_shrinkage_parameter = global_shrinkage_parameter,
                sigma2 = sigma_parameters,
                active_set = active_sets,
                meff = meffs,
                spand_time = step_checks))

  } else {

    return(list(beta = beta[-1, ],
                local_shrinkage_parameter = local_shrinkage_parameters,
                global_shrinkage_parameter = global_shrinkage_parameter,
                sigma2 = sigma_parameters,
                meff = meffs,
                active_set = active_sets))

  }

}

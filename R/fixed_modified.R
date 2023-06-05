# Run modified approximate algorithm with fixed global shrinkage parameter and dependence prior
#' @importFrom invgamma rinvgamma
#' @export
fixed_modified <- function(W, z, xi = 1, sigma = 1, iteration = 5000,
                           a = 1/5, b = 10, max_Epsilon = 10^(30), w = 1,
                           t = 50, alpha0 = 0, alpha1 = -9*10^(-4),
                           step_check = FALSE) {

  # data size
  N <- nrow(W)
  p <- ncol(W)

  # initial values
  beta <- matrix(0, nrow = iteration+1, ncol = p)
  Q <- t(W) %*% W
  s2.vec <- diag(Q)
  m_eff <- p

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
    lambda_star <- sqrt(1/xi) * 1
    L <- chol((1/sigma) * (Q + diag(1/as.numeric(lambda_star^2), p, p)))
    v <- solve(t(L), t(t(z) %*% W)/sigma)
    mu <- solve(L, v)
    u <- solve(L, stats::rnorm(p))
    beta[1, ] <- mu + u
  }
  # 임시 테스트용 -------------------------------------------------------------

  # parameters
  local_shrinkage_parameters <- matrix(0, nrow = iteration, ncol = p)
  global_shrinkage_parameter <- rep(0, iteration)
  sigma_parameters <- rep(0, iteration)
  active_sets <- rep(0, iteration)
  meffs <- rep(0, iteration)

  if (step_check == TRUE) {

    step_checks <- data.frame(matrix(rep(0, 5), nrow = 1))
    colnames(step_checks) <- c("number of active set columns",
                               "step1", "step2", "step3", "total_time")

  }

  # sampling start
  for(i in 1:iteration) {

    if (step_check == TRUE)
      iteration_start_time <- Sys.time()

    # 1. eta sampling
    eta <- rejection_sampler((beta[i, ]^2)*xi/(2 * sigma), a, b, max_Epsilon)

    if(step_check == TRUE)
      step1_time <- Sys.time()

    # meff 계산 및 threshold 수정
    if (i %% t == 0) {

      u_i <- runif(1,0,1)
      p_i <- exp(alpha0 + alpha1 * i)

      if (u_i < p_i)
        m_eff <- sum(1/(eta*xi/s2.vec + 1))

    }

    threshold <- sort(eta)[ceiling(m_eff)]
    active_set_column_index <- which(eta <= threshold)
    S <- length(active_set_column_index)

    diagonal <- (eta*xi)
    diagonal_delta <- 1/diagonal
    diagonal_delta[-active_set_column_index] <- 0

    # active W matrix
    W_s <- W[, active_set_column_index, drop = FALSE]

    # set inverse_M %*% z, v_star
    u <- rnorm(n = p, mean = 0, sd = sqrt(1/diagonal))
    f <- rnorm(n = N, mean = 0, sd = 1)
    v <- W %*% u + f
    U <- diagonal_delta * t(W)

    if (S > N) {
      Q <- W_s %*% U[active_set_column_index, ]
      M <- diag(N) + Q
      inv_mz <- solve(M, z)
    } else {
      Q_star <- Q[active_set_column_index, active_set_column_index, drop = FALSE] + diag(diagonal[active_set_column_index], nrow = S)
      inv_mz <- z - W_s %*% solve(Q_star, t(W_s) %*% z)
    }

    # 2. sigma sampling
    sigma <- rinvgamma(1,
                       shape = (w+N)/2,
                       rate = (w + z %*% inv_mz)/2)

    if(step_check == TRUE)
      step2_time <- Sys.time()

    # 3. beta sampling
    if (S > N) {
      v_star <- inv_mz / sqrt(sigma) - solve(M, v)
    } else {
      v_star <- inv_mz / sqrt(sigma) - v + W_s %*% solve(Q_star, t(W_s) %*% v)
    }

    new_beta <- sqrt(sigma) * (u + U %*% v_star)

    if(step_check == TRUE)
      step3_time <- Sys.time()

    # save the sampled value
    beta[i+1, ] <- new_beta
    local_shrinkage_parameters[i, ] <- eta
    global_shrinkage_parameter[i] <- xi
    sigma_parameters[i] <- sigma
    active_sets[i] <- S
    meffs[i] <- m_eff

    if ((i %% 50) == 0) {

      cat("iteration : ", i, ", active : ", active_sets[i], "\n")

    }

    # 상세 시간 체크 옵션
    if (step_check == TRUE) {

      step1 <- step1_time - iteration_start_time
      step2 <- step2_time - step1_time
      step3 <- step3_time - step2_time
      total_time <- step3_time - iteration_start_time

      step_checks[i, ] <- c(S, step1, step2, step3, total_time)

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
                spend_time = step_checks))

  } else {

    return(list(beta = beta[-1, ],
                local_shrinkage_parameter = local_shrinkage_parameters,
                global_shrinkage_parameter = global_shrinkage_parameter,
                sigma2 = sigma_parameters,
                meff = meffs,
                active_set = active_sets))

  }

}

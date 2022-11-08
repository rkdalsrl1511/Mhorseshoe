#' Run modified approximate algorithm with dependence prior
#' @title new approximate algorithm_c
#' @param W predictor matrix
#' @param z response variance
#' @param iteration Number of iterations
#' @param a,b rejection sampler parameters
#' @param w Sigma's hyperparameter
#' @param threshold Parameter that determines the active set
#' @param iteration_check Output iteration
#' @param alpha0 Parameters of adaptive probability
#' @param alpha1 Parameters of adaptive probability
#' @return beta posterior samples
#' @importFrom invgamma rinvgamma
#' @export
modified_approximate_algorithm <- function(W, z, iteration = 5000,
                                           a = 1/5, b = 10, w = 1,
                                           alpha0 = -0.5, alpha1 = -7*10^(-4),
                                           step_check = FALSE) {

  # data size
  N <- nrow(W)
  p <- ncol(W)

  # initial values
  beta <- matrix(0, nrow = iteration+1, ncol = p)
  xi <- 1
  sigma <- 1
  m_eff <- 0
  threshold <- N*p

  # parameters
  local_shrinkage_parameters <- matrix(0, nrow = iteration, ncol = p)
  global_shrinkage_parameter <- rep(0, iteration)
  sigma_parameters <- rep(0, iteration)
  active_sets <- rep(0, iteration)

  if (step_check == TRUE) {

    step_checks <- data.frame(matrix(rep(0, 7), nrow = 1))
    colnames(step_checks) <- c("number of active set columns",
                               "step1", "step2", "step3", "total_time")

  }

  # sampling start
  for(i in 1:iteration) {

    if (step_check == TRUE)
      iteration_start_time <- Sys.time()

    if (i == 2) {

      xi <- sqrt(N)*p^2

    }

    # 1. eta sampling
    eta <- rejection_sampler((beta[i, ]^2)*xi/(2 * sigma), a, b)

    # main concept
    if (i > 200 & m_eff != 0) {

      threshold <- sort(eta)[ceiling(m_eff)+1]
      active_set_column_index <- which(eta <= threshold)

    } else {

      active_set_column_index <- which(eta*xi < threshold)

    }

    S <- length(active_set_column_index)
    diagonal <- (eta*xi)
    diagonal_delta <- 1/diagonal
    diagonal_delta[-active_set_column_index] <- 0

    # active W matrix
    W_s <- W[, active_set_column_index]

    if(step_check == TRUE)
      step1_time <- Sys.time()

    # 2. sigma/beta sampling
    u <- rnorm(n = p, mean = 0, sd = sqrt(1/diagonal))
    f <- rnorm(n = N, mean = 0, sd = 1)
    v <- W %*% u + f
    U <- diagonal_delta * t(W)

    # set inverse_M %*% z, v_star
    if (S > N) {

      Q <- W_s %*% U[active_set_column_index, ]
      M <- diag(N) + Q
      inv_mz <- solve(M, z)
      v_star <- inv_mz / sqrt(sigma) - solve(M, v)

    } else {

      Q <- t(W_s) %*% W_s
      Q_star <- Q + diag(diagonal[active_set_column_index])
      inv_mz <- z - W_s %*% solve(Q_star, t(W_s) %*% z)
      v_star <- inv_mz / sqrt(sigma) - v + W_s %*% solve(Q_star, t(W_s) %*% v)

    }

    sigma <- rinvgamma(1,
                       shape = (w+N)/2,
                       rate = (w + z %*% inv_mz)/2)
    new_beta <- sqrt(sigma) * (u + U %*% v_star)

    if(step_check == TRUE)
      step2_time <- Sys.time()

    # meff 계산
    if (i %% 20 == 0) {

      u_i <- runif(1,0,1)
      p_i <- exp(alpha0 + alpha1 * i)

      if (u_i < p_i) {

        if (S > N) {

          k_star <- solve(M, Q)

        } else {

          k_star <- solve(Q_star, Q)

        }

        m_eff <- sum(diag(k_star))

      }

    }

    if(step_check == TRUE)
      step3_time <- Sys.time()

    # save the sampled value
    beta[i+1, ] <- new_beta
    local_shrinkage_parameters[i, ] <- eta
    global_shrinkage_parameter[i] <- xi
    sigma_parameters[i] <- sigma
    active_sets[i] <- S

    if ((i %% 500) == 0) {

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
                spand_time = step_checks))

  } else {

    return(list(beta = beta[-1, ],
                local_shrinkage_parameter = local_shrinkage_parameters,
                global_shrinkage_parameter = global_shrinkage_parameter,
                sigma2 = sigma_parameters,
                active_set = active_sets))

  }

}

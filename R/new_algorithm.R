#' Run modified approximate algorithm
#' @title new approximate algorithm
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
                                           a = 1/5, b = 10, w = 0, threshold = 0,
                                           alpha0 = -0.5, alpha1 = -7*10^(-4),
                                           step_check = FALSE) {

  # data size
  N <- nrow(W)
  p <- ncol(W)

  # initial value : beta, global shrinkage parameter, sigma
  beta <- matrix(0, nrow = iteration+1, ncol = p)
  xi <- 1
  sigma <- 1
  m_eff <- 0

  # parameters
  local_shrinkage_parameters <- matrix(0, nrow = iteration, ncol = p)
  global_shrinkage_parameter <- rep(0, iteration)
  sigma_parameters <- rep(0, iteration)
  m_effs <- rep(0, iteration)
  active_sets <- rep(0, iteration)

  if (step_check == TRUE) {

    step_checks <- data.frame(matrix(rep(0, 7), nrow = 1))
    colnames(step_checks) <- c("total_active_column",
                               "step1", "step2", "step3",
                               "step4", "step5", "total_time")

  }

  # sampling start
  for(i in 1:iteration) {

    if (step_check == TRUE)
      iteration_start_time <- Sys.time()

    # For error prevention
    if (i == 20) {

      xi <- N*p

    }

    # 1. eta sampling
    eta <- rejection_sampler((beta[i, ]^2)*xi/2, a, b)

    # active_set
    active_set_column_index <- which(eta*xi*sigma < threshold)
    S <- length(active_set_column_index)

    # active W matrix
    W_s <- W[, active_set_column_index]

    if(step_check == TRUE)
      step1_time <- Sys.time()

    # 2. sigma sampling
    e <- z - W_s %*% beta[i, active_set_column_index]
    sigma <- rinvgamma(1,
                       shape = (w + N)/2,
                       rate = (w + t(e) %*% e)/2)

    if(step_check == TRUE)
      step2_time <- Sys.time()

    # D matrix
    diagonal <- (eta*xi*sigma)
    diagonal_delta <- 1/diagonal
    diagonal_delta[-active_set_column_index] <- 0

    if(step_check == TRUE)
      step3_time <- Sys.time()

    # 3. beta sampling : using u, f, and v
    u <- rnorm(n = p, mean = 0, sd = sqrt(1/diagonal))
    f <- rnorm(n = N, mean = 0, sd = 1)
    v <- W %*% u + f
    U <- diagonal_delta * t(W)

    if (S > N) {

      Q <- W_s %*% U[active_set_column_index, ]
      v_star <- solve((Q + diag(N)), ((z / sqrt(sigma)) - v))
      new_beta <- sqrt(sigma) * (u + U %*% v_star)

    } else {

      D <- diag(diagonal[active_set_column_index])
      Q_star <- t(W_s) %*% W_s
      zv <- ((z / sqrt(sigma)) - v)
      wzv <- t(W_s) %*% zv
      m <- solve(Q_star + D, wzv)
      m_star <- zv - W_s %*% m
      new_beta <- sqrt(sigma) * (u + U %*% m_star)

    }

    if(step_check == TRUE)
      step4_time <- Sys.time()

    # meff 계산
    if (i %% 20 == 0) {

      u_i <- runif(1,0,1)
      p_i <- exp(alpha0 + alpha1 * i)

      if (u_i < p_i) {

        if (S > N) {

          k_star <- solve(Q + diag(N), Q)
          m_eff <- sum(diag(k_star))

        } else {

          k_star <- solve(Q_star + D, Q_star)
          m_eff <- sum(diag(k_star))

        }

        # main concept
        if (i > 200) {

        threshold <- sort(diagonal)[ceiling(m_eff) + 10]

        }

      }

    }

    # save the sampled value
    beta[i+1, ] <- new_beta
    local_shrinkage_parameters[i, ] <- eta
    global_shrinkage_parameter[i] <- xi
    sigma_parameters[i] <- sigma
    active_sets[i] <- S
    m_effs[i] <- m_eff

    if(step_check == TRUE)
      step5_time <- Sys.time()

    if ((i %% 500) == 0) {

      cat("iteration : ", i, ", active : ", active_sets[i], "\n")

    }

    # 상세 시간 체크 옵션
    if (step_check == TRUE) {

      step1 <- step1_time - iteration_start_time
      step2 <- step2_time - step1_time
      step3 <- step3_time - step2_time
      step4 <- step4_time - step3_time
      step5 <- step5_time - step4_time
      total_time <- step5_time - iteration_start_time

      step_checks[i, ] <- c(S, step1, step2, step3, step4, step5, total_time)

    }

  }

  # return
  if(step_check == TRUE){

    return(list(beta = beta[-1, ],
                local_shrinkage_parameter = local_shrinkage_parameters,
                global_shrinkage_parameter = global_shrinkage_parameter,
                sigma2 = sigma_parameters,
                m_eff = m_effs,
                active_set = active_sets,
                spand_time = step_checks))

  } else {

    return(list(beta = beta[-1, ],
                local_shrinkage_parameter = local_shrinkage_parameters,
                global_shrinkage_parameter = global_shrinkage_parameter,
                sigma2 = sigma_parameters,
                m_eff = m_effs,
                active_set = active_sets))

  }

}

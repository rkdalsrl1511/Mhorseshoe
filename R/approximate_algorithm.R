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
                                  time_check = FALSE, iteration_check = TRUE,
                                  step_check = FALSE) {

  # time_check
  if (time_check == TRUE) {

    start_time <- Sys.time()

  }

  # data size
  N <- nrow(W)
  p <- ncol(W)
  beta <- matrix(data = rep(1, times = p), nrow = 1)

  # threshold
  if (p >= N) {

    threshold <- 1/p

  }else {

    threshold <- 1/sqrt(N*p)

  }

  # parameters
  global_shrinkage_parameters <- vector(mode = "numeric", length = 0)
  local_shrinkage_parameters <- matrix(data = rep(1, times = p), nrow = 1)
  sigma_parameters <- vector(mode = "numeric", length = 0)
  beta_variances <- vector(mode = "list", length = iteration)

  # sampling start
  i <- 1
  while(i <= iteration) {

    # 시작 시간
    iteration_start_time <- Sys.time()

    # 1. eta sampling
    eta <- rejection_sampler((beta[i, ]^2)*xi/(2 * sigma), a, b)
    inverse_eta <- 1/eta

    # active_set
    if(i == 1)
      max_xi <- xi

    active_set_column_index <- which(inverse_eta/max_xi > threshold)

    # active W matrix
    W_s <- W[, active_set_column_index]

    # step1 끝낸 시간
    step1_time <- Sys.time()

    # 2. xi sampling
    log_xi <- rnorm(1, mean = log(xi), sd = sqrt(s))
    new_xi <- exp(log_xi)
    max_xi <- max(xi, new_xi)

    # s < N인 경우 inverse_M 계산
    if (length(active_set_column_index) < N) {

      WTW <- t(W_s) %*% W_s
      WW <- xi*diag(eta[active_set_column_index]) + WTW
      new_WW <- solve(new_xi*diag(eta[active_set_column_index]) + WTW)
      inverse_M <- diag(N) - W_s %*% solve(WW) %*% t(W_s)
      new_inverse_M <- diag(N) - W_s %*% new_WW %*% t(W_s)

      # acceptance probability
      acceptance_probability <- probability_a2(N,
                                               xi, new_xi,
                                               WW, new_WW,
                                               inverse_M, new_inverse_M,
                                               z, w)

      # s >= N인 경우 inverse_M 계산
    } else {

      # M matrix
      WDW <- W_s %*% diag(inverse_eta[active_set_column_index]) %*% t(W_s)
      M <- diag(N) + WDW/xi
      inverse_M <- solve(M)
      new_inverse_M <- solve(diag(N) + WDW/new_xi)

      # acceptance probability
      acceptance_probability <- probability_a(N,
                                              xi, new_xi,
                                              M, inverse_M, new_inverse_M,
                                              z, w)

    }

    u <- runif(n = 1, min = 0, max = 1)

    # new xi accept/reject process
    if (u < acceptance_probability) {

      xi <- new_xi
      inverse_M <- new_inverse_M

    }

    step2_time <- Sys.time()

    # 3. sigma sampling
    sigma <- rinvgamma(1,
                       shape = (w+N)/2,
                       rate = (w + t(z) %*% inverse_M %*% z)/2)

    step3_time <- Sys.time()

    # D matrix
    diagonal <- ifelse(inverse_eta/xi == Inf, 10^30, inverse_eta/xi)
    diagonal_delta <- diagonal
    diagonal_delta[-active_set_column_index] <- 0
    beta_variances[[i]] <- sigma * solve(t(W) %*% W + diag(eta * xi))

    # 4. beta sampling : using u, f, and v
    u <- rnorm(n = p, mean = 0, sd = sqrt(diagonal))
    f <- rnorm(n = N, mean = 0, sd = 1)
    v <- W %*% u + f
    v_star <- inverse_M %*% (z / sqrt(sigma) - v)
    wv <- t(W) %*% v_star

    for (j in 1:p) {

      wv[j, ] <- wv[j, ] * diagonal_delta[j]

    }

    new_beta <- sqrt(sigma) * (u + wv)

    # save the sampled value
    beta <- rbind(beta, t(new_beta))
    global_shrinkage_parameters <- append(global_shrinkage_parameters, xi)
    local_shrinkage_parameters <- rbind(local_shrinkage_parameters, eta)
    sigma_parameters <- append(sigma_parameters, sigma)

    step4_time <- Sys.time()

    # iteration 출력 옵션
    if (iteration_check == TRUE && (i %% 50) == 0) {

      cat("iteration : ", i,
          "active : ", length(active_set_column_index),
          "global : ", xi, "\n")

    }

    # 상세 시간 체크 옵션
    if (step_check == TRUE) {

      cat("total active column : ", length(active_set_column_index), "\n")
      cat("Time spent in step 1 : ", step1_time - iteration_start_time, "\n")
      cat("Time spent in step 2 : ", step2_time - step1_time, "\n")
      cat("Time spent in step 3 : ", step3_time - step2_time, "\n")
      cat("Time spent in step 4 : ", step4_time - step3_time, "\n")
      cat("Total time required : ", step4_time - iteration_start_time, "\n")

    }

    i <- i + 1

  }

  # time_check
  if (time_check == TRUE) {

    cat("총 소요시간 : ", Sys.time() - start_time, "\n")

  }

  # return
  list(beta, global_shrinkage_parameters, local_shrinkage_parameters, sigma_parameters, beta_variances)

}

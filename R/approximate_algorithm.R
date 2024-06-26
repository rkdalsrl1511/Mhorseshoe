#' Run approximate MCMC algorithm for horseshoe prior
#'
#' The approximate MCMC algorithm for the horseshoe prior
#'
#' This function implements the approximate algorithm introduced in Section
#' 2.2 of Johndrow et al. (2020) and the method proposed in this package, which
#' improves computation speed when p >> N. The approximate algorithm introduces
#' a threshold and uses only a portion of the total \eqn{p} columns for matrix
#' multiplication, reducing the computational cost compared to the existing
#' MCMC algorithms for the horseshoe prior. The "auto.threshold" argument
#' determines whether the threshold used in the algorithm will be updated by
#' the adaptive method proposed in this package.
#'
#' @references Johndrow, J., Orenstein, P., & Bhattacharya, A. (2020).
#' Scalable Approximate MCMC Algorithms for the Horseshoe Prior. In Journal
#' of Machine Learning Research, 21, 1-61.
#'
#' @inheritParams exact_horseshoe
#' @param auto.threshold Argument for setting whether to use an algorithm that
#'  automatically updates the threshold using adaptive probability.
#' @param threshold Threshold to be used in the approximate MCMC algorithm.
#'  This argument is ignored when auto.threshold=TRUE. If you select
#'  auto.threshold = FALSE and threshold = 0 (This is the default value for the
#'  threshold argument), the threshold is set to
#'  \eqn{\sqrt{p \times min(N, p)}} as suggested in Johndrow et al. (2020). Or,
#'  you can set your custom value directly through this argument. For more
#'  information about \eqn{\delta}, browseVignettes("Mhorseshoe") and 4.1 of
#'  Johndrow et al. (2020).
#' @param t Threshold update cycle for adaptive probability algorithm when
#'  auto.threshold is set to TRUE. The default is 10.
#' @param adapt_p0 A tuning parameter \eqn{p_{0}} of the adaptive probability,
#'  \eqn{p(t) = exp[p_{0} + p_{1}t]}. The default is \eqn{0}.
#' @param adapt_p1 A tuning parameter \eqn{a_{1}} of the adaptive probability,
#'  \eqn{p(t) = exp[p_{0} + p_{1}t]}. The default is \eqn{-4.6 \times 10^{-4}}.
#' @return \item{BetaHat}{Posterior mean of \eqn{\beta}.}
#' \item{LeftCI}{Lower bound of \eqn{100(1-\alpha)\%} credible interval for
#'  \eqn{\beta}.}
#' \item{RightCI}{Upper bound of \eqn{100(1-\alpha)\%} credible interval for
#'  \eqn{\beta}.}
#' \item{Sigma2Hat}{Posterior mean of \eqn{\sigma^{2}}.}
#' \item{TauHat}{Posterior mean of \eqn{\tau}.}
#' \item{LambdaHat}{Posterior mean of \eqn{\lambda_{j},\ j=1,2,...p.}.}
#' \item{ActiveMean}{Average number of elements in the active set per iteration
#'  in this algorithm.}
#' \item{BetaSamples}{Posterior samples of \eqn{\beta}.}
#' \item{LambdaSamples}{Posterior samples of local shrinkage parameters.}
#' \item{TauSamples}{Posterior samples of global shrinkage parameter.}
#' \item{Sigma2Samples}{Posterior samples of \eqn{sigma^{2}}.}
#' \item{ActiveSet}{\eqn{\mathbb{R}^{iter \times p}} Matrix indicating active
#'  elements as 1 and non-active elements as 0 per iteration of the MCMC
#'  algorithm.}
#'
#' @examples
#' # Making simulation data.
#' set.seed(123)
#' N <- 200
#' p <- 100
#' true_beta <- c(rep(1, 10), rep(0, 90))
#'
#' X <- matrix(1, nrow = N, ncol = p) # Design matrix X.
#' for (i in 1:p) {
#'   X[, i] <- stats::rnorm(N, mean = 0, sd = 1)
#' }
#'
#' y <- vector(mode = "numeric", length = N) # Response variable y.
#' e <- rnorm(N, mean = 0, sd = 2) # error term e.
#' for (i in 1:10) {
#'   y <- y + true_beta[i] * X[, i]
#' }
#' y <- y + e
#'
#' # Run with auto.threshold set to TRUE
#' result1 <- approx_horseshoe(y, X, burn = 0, iter = 100,
#'                             auto.threshold = TRUE)
#'
#' # Run with fixed custom threshold
#' result2 <- approx_horseshoe(y, X, burn = 0, iter = 100,
#'                             auto.threshold = FALSE, threshold = 1/(5 * p))
#'
#' # posterior mean
#' betahat <- result1$BetaHat
#'
#' # Lower bound of the 95% credible interval
#' leftCI <- result1$LeftCI
#'
#' # Upper bound of the 95% credible interval
#' RightCI <- result1$RightCI
#'
#' @export
approx_horseshoe <- function(y, X, burn = 1000, iter = 5000,
                             auto.threshold = TRUE, threshold = 0, tau = 1,
                             s = 0.8, sigma2 = 1, w = 1, alpha = 0.05, a = 0.2,
                             b = 10, t = 10, adapt_p0 = 0,
                             adapt_p1 = -4.6*10^(-4)) {
  N <- nrow(X)
  p <- ncol(X)
  eta <- rep(1, p)
  xi <- tau^(-2)
  Q <- t(X) %*% X
  nmc <- burn + iter
  if (auto.threshold == TRUE) {
    S <- p
    active_index <- 1:p
    m_eff <- p
    s2.vec <- diag(Q)
  } else if (threshold == 0) {
    threshold <- ifelse(p >= N, 1/p, 1/sqrt(N*p))
    message("You chose FALSE for the auto.threshold argument. ",
            "and since threshold = 0, set it to the default value of ",
            threshold, ".")
  }
  betaout <- matrix(0, nrow = nmc, ncol = p)
  etaout <- matrix(0, nrow = nmc, ncol = p)
  xiout <- rep(0, nmc)
  sigma2out <- rep(0, nmc)
  activeout <- matrix(0, nrow = nmc, ncol = p)
  # run
  for(i in 1:nmc) {
    log_xi <- stats::rnorm(1, mean = log(xi), sd = s)
    new_xi <- exp(log_xi)
    # when to use a fixed threshold
    if (auto.threshold == FALSE) {
      max_xi <- max(xi, new_xi)
      active_index <- which((1/(eta * max_xi) > threshold))
      S <- length(active_index)
    }
    if (S == 0) {
      warning("All coefficients in the linear model are estimated to be 0, ",
              "so the algorithm terminates and the sampling results are ",
              "returned. This problem may be caused by the threshold being ",
              "set too large. To solve this, Use the auto.threshold option, ",
              "or if you set auto.threshold == FALSE, set the threshold ",
              "argument smaller. Alternatively, run the exact_horseshoe ",
              "function instead of approx_horseshoe and check the results.")
      burn <- 0
      nmc <- i
      break
    }
    Q_s <- Q[active_index, active_index, drop = FALSE]
    X_s <- X[, active_index, drop = FALSE]
    # xi update
    if (S < N) {
      Xy <- t(X_s) %*% y
      y_square <- t(y) %*% y
      Q_star <- xi * diag(eta[active_index], nrow = S) + Q_s
      m <- solve(Q_star, Xy)
      ymy <- y_square - t(y) %*% X_s %*% m
      if (s != 0) {
        new_Q_star <- new_xi * diag(eta[active_index], nrow = S) + Q_s
        new_m <- solve(new_Q_star, Xy)
        new_ymy <- y_square - t(y) %*% X_s %*% new_m
        cM <- (diag(chol(Q_star))^2) / xi
        new_cM <- (diag(chol(new_Q_star))^2) / new_xi
        curr_ratio <- -sum(log(cM)) / 2 - ((N + w) / 2) * log(w + ymy) -
          log(sqrt(xi) * (1 + xi))
        new_ratio <- -sum(log(new_cM)) / 2 - ((N + w) / 2) * log(w + new_ymy) -
          log(sqrt(new_xi) * (1 + new_xi))
        acceptance_probability <- exp(new_ratio - curr_ratio + log(new_xi) -
                                        log(xi))
        u <- stats::runif(n = 1, min = 0, max = 1)
        if (u < acceptance_probability) {
          xi <- new_xi
          ymy <- new_ymy
          Q_star <- new_Q_star
        }
      }
    } else {
      DX <- (1/eta[active_index]) * t(X_s)
      XDX <- X_s %*% DX
      M <- diag(N) + XDX/xi
      m <- solve(M, y)
      ymy <- t(y) %*% m
      if (s != 0) {
        new_M <- diag(N) + XDX/new_xi
        new_m <- solve(new_M, y)
        new_ymy <- t(y) %*% new_m
        cM <- diag(chol(M))^2
        new_cM <- diag(chol(new_M))^2
        curr_ratio <- -sum(log(cM))/2 - ((N + w)/2) * log(w + ymy) -
          log(sqrt(xi) * (1 + xi))
        new_ratio <- -sum(log(new_cM))/2 - ((N + w)/2) * log(w + new_ymy) -
          log(sqrt(new_xi) * (1 + new_xi))
        acceptance_probability <- exp(new_ratio - curr_ratio + log(new_xi) -
                                        log(xi))
        u <- stats::runif(n = 1, min = 0, max = 1)
        if (u < acceptance_probability) {
          xi <- new_xi
          ymy <- new_ymy
          M <- new_M
        }
      }
    }
    # sigma update
    sigma2 <- 1/stats::rgamma(1, shape = (w + N)/2, rate = (w + ymy)/2)
    # beta update
    diag_D <- 1 / (eta * xi)
    u <- stats::rnorm(n = p, mean = 0, sd = sqrt(diag_D))
    f <- stats::rnorm(n = N, mean = 0, sd = 1)
    v <- X %*% u + f
    diag_D[-active_index] <- 0
    U <- diag_D * t(X)
    if (S < N) {
      yv <- y/sqrt(sigma2) - v
      Xyv <- t(X_s) %*% yv
      m <- solve(Q_star, Xyv)
      m_star <- yv - X_s %*% m
      new_beta <- sqrt(sigma2) * (u + U %*% m_star)
    } else {
      v_star <- solve(M, (y/sqrt(sigma2) - v))
      new_beta <- sqrt(sigma2) * (u + U %*% v_star)
    }
    # eta update
    eta <- rejection_sampler((new_beta^2)*xi/(2*sigma2), a, b)
    eta <- ifelse(eta <= 2.220446e-16, 2.220446e-16, eta)
    # save results
    betaout[i, ] <- new_beta
    etaout[i, ] <- eta
    xiout[i] <- xi
    sigma2out[i] <- sigma2
    activeout[i, active_index] <- 1
    # when use a auto threshold
    if (auto.threshold == TRUE) {
      if (i %% t == 0) {
        u_i <- stats::runif(1, 0, 1)
        p_i <- exp(adapt_p0 + adapt_p1 * i)
        if (u_i < p_i) {
          m_eff <- sum(1/((eta*xi)/s2.vec + 1))
        }
      }
      threshold <- sort(eta)[ceiling(m_eff)]
      active_index <- which(eta <= threshold)
      S <- length(active_index)
    }
  }
  betaout <- betaout[(burn+1):nmc, , drop = FALSE]
  lambdaout <- 1/sqrt(etaout[(burn+1):nmc, , drop = FALSE])
  activeout <- activeout[(burn+1):nmc, , drop = FALSE]
  tauout <- 1/sqrt(xiout[(burn+1):nmc])
  sigma2out <- sigma2out[(burn+1):nmc]
  betahat <- apply(betaout, 2, mean)
  lambdahat <- apply(lambdaout, 2, mean)
  tauhat <- mean(tauout)
  sigma2hat <- mean(sigma2out)
  activemean <- sum(activeout)/nrow(activeout)
  leftci <- apply(betaout, 2, stats::quantile, probs = alpha/2)
  rightci <- apply(betaout, 2, stats::quantile, probs = 1-alpha/2)
  result <- list(BetaHat = betahat, LeftCI = leftci, RightCI = rightci,
                 Sigma2Hat = sigma2hat, TauHat = tauhat, LambdaHat = lambdahat,
                 ActiveMean = activemean, BetaSamples = betaout,
                 LambdaSamples = lambdaout, TauSamples = tauout,
                 Sigma2Samples = sigma2out, ActiveSet = activeout)
  return(result)
}

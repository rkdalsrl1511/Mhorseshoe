#' Run exact MCMC algorithm for horseshoe prior
#'
#' The exact MCMC algorithm introduced in Section 2.1 of Johndrow et al. (2020)
#' was implemented in this function. This algorithm is the horseshoe estimator
#' that updates the global shrinkage parameter \eqn{\tau} using
#' Metropolis-Hastings algorithm, and uses blocked-Gibbs sampling for
#' \eqn{(\tau, \beta, \sigma)}. The local shrinkage parameter
#' \eqn{\lambda_{j},\ j = 1,2,...,p} is updated by the rejection sampler.
#'
#' See \code{\link{Mhorseshoe}} or browseVignettes("Mhorseshoe").
#'
#' @references Johndrow, J., Orenstein, P., & Bhattacharya, A. (2020). Scalable
#' Approximate MCMC Algorithms for the Horseshoe Prior. In Journal of Machine
#' Learning Research (Vol. 21).
#'
#' @param X Design matrix, \eqn{X \in \mathbb{R}^{N \times p}}.
#' @param y Response vector, \eqn{y \in \mathbb{R}^{N}}.
#' @param burn Number of burn-in samples. Default is 1000.
#' @param iter Number of samples to be drawn from the posterior. Default is
#'  5000.
#' @param a Parameter of the rejection sampler, and it is recommended to leave
#'  it at the default value, \eqn{a = 1/5}.
#' @param b Parameter of the rejection sampler, and it is recommended to leave
#'  it at the default value, \eqn{b = 10}.
#' @param s \eqn{s^{2}} is the variance of tau's MH proposal distribution.
#'  0.8 is a good default. If set to 0, the algorithm proceeds by
#'  fixing the global shrinkage parameter \eqn{\tau} to the initial setting
#'  value.
#' @param tau Initial value of the global shrinkage parameter \eqn{\tau} when
#'  starting the algorithm. Default is 1.
#' @param sigma2 error variance \eqn{\sigma^{2}}. Default is 1.
#' @param w Parameter of gamma prior for \eqn{\sigma^{2}}. Default is 1.
#' @param alpha \eqn{100(1-\alpha)\%} credible interval setting argument.
#' @return \item{BetaHat}{Posterior mean of \eqn{\beta}.}
#' \item{LeftCI}{Lower bound of \eqn{100(1-\alpha)\%} credible interval for
#'  \eqn{\beta}.}
#' \item{RightCI}{Upper bound of \eqn{100(1-\alpha)\%} credible interval for
#'  \eqn{\beta}.}
#' \item{Sigma2Hat}{Posterior mean of \eqn{\sigma^{2}}.}
#' \item{TauHat}{Posterior mean of \eqn{\tau}.}
#' \item{LambdaHat}{Posterior mean of \eqn{\lambda_{j},\ j=1,2,...p.}.}
#' \item{BetaSamples}{Samples from the posterior of \eqn{\beta}.}
#' \item{LambdaSamples}{Lambda samples through rejection sampling.}
#' \item{TauSamples}{Tau samples through MH algorithm.}
#' \item{Sigma2Samples}{Samples from the posterior of the parameter
#'  \eqn{sigma^{2}}.}
#'
#' @examples
#' # Making simulation data.
#' set.seed(123)
#' N <- 50
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
#' # Run
#' result <- exact_horseshoe(X, y, burn = 0, iter = 100)
#'
#' # posterior mean
#' betahat <- result$BetaHat
#'
#' # Lower bound of the 95% credible interval
#' leftCI <- result$LeftCI
#'
#' # Upper bound of the 95% credible interval
#' RightCI <- result$RightCI
#'
#' @export
exact_horseshoe <- function(X, y, burn = 1000, iter = 5000, a = 1/5, b = 10,
                            s = 0.8, tau = 1, sigma2 = 1, w = 1, alpha = 0.05) {
  N <- nrow(X)
  p <- ncol(X)
  eta <- rep(1, p)
  xi <- tau^(-2)
  Q <- t(X) %*% X
  nmc <- burn + iter
  betaout <- matrix(0, nrow = nmc, ncol = p)
  etaout <- matrix(0, nrow = nmc, ncol = p)
  xiout <- rep(0, nmc)
  sigma2out <- rep(0, nmc)
  # run
  for(i in 1:nmc) {
    log_xi <- stats::rnorm(1, mean = log(xi), sd = sqrt(s))
    new_xi <- exp(log_xi)
    # xi update
    if (p < N) {
      Xy <- t(X) %*% y
      y_square <- t(y) %*% y
      Q_star <- xi * diag(eta) + Q
      m <- solve(Q_star, Xy)
      ymy <- y_square - t(y) %*% X %*% m
      if (s != 0) {
        new_Q_star <- new_xi * diag(eta) + Q
        new_m <- solve(new_Q_star, Xy)
        new_ymy <- y_square - t(y) %*% X %*% new_m
        cM <- (diag(chol(Q_star))^2)/xi
        new_cM <- (diag(chol(new_Q_star))^2)/new_xi
        curr_ratio <- -sum(log(cM))/2 - ((N + w)/2) * log(w+ymy) -
          log(sqrt(xi) * (1 + xi))
        new_ratio <- -sum(log(new_cM))/2 - ((N + w)/2) * log(w + new_ymy) -
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
      DX <- (1/eta) * t(X)
      XDX <- X %*% DX
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
    diag_D <- 1/(eta * xi)
    u <- stats::rnorm(n = p, mean = 0, sd = sqrt(diag_D))
    f <- stats::rnorm(n = N, mean = 0, sd = 1)
    v <- X %*% u + f
    U <- diag_D * t(X)
    if (p < N) {
      yv <- y/sqrt(sigma2) - v
      Xyv <- t(X) %*% yv
      m <- solve(Q_star, Xyv)
      m_star <- yv - X %*% m
      new_beta <- sqrt(sigma2) * (u + U %*% m_star)
    } else {
      v_star <- solve(M, (y / sqrt(sigma2) - v))
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
  }
  betaout <- betaout[(burn+1):nmc, ]
  lambdaout <- 1/sqrt(etaout[(burn+1):nmc, ])
  tauout <- 1/sqrt(xiout[(burn+1):nmc])
  sigma2out <- sigma2out[(burn+1):nmc]
  betahat <- apply(betaout, 2, mean)
  lambdahat <- apply(lambdaout, 2, mean)
  tauhat <- mean(tauout)
  sigma2hat <- mean(sigma2out)
  leftci <- apply(betaout, 2, stats::quantile, probs = alpha/2)
  rightci <- apply(betaout, 2, stats::quantile, probs = 1-alpha/2)
  result <- list(BetaHat = betahat, LeftCI = leftci, RightCI = rightci,
                 Sigma2Hat = sigma2hat, TauHat = tauhat, LambdaHat = lambdahat,
                 BetaSamples = betaout, LambdaSamples = lambdaout,
                 TauSamples = tauout, Sigma2Samples = sigma2out)
  return(result)
}

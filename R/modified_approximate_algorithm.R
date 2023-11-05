#' Run modified approximate MCMC algorithm for horseshoe prior
#'
#' This function implements an approximate MCMC algorithm that reduces
#' computational cost in the same way as Johndrow et al. (2020). See
#' \code{\link{approx_horseshoe}}. Unlike approx_horseshoe, instead of using a
#' fixed threshold, it uses an adaptive probability method to update the
#' threshold estimated by applying the shrinkage weight of Piironen and Vehtari
#' (2017).
#'
#' In this algorithm, the process of updating a new threshold is added by
#' applying the properties of the shrinkage weight \eqn{k_{j},\ j=1,2,...,p}
#' proposed by Piironen and Vehtari (2017). In the prior of \eqn{\beta_{j} \sim
#' N(0, \sigma^{2}\tau^{2}\lambda_{j}^{2}) = N(0, \sigma^{2}\xi^{-1}
#' \eta_{j}^{-1})}, the variable \eqn{m_{eff}} is defined as follows.
#'
#' \deqn{k_{j} = 1/\left(1+n\xi^{-1}s_{j}^{2}\eta_{j}^{-1} \right),
#' \quad j=1,2,...,p, \\ m_{eff} = \sum_{j=1}^{p}{\left(1-k_{j} \right)}.}
#'
#' The assumptions and notations for the model are the same as those in
#' \code{\link{Mhorseshoe}}, and \eqn{s_{j},\ j=1,2,...,p} are the diagonal
#' components of \eqn{X^{T}X}. For the zero components of \eqn{\beta},
#' \eqn{k_{j}} is derived close to 1, and nonzero's \eqn{k_{j}} is derived
#' close to 0, so the variable \eqn{m_{eff}} is called the effective number of
#' nonzero coefficients. In this algorithm, the threshold \eqn{\delta} is
#' updated to set \eqn{s_{\delta} = ceiling(m_{eff})}. \eqn{s_{\delta}} is the
#' same as \eqn{length(S)} in details of \code{\link{approx_horseshoe}}.
#'
#' Adaptive probability is defined to satisfy
#' Theorem 5(diminishing adaptation condition) of Roberts and Rosenthal (2007).
#' at \eqn{T}th iteration,
#'
#' \deqn{p(T) = exp[a_{0} + a_{1}T],\quad a_{1} < 0,
#' \quad u \sim U(0, 1), \\ if\ u < p(T),\ update\ \delta\ so\ that\
#' s_{\delta} = ceiling(m_{eff}).}
#'
#' The default is \eqn{a_{0} = 0}, \eqn{a_{1} = -4.6 \times 10^{-4}}, and
#' under this condition, \eqn{p(10000) < 0.01} is satisfied.
#'
#' @references Johndrow, J., Orenstein, P., & Bhattacharya, A. (2020).
#' Scalable Approximate MCMC Algorithms for the Horseshoe Prior. In Journal
#' of Machine Learning Research (Vol. 21).
#'
#' Piironen, J., & Vehtari, A. (2017). Sparsity information and regularization
#' in the horseshoe and other shrinkage priors. Electronic Journal of
#' Statistics, 11, 5018-5051.
#'
#' Roberts G, Rosenthal J. Coupling and ergodicity of adaptive Markov chain
#' Monte Carlo algorithms. J Appl Prob. 2007;44:458–475.
#'
#' @inheritParams exact_horseshoe
#' @param t Threshold \eqn{\delta} update cycle for adaptive probability
#'  method. default is 10.
#' @param alpha0 Parameter \eqn{a_{0}} of adaptive probability,
#'  \eqn{p(t) = exp[a_{0} + a_{1}t]}. default is 0.
#' @param alpha1 Parameter \eqn{a_{1}} of adaptive probability,
#'  \eqn{p(t) = exp[a_{0} + a_{1}t]}. default is \eqn{-4.6 \times 10^{-4}}.
#' @return \item{BetaSamples}{Samples from the posterior of the parameter
#' \eqn{\beta}.}
#' \item{LambdaSamples}{Lambda samples through rejection sampling.}
#' \item{TauSamples}{Tau samples through MH algorithm.}
#' \item{Sigma2Samples}{Samples from the posterior of the parameter
#' \eqn{sigma^{2}}.}
#' \item{ActiveSamples}{Samples recording the number of columns that satisfy
#' condition \eqn{\tau^{2}\lambda^{2}_{j} > \delta,\ j=1,2,...,p} in this
#' algorithm.}
#'
#' @examples
#' \dontrun{
#' # Making simulation data.
#' set.seed(123)
#' N <- 100
#' p <- 200
#' true_beta <- c(rep(1, 10), rep(0, 190))
#'
#' X <- matrix(1, nrow = N, ncol = p) # Design matrix X.
#' for (i in 1:p) {
#'   X[, i] <- stats::rnorm(N, mean = 0, sd = 1)
#' }
#'
#' y <- vector(mode = "numeric", length = N) # Response variable y.
#' e <- rnorm(N, mean = 0, sd = 2) # error term e.
#' for (i in 1:50) {
#'   y <- y + true_beta[i] * X[, i]
#' }
#' y <- y + e
#'
#' # Run as default
#' result <- mapprox_horseshoe(X, y, iteration = 1000)
#'
#' # posterior mean
#' post_mean <- apply(result$BetaSamples, MARGIN = 2, mean)
#'
#' # Lower bound of the 95% credible interval
#' post_leftCI <- apply(result$BetaSamples, MARGIN = 2,
#'                      quantile, probs = 0.025)
#'
#' # Upper bound of the 95% credible interval
#' post_rightCI <- apply(result$BetaSamples, MARGIN = 2,
#'                       quantile, probs = 0.975)}
#'
#' @export
mapprox_horseshoe <- function(X, y, iteration = 5000, a = 1/5, b = 10,
                              s = 0.8, tau = 1, sigma2 = 1, w = 1,
                              t = 10, alpha0 = 0, alpha1 = -4.6*10^(-4)) {
  N <- nrow(X)
  p <- ncol(X)
  S <- p
  active_index <- 1:p
  m_eff <- p
  eta <- rep(1, p)
  xi <- tau^(-2)
  Q <- t(X) %*% X
  s2.vec <- diag(Q)
  betaout <- matrix(0, nrow = iteration, ncol = p)
  etaout <- matrix(0, nrow = iteration, ncol = p)
  xiout <- rep(0, iteration)
  sigma2out <- rep(0, iteration)
  activeout <- rep(0, iteration)
  for(i in 1:iteration) {
    Q_s <- Q[active_index, active_index, drop = FALSE]
    X_s <- X[, active_index, drop = FALSE]
    log_xi <- stats::rnorm(1, mean = log(xi), sd = s)
    new_xi <- exp(log_xi)
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
      DX <- (1 / eta[active_index]) * t(X_s)
      XDX <- X_s %*% DX
      M <- diag(N) + XDX / xi
      m <- solve(M, y)
      ymy <- t(y) %*% m
      if (s != 0) {
        new_M <- diag(N) + XDX / new_xi
        new_m <- solve(new_M, y)
        new_ymy <- t(y) %*% new_m
        cM <- diag(chol(M))^2
        new_cM <- diag(chol(new_M))^2
        curr_ratio <- -sum(log(cM)) / 2 - ((N + w) / 2) * log(w + ymy) -
          log(sqrt(xi) * (1 + xi))
        new_ratio <- -sum(log(new_cM)) / 2 - ((N + w) / 2) * log(w+new_ymy) -
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
    sigma2 <- 1 / stats::rgamma(1, shape = (w + N) / 2, rate = (w + ymy) / 2)
    diag_D <- 1 / (eta * xi)
    u <- stats::rnorm(n = p, mean = 0, sd = sqrt(diag_D))
    f <- stats::rnorm(n = N, mean = 0, sd = 1)
    v <- X %*% u + f
    diag_D[-active_index] <- 0
    U <- diag_D * t(X)
    if (S < N) {
      yv <- y / sqrt(sigma2) - v
      Xyv <- t(X_s) %*% yv
      m <- solve(Q_star, Xyv)
      m_star <- yv - X_s %*% m
      new_beta <- sqrt(sigma2) * (u + U %*% m_star)
    } else {
      v_star <- solve(M, (y / sqrt(sigma2) - v))
      new_beta <- sqrt(sigma2) * (u + U %*% v_star)
    }
    betaout[i, ] <- new_beta
    etaout[i, ] <- eta
    xiout[i] <- xi
    sigma2out[i] <- sigma2
    activeout[i] <- S
    eta <- rejection_sampler((new_beta^2) * xi / (2 * sigma2), a, b)
    eta <- ifelse(eta <= 2.220446e-16, 2.220446e-16, eta)
    if (i %% t == 0) {
      u_i <- stats::runif(1, 0, 1)
      p_i <- exp(alpha0 + alpha1 * i)
      if (u_i < p_i)
        m_eff <- sum(1 / ((eta*xi) / s2.vec + 1))
    }
    threshold <- sort(eta)[ceiling(m_eff)]
    active_index <- which(eta <= threshold)
    S <- length(active_index)
  }
  result <- list(BetaSamples = betaout,
                 LambdaSamples = 1 / sqrt(etaout),
                 TauSamples = 1 / sqrt(xiout),
                 Sigma2Samples = sigma2out,
                 ActiveSamples = activeout)
  return(result)
}

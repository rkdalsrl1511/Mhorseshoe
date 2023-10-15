#' Run approximate MCMC algorithm for horseshoe prior
#'
#' In this function, The algorithm introduced in Section 2.2 of Johndrow et al.
#' (2020) is implemented, and is a horseshoe estimator that generally considers
#' the case where \eqn{p>>N}. The assumptions and notations for the model are
#' the same as those in \code{\link{Mhorseshoe}}. This algorithm introduces a
#' threshold and uses only a portion of the total \eqn{p} columns for matrix
#' multiplication, lowering the computational cost compared to the existing
#' horseshoe estimator. According to Section 3.2 of Johndrow et al. (2020), the
#' approximate MCMC algorithm applying the methodology constructs an
#' approximate Markov chain \eqn{P_{\epsilon}} that can converge to an exact
#' Markov chain \eqn{P}, and acceptable results were confirmed through
#' empirical analysis of simulation and real data.
#'
#' This algorithm has the following changes compared to the exact algorithm:
#'
#' \deqn{D_{\delta} = diag\left(\eta_{j}^{-1}1\left(\xi^{-1}\eta_{j}^{-1}
#' > \delta,\ j=1,2,...,p. \right) \right),}
#' \deqn{M_{\xi} \approx M_{\xi, \delta} = I_{N} + \xi^{-1}XD_{\delta}X^{T}.}
#'
#' The set of columns that satisfies the condition(\eqn{\xi^{-1}\eta_{j}^{-1}
#' > \delta}) is defined as the active set, and let's define \eqn{S} as the
#' index set of the following columns.
#'
#' \deqn{S = \{j\ |\ \xi^{-1}\eta_{j}^{-1} > \delta,\ j=1,2,...,p. \}.}
#'
#' If \eqn{\xi^{-1}\eta_{j}^{-1}} is very small, the posterior of \eqn{\beta}
#' will have a mean and variance close to 0. Therefore, let's set
#' \eqn{\xi^{-1}\eta_{j}^{-1}} smaller than \eqn{\delta} to 0 and the size of
#' inverse \eqn{M_{\xi, \delta}} matrix is reduced as follows.
#'
#' \deqn{length(S)=s_{\delta} \le p, \\ X_{S} \in R^{N \times s_{\delta}},
#' \quad D_{S} \in R^{s_{\delta} \times s_{\delta}}, \\ M_{\xi, \delta}^{-1} =
#' \left(I_{N} + \xi^{-1}X_{S}D_{S}X_{S}^{T} \right)^{-1}.}
#'
#' \eqn{M_{\xi, \delta}^{-1}} can be expressed using the Woodbury identity
#' as follows.
#'
#' \deqn{M_{\xi, \delta}^{-1} = I_{N} - X_{S}\left(\xi D_{S}^{-1} +
#' X_{S}^{T}X_{S} \right)^{-1}X_{S}^{T}.}
#'
#' \eqn{M_{\xi, \delta}^{-1}}, which reduces the computational cost, is
#' applied to all parts of this algorithm, \eqn{\beta} samples are extracted
#' from the posterior using fast sampling(Bhattacharya et al.,2016) as follows.
#'
#' \deqn{u \sim N_{p}(0, \xi^{-1}D),\quad f \sim N_{N}(0, I_{N}), \\
#' v = Xu + f,\quad v^{\star} = M_{\xi, \delta}^{-1}(y/\sigma - v), \\
#' \beta = \sigma(u + \xi^{-1}D_{\delta}X^{T}v^{\star}).}
#'
#' @references Bhattacharya, A., Chakraborty, A., & Mallick, B. K. (2016).
#' Fast sampling with Gaussian scale mixture priors in high-dimensional
#' regression. Biometrika, asw042.
#'
#' Johndrow, J., Orenstein, P., & Bhattacharya, A. (2020).
#' Scalable Approximate MCMC Algorithms for the Horseshoe Prior. In Journal
#' of Machine Learning Research (Vol. 21).
#'
#' @inheritParams exact_horseshoe
#' @param threshold Threshold \eqn{\delta} to be used in the approximate MCMC
#'  algorithm. If you select NULL(this is the default), the default is set
#'  according to the sizes of N and p. if \eqn{p < N},
#'  \eqn{\delta = 1/\sqrt{Np}}, else \eqn{\delta = 1/p}. Or, you can set the
#'  value directly through this argument. For more information about
#'  \eqn{\delta}, see \code{\link{Mhorseshoe}} and 4.1 of Johndrow et al.
#'  (2020).
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
#' result <- approx_horseshoe(X, y, iteration = 1000)
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
approx_horseshoe <- function(X, y, iteration = 5000, threshold = NULL,
                             a = 1 / 5, b = 10, s = 0.8, tau = 1, sigma2 = 1,
                             w = 1) {
  N <- nrow(X)
  p <- ncol(X)
  if (is.null(threshold))
    threshold <- ifelse(p >= N, 1 / p, 1 / sqrt(N * p))
  eta <- rep(1, p)
  xi <- tau^(-2)
  Q <- t(X) %*% X
  betaout <- matrix(0, nrow = iteration, ncol = p)
  etaout <- matrix(0, nrow = iteration, ncol = p)
  xiout <- rep(0, iteration)
  sigma2out <- rep(0, iteration)
  activeout <- rep(0, iteration)
  for(i in 1:iteration) {
    log_xi <- stats::rnorm(1, mean = log(xi), sd = s)
    new_xi <- exp(log_xi)
    max_xi <- max(xi, new_xi)
    active_index <- which((1 / (eta * max_xi) > threshold))
    S <- length(active_index)
    Q_s <- Q[active_index, active_index, drop = FALSE]
    X_s <- X[, active_index, drop = FALSE]
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
        new_ratio <- -sum(log(new_cM)) / 2 - ((N + w) / 2)*log(w + new_ymy) -
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
    sigma2 <- 1/stats::rgamma(1, shape = (w + N) / 2, rate = (w + ymy) / 2)
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
    if ((i %% 1000) == 0)
      cat("Number of iterations : ", i, "\n")
  }
  result <- list(BetaSamples = betaout,
                 LambdaSamples = 1 / sqrt(etaout),
                 TauSamples = 1 / sqrt(xiout),
                 Sigma2Samples = sigma2out,
                 ActiveSamples = activeout)
  return(result)
}

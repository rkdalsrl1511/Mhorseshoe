#' @description Mhorseshoe is a package for a high-dimensional Bayesian linear
#' modeling algorithm using a horseshoe prior. This package provides two
#' different algorithm functions : \code{\link{exact_horseshoe}},
#' \code{\link{approx_horseshoe}}. approx_horseshoe is version that can lower
#' the computational cost than the existing horseshoe estimator through the
#' approximate MCMC algorithm in the case of \eqn{p >> N} for \eqn{p}
#' predictors and \eqn{N} observations. You can see examples of the use of the
#' two algorithms through the vignette, `browseVignettes("Mhorseshoe")`.
#'
#' @keywords internal
#' @useDynLib Mhorseshoe, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import stats
#' @docType package
#' @name Mhorseshoe
"_PACKAGE"

NULL

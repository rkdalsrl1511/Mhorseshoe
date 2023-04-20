#' Run modified approximate algorithm
#' @title modified approximate algorithm
#' @param W Design matrix
#' @param z Response variable
#' @param iteration Number of iterations
#' @param independence_sigma Determine whether to remove sigma in the beta prior distribution
#' @param sigma error term's variance
#' @param w sigma prior's hyperparameter
#' @param a,b rejection sampler parameters
#' @param fixed_xi Determine whether to fix global shrinkage parameter
#' @param xi Global shrinkage parameter
#' @param s Variance of log xi when sampling xi
#' @param t Interval to run the adaptive probability
#' @param alpha0 Parameters of the adaptive probability
#' @param alpha1 Parameters of the adaptive probability
#' @param step_check Time measurement at each step of the algorithm
#' @return beta posterior samples
#' @export
modified_approximate_algorithm <- function(W, z, iteration = 5000,
                                           independence_sigma = FALSE, sigma = 1, w = 0,
                                           a = 1/5, b = 10, fixed_xi = FALSE, xi = 1, s=0.01,
                                           t = 50, alpha0 = -0.5, alpha1 = -7*10^(-4),
                                           step_check = FALSE) {


  if (independence_sigma == TRUE) {

    result <- modified_approximate_algorithm2(W, z, xi, sigma, iteration, a, b, w, t,
                                              alpha0, alpha1, step_check)

  } else if (fixed_xi == TRUE) {

    result <- modified_approximate_algorithm3(W, z, xi, sigma, iteration, a, b, w, t,
                                              alpha0, alpha1, step_check)

  } else {

    result <- modified_approximate_algorithm1(W, z, xi, sigma, iteration, a, b, s, w, t,
                                              alpha0, alpha1, step_check)

  }

  return(result)

}

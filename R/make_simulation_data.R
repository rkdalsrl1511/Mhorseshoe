#' A response variable of the simulation data is generated,
#' and a value of non-zero coefficients is 10.
#' @title make_response
#' @param W W matrix
#' @param non_zero Number of non_zero signals
#' @export
make_simulation_data <- function(N = 500, p = 1000, W_mean = 0, W_sd = 1, W_scaling = NULL,
                                 non_zero = 30, error_sd = 1, fixed_coefficients = NULL) {

  W <- matrix(1, nrow = N, ncol = p)



  if (is.null(W_scaling) == TRUE) {

    for (i in 1:p)
      W[, i] <- stats::rnorm(N, mean = W_mean, sd = W_sd)

  } else if (W_scaling == "min_max") {

    for (i in 1:p) {

      W[, i] <- stats::rnorm(N, mean = W_mean, sd = W_sd)
      W[, i] <- (W[, i] - min(W[, i]))/(max(W[, i]) - min(W[, i]))

    }

  } else if (W_scaling == "standard") {

    for (i in 1:p) {

      W[, i] <- stats::rnorm(N, mean = W_mean, sd = W_sd)
      W[, i] <- scale(W[, i], center = TRUE, scale = TRUE)

    }

  } else {

    stop("please enter 'min_max', 'standard' or NULL \n")

  }

  # non_zero column을 랜덤으로 할당
  non_zero_index <- sample(1:ncol(W), size = non_zero, replace = FALSE)
  non_zero_coefficients <- vector(mode = "numeric", length = non_zero)

  # z : response variable
  z <- vector(mode = "numeric", length = nrow(W))
  e <- rnorm(N, mean = 0, sd = error_sd)

  # non_zero coefficients을 -10~10 사이 값으로 설정
  for(i in 1:length(non_zero_index)) {

    if (is.null(fixed_coefficients) == TRUE) {

      coef <- round(runif(1,-10,10), digits = 0)

    } else if (length(fixed_coefficients) == 1) {

      coef <- fixed_coefficients

    } else {

      coef <- fixed_coefficients[i]

    }

    non_zero_coefficients[i] <- coef
    z <- z + coef * W[, non_zero_index[i]]

  }

  z <- z + e

  # 정보 출력
  for (i in 1:length(non_zero_index)) {

    cat("non zero coefficient : ", non_zero_index[i],
        "value : ", non_zero_coefficients[i], "\n")

  }

  # 반환
  return(list(W, z, non_zero_index, non_zero_coefficients))

}

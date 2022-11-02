#' make simulation data for Scalable Approximate MCMC Algorithms for the Horseshoe Prior
#' @title make_W
#' @param N Number of rows
#' @param p Number of cols
#' @export
make_W <- function(N = 100, p = 10, cov_matrix = NULL,
                   level = vector(mode = "numeric", length = 0)) {

  # numerical data
  if (is.null(cov_matrix)) {

    cov_matrix <- diag(p)

  }

  W <- MASS::mvrnorm(n = N, mu = rep(0, p), Sigma = cov_matrix)
  W <- cbind(rep(1, N), W)

  # categorical data
  if (length(level) > 0) {

    # 만약 level = c(2,3)이라면 각각 1개, 2개의 더미변수를 생성하고 W와 병합.
    for (i in 1:length(level)) {

      C <- matrix(0, nrow = N, ncol = level[i])

      for (j in 1:N) {

        u <- sample(1:level[i], 1, prob = rep(1/level[i], level[i]))
        C[j, u] <- 1
      }

      W <- cbind(W, C[, 1:(level[i]-1)])

    }

  }

  W

}


#' A response variable of the simulation data is generated,
#' and a value of non-zero coefficients is 10.
#' @title make_response
#' @param W W matrix
#' @param non_zero Number of non_zero signals
#' @export
make_response <- function(W, non_zero = 2,
                          SD = 5, Heteroskedasticity = FALSE,
                          fixed_coefficients = NULL) {

  # error
  N <- nrow(W)
  if (Heteroskedasticity == TRUE) {

    SD <- runif(n = N, min = 1, max = SD)
    e <- rnorm(n = N, mean = 0, sd = SD)

  } else {

    e <- rnorm(n = N, mean = 0, sd = SD)

  }

  # non_zero column을 랜덤으로 할당
  non_zero_index <- sample(1:ncol(W), size = non_zero, replace = FALSE)
  non_zero_coefficients <- vector(mode = "numeric", length = non_zero)

  # z : response variable
  z <- vector(mode = "numeric", length = nrow(W))

  # non_zero coefficients을 -100~100 사이 값으로 설정
  for(i in 1:length(non_zero_index)) {

    if (is.null(fixed_coefficients) == TRUE) {

      coef <- round(runif(1,-100,100), digits = 0)

    } else if (length(fixed_coefficients) == 1) {

      coef <- fixed_coefficients

    } else if (is.na(fixed_coefficients[i]) == TRUE){

      stop("올바른 형식을 입력하라.")

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
  if (Heteroskedasticity == TRUE) {

    return(list(z, non_zero_index, non_zero_coefficients, SD))

  } else {

    return(list(z, non_zero_index, non_zero_coefficients))

  }

}

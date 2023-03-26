# f(x) := f_e(x) 계산 함수
f_e_x <- function(Epsilon, x){

  Epsilon * x + log(1+x)

}

# lambda 계산 함수
lambda_calculation <- function(Epsilon, x1, x2){

  (f_e_x(Epsilon, x2/Epsilon) - f_e_x(Epsilon, x1/Epsilon))/(x2/Epsilon - x1/Epsilon)

}

# v 계산 함수
v_calculation <- function(Epsilon, lambda2, lambda3, a, b){

  A <- f_e_x(Epsilon, a/Epsilon)
  B <- f_e_x(Epsilon, b/Epsilon)
  I <- f_e_x(Epsilon, 1/Epsilon)

  v1 <- log(1+a/Epsilon)
  v2 <- (1/lambda2) * exp(-A) * (1 - exp(-(I-A)))
  v3 <- (1/lambda3) * exp(-I) * (1 - exp(-(B-I)))
  v4 <- (1/Epsilon) * exp(-B)

  v <- cbind(v1, v2, v3, v4)

  v

}

# f_L(z) 계산 함수
f_L_z <- function(selected_interval, Epsilon, z, a, b, lambda2, lambda3){

  result <- ifelse(selected_interval == 1, log(1+z),
                   ifelse(selected_interval == 2,
                          f_e_x(Epsilon, a/Epsilon) + lambda2 * (z - a/Epsilon),
                          ifelse(selected_interval == 3,
                                 f_e_x(Epsilon, 1/Epsilon) + lambda3 * (z - b/Epsilon), f_e_x(Epsilon, b/Epsilon) + Epsilon * (z - b/Epsilon))))

  result

}

# cdf값인 u의 값에 따라서 z를 설정해주는 함수
sample_z <- function(selected_interval, Epsilon,
                     lambda2, lambda3, u, a, b){

  z <- ifelse(selected_interval == 1, (1+a/Epsilon)^u - 1,
              ifelse(selected_interval == 2,
                     a/Epsilon - log(1 - u + u * exp(-(f_e_x(Epsilon, 1/Epsilon) - f_e_x(Epsilon, a/Epsilon))))/lambda2,
                     ifelse(selected_interval == 3,
                            1/Epsilon - log(1 - u + u * exp(-(f_e_x(Epsilon, b/Epsilon) - f_e_x(Epsilon, 1/Epsilon))))/lambda3,
                            b/Epsilon - log(1 - u)/Epsilon)))

  z

}

sample_eta <- function(eta, rejected_index, prob_v,
                       Epsilon, lambda2, lambda3, a, b) {

  p <- length(rejected_index)

  if (p == 0) {

    return(eta)

  } else {

    # random select interval
    u <- runif(p, min = 0, max = 1)
    selected_interval <- ifelse(u < prob_v[rejected_index, 1], 1,
                                ifelse(u < prob_v[rejected_index, 1] + prob_v[rejected_index, 2], 2,
                                       ifelse(u < prob_v[rejected_index, 1] + prob_v[rejected_index, 2] + prob_v[rejected_index, 3], 3, 4)))

    # sample z with cdf method
    u_z <- runif(p, min = 0, max = 1)
    z <- sample_z(selected_interval, Epsilon[rejected_index],
                  lambda2[rejected_index], lambda3[rejected_index],
                  u_z, a, b)

    f_l <- f_L_z(selected_interval, Epsilon[rejected_index],
                 z, a, b, lambda2[rejected_index], lambda3[rejected_index])
    f <- f_e_x(Epsilon[rejected_index], z)

    # acceptance / rejection
    u <- runif(p, min = 0, max = 1)
    eta[rejected_index] <- ifelse(u < exp(-(f- f_l)), z, -1)
    rejected_index <- which(eta == -1 | eta == 0)

    return(sample_eta(eta, rejected_index, prob_v, Epsilon, lambda2, lambda3, a, b))

  }

}

sample_eta2 <- function(eta, rejected_index, prob_v,
                       Epsilon, lambda2, lambda3, a, b) {

  while(p != 0){

    # random select interval
    u <- runif(p, min = 0, max = 1)
    selected_interval <- ifelse(u < prob_v[rejected_index, 1], 1,
                                ifelse(u < prob_v[rejected_index, 1] + prob_v[rejected_index, 2], 2,
                                       ifelse(u < prob_v[rejected_index, 1] + prob_v[rejected_index, 2] + prob_v[rejected_index, 3], 3, 4)))

    # sample z with cdf method
    u_z <- runif(p, min = 0, max = 1)
    z <- sample_z(selected_interval, Epsilon[rejected_index],
                  lambda2[rejected_index], lambda3[rejected_index],
                  u_z, a, b)

    f_l <- f_L_z(selected_interval, Epsilon[rejected_index],
                 z, a, b, lambda2[rejected_index], lambda3[rejected_index])
    f <- f_e_x(Epsilon[rejected_index], z)

    # acceptance / rejection
    u <- runif(p, min = 0, max = 1)
    eta[rejected_index] <- ifelse(u < exp(-(f- f_l)), z, -1)
    rejected_index <- which(eta == -1 | eta == 0)
    p <- length(rejected_index)

  }

  return(eta)

}


# foreach를 이용한 병렬처리 : 이건 foreach로 병렬연산이 아니라 하나하나 해야하네
# 시작도 안했지만 별로 효율은 안 좋을 것 같다.
sample_eta3 <- function(eta, rejected_index, prob_v,
                        Epsilon, lambda2, lambda3, a, b) {

  while(p != 0){

    # random select interval
    u <- runif(p, min = 0, max = 1)
    selected_interval <- ifelse(u < prob_v[rejected_index, 1], 1,
                                ifelse(u < prob_v[rejected_index, 1] + prob_v[rejected_index, 2], 2,
                                       ifelse(u < prob_v[rejected_index, 1] + prob_v[rejected_index, 2] + prob_v[rejected_index, 3], 3, 4)))

    # sample z with cdf method
    u_z <- runif(p, min = 0, max = 1)
    z <- sample_z(selected_interval, Epsilon[rejected_index],
                  lambda2[rejected_index], lambda3[rejected_index],
                  u_z, a, b)

    f_l <- f_L_z(selected_interval, Epsilon[rejected_index],
                 z, a, b, lambda2[rejected_index], lambda3[rejected_index])
    f <- f_e_x(Epsilon[rejected_index], z)

    # acceptance / rejection
    u <- runif(p, min = 0, max = 1)
    eta[rejected_index] <- ifelse(u < exp(-(f- f_l)), z, -1)
    rejected_index <- which(eta == -1 | eta == 0)
    p <- length(rejected_index)

  }

  return(eta)

}

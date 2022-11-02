
# xi's acceptance probability
probability_a <- function(N,
                          xi, new_xi,
                          M, inverse_M, new_inverse_M,
                          z, w){

  a <- sqrt(det(new_inverse_M %*% M)) * ((w/2 + t(z) %*% new_inverse_M %*% z / 2)/(w/2 + t(z) %*% inverse_M %*% z / 2))^(-(N+w)/2) * ((sqrt(xi)*(1+xi))/(sqrt(new_xi)*(1+new_xi))) * (new_xi/xi)

  if (a == "NaN") {

    a <- 0

  }

  a

}

# xi's acceptance probability
probability_a2 <- function(N,
                           xi, new_xi,
                           WW, new_WW,
                           inverse_M, new_inverse_M,
                           z, w){

  a <- sqrt(det(new_WW  %*% WW * new_xi / xi)) * ((w/2 + t(z) %*% new_inverse_M %*% z / 2)/(w/2 + t(z) %*% inverse_M %*% z / 2))^(-(N+w)/2) * ((sqrt(xi)*(1+xi))/(sqrt(new_xi)*(1+new_xi))) * (new_xi/xi)

  if (a == "NaN") {

    a <- 0

  }

  a

}

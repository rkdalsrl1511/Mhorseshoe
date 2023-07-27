
probability_a <- function(N, xi, new_xi, k, zmz, new_zmz, w){

  a <- k*(((w/2 + new_zmz/2)/(w/2 + zmz/2))^(-(N+w)/2))*((sqrt(xi)*(1+xi))/(sqrt(new_xi)*(1+new_xi)))*(new_xi/xi)

  if (is.nan(a)) {

    return(0)

  } else {

    return(a)

  }



}

lmh_latio <- function(N, xi, cM, zmz, w){

  a <- -sum(log(cM))/2 - ((N+w)/2)*log(w+zmz) - log(sqrt(xi)*(1+xi))

  return(a)

}

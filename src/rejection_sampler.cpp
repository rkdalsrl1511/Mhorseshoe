//' @useDynLib Mhorseshoe, .registration = TRUE
//' @importFrom Rcpp sourceCpp evalCpp
#include <Rcpp.h>
using namespace Rcpp;

NumericVector fex(NumericVector eps, NumericVector x){
  NumericVector f = eps*x+log(1+x);
  return f;
}

double dfex(double eps, double x){
  double df = eps*x+log(1+x);
  return df;
}

NumericVector calculate_lambda(NumericVector eps, double x, double y){
  NumericVector lambda = (fex(eps, y/eps)-fex(eps, x/eps))/(y/eps-x/eps);
  return lambda;
}

double sample_z(int selected_interval, double eps, double l2, double l3,
                double u, double a, double b){
  double z;
  switch(selected_interval){
  case 1:
    z = exp(u*log(1+a/eps))-1;
    break;
  case 2:
    z = a/eps-log(1-u+u*exp(-(dfex(eps, 1/eps)-dfex(eps, a/eps))))/l2;
    break;
  case 3:
    z = 1/eps-log(1-u+u*exp(-(dfex(eps, b/eps)-dfex(eps, 1/eps))))/l3;
    break;
  default:
    z = b/eps-log(1-u)/eps;
  }
  return z;
}

double fl(int selected_interval, double eps, double z, double a, double b,
          double l2, double l3){
  double f;
  switch(selected_interval){
  case 1:
    f = log(1+z);
    break;
  case 2:
    f = dfex(eps, a/eps)+l2*(z-a/eps);
    break;
  case 3:
    f = dfex(eps, 1/eps)+l3*(z-b/eps);
    break;
  default:
    f = dfex(eps, b/eps)+eps*(z-b/eps);
  }
  return f;
}

// [[Rcpp::export]]
NumericVector rejection_sampler(NumericVector eps, double a, double b){
  int p = eps.length();
  NumericVector lambda2 = calculate_lambda(eps, a, 1);
  NumericVector lambda3 = calculate_lambda(eps, 1, b);
  NumericVector A = fex(eps, a/eps);
  NumericVector I = fex(eps, 1/eps);
  NumericVector B = fex(eps, b/eps);
  NumericMatrix total_prob (p, 4);
  total_prob(_, 0) = log(1+a/eps);
  total_prob(_, 1) = total_prob(_, 0) + (1/lambda2) * exp(-A) * (1 - exp(-(I-A)));
  total_prob(_, 2) = total_prob(_, 1) + (1/lambda3) * exp(-I) * (1 - exp(-(B-I)));
  total_prob(_, 3) = total_prob(_, 2) + (1/eps) * exp(-B);
  total_prob(_, 0) = total_prob(_, 0) / total_prob(_, 3);
  total_prob(_, 1) = total_prob(_, 1) / total_prob(_, 3);
  total_prob(_, 2) = total_prob(_, 2) / total_prob(_, 3);
  std::fill(total_prob(_, 3).begin(), total_prob(_, 3).end(), 1);
  NumericVector eta(p);
  for(int i=0; i<p; ++i){
    if(eps[i] > 1){
      NumericVector v = Rcpp::runif(1, 0, 1);
      double z = -log(1-v[0])/eps[i];
      NumericVector u = Rcpp::runif(1, 0, 1);
      if(u[0] < (1/(1+z))){
        eta[i] = z;
      } else{
        --i;
      }
    } else {
      NumericVector u = Rcpp::runif(1, 0, 1);
      int selected_interval;
      if(u[0] < total_prob(i, 0)){
        selected_interval = 1;
      } else if(u[0] < total_prob(i, 1)){
        selected_interval = 2;
      } else if(u[0] < total_prob(i, 2)){
        selected_interval = 3;
      } else{
        selected_interval = 4;
      }
      NumericVector u_z = Rcpp::runif(1, 0, 1);
      double z = sample_z(selected_interval, eps[i], lambda2[i], lambda3[i], u_z[0], a, b);
      double f_l = fl(selected_interval, eps[i], z, a, b, lambda2[i], lambda3[i]);
      double f = dfex(eps[i], z);
      u = Rcpp::runif(1, 0, 1);
      if(u[0] < exp(-(f - f_l))){
        eta[i] = z;
      } else{
        --i;
      }
    }
  }
  return eta;
}

#include <Rcpp.h>
using namespace Rcpp;

//' rejection sampler to update local shrinkage parameters lambda.
//' @keywords internal
// [[Rcpp::export]]
NumericVector rejection_sampler(NumericVector eps, double a, double b){
  int p = eps.length();
  NumericVector a_eps = a/eps;
  NumericVector i_eps = 1/eps;
  NumericVector b_eps = b/eps;
  NumericVector A = a + log(1 + a_eps);
  NumericVector I = 1 + log(1 + i_eps);
  NumericVector B = b + log(1 + b_eps);
  NumericVector lambda2 = (I - A) / (i_eps - a_eps);
  NumericVector lambda3 = (B - I) / (b_eps - i_eps);
  NumericMatrix total_prob (p, 3);
  total_prob(_, 0) = log(1+a_eps);
  total_prob(_, 1) = total_prob(_, 0) + (1/lambda2) * exp(-A) * (1 - exp(-(I-A)));
  total_prob(_, 2) = total_prob(_, 1) + (1/lambda3) * exp(-I) * (1 - exp(-(B-I)));
  NumericVector total_volume = total_prob(_, 2) + (i_eps) * exp(-B);
  total_prob(_, 0) = total_prob(_, 0)/total_volume;
  total_prob(_, 1) = total_prob(_, 1)/total_volume;
  total_prob(_, 2) = total_prob(_, 2)/total_volume;
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
      NumericVector u_z = Rcpp::runif(1, 0, 1);
      double z;
      double fl;
      if(u[0] < total_prob(i, 0)){
        z = exp(u_z[0] * log(1 + a_eps[i]))-1;
        fl = log(1 + z);
      } else if(u[0] < total_prob(i, 1)){
        z = a_eps[i] - log(1 - u_z[0] + u_z[0]*exp(-(I[i]-A[i])))/lambda2[i];
        fl = A[i] + lambda2[i] * (z - a_eps[i]);
      } else if(u[0] < total_prob(i, 2)){
        z = i_eps[i] - log(1 - u_z[0] + u_z[0]*exp(-(B[i]-I[i])))/lambda3[i];
        fl = I[i] + lambda3[i] * (z - b_eps[i]);
      } else{
        z = b_eps[i] - log(1 - u_z[0])/eps[i];
        fl = B[i] + eps[i] * (z - b_eps[i]);
      }
      double f = eps[i] * z + log(1 + z);
      u = Rcpp::runif(1, 0, 1);
      if(u[0] < exp(-(f - fl))){
        eta[i] = z;
      } else{
        --i;
      }
    }
  }
  return eta;
}

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' estimate the regression estimates based on given the number of repetitions
//'
//' @param X matrix
//' @param y numeric vector
//' @param Z numeric vector for weight
// [[Rcpp::export]]
Rcpp::List fastlm(const arma::mat& X, const arma::colvec& y, const arma::colvec& Z) {
  int n = X.n_rows, k = X.n_cols;

  // slope coefficients
  // X^T * (W)^1/2 * (W)^1/2 * X = X^T * (W)^1/2 * (W)^1/2 * y
  arma::mat diag_Z    = diagmat(sqrt(Z));
  arma::mat Xasterisk = diag_Z * X;
  arma::mat Yasterisk = diag_Z * y;
  arma::colvec coef   = arma::solve(Xasterisk,Yasterisk);


  arma::colvec res    = y - X*coef; // residuals

  // std.errors of coefficients
  double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(n - k);
  arma::colvec std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X)*X)));


  return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
                            Rcpp::Named("rank")         = coef.n_rows,
                            Rcpp::Named("stderr")       = std_err,
                            Rcpp::Named("df.residual")  = n - k);
}
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
double getDeterminant(NumericMatrix AA){
  using Eigen::Map;
  using Eigen::MatrixXd;
  using Eigen::Lower;
  const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
  return  A.determinant();
}

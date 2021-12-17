// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::export(rng = false)]]
Eigen::SparseMatrix<double> IntegrateDataC(
    Eigen::SparseMatrix<double> integration_matrix,
    Eigen::SparseMatrix<double> weights,
    Eigen::SparseMatrix<double> expression_cells2
) {
  return(expression_cells2 - integration_matrix * weights);
}

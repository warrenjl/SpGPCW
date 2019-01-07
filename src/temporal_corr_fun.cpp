#include "RcppArmadillo.h"
#include "SpGPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List temporal_corr_fun(int m,
                             double phi){

double log_deter = 0.00; 
double sign = 0.00;     
arma::mat temporal_corr(m, m);
for(int j = 0; j < m; ++j){
   for(int k = 0; k < m; ++k){
      temporal_corr(j,k) = exp(-phi*abs(j - k));
      }
   }

arma::mat temporal_corr_inv = inv_sympd(temporal_corr);
log_det(log_deter, sign, temporal_corr);

return Rcpp::List::create(Rcpp::Named("temporal_corr_inv") = temporal_corr_inv,
                          Rcpp::Named("log_deter") = log_deter);

}


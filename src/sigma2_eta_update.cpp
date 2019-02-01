#include "RcppArmadillo.h"
#include "SpGPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_eta_update(arma::vec eta,
                         arma::mat corr_inv,
                         double alpha_sigma2_eta,
                         double beta_sigma2_eta){

int m = eta.size();
double alpha_sigma2_eta_update = 0.50*m + 
                                 alpha_sigma2_eta;

double beta_sigma2_eta_update = 0.50*dot(eta, ((corr_inv)*eta)) + 
                                beta_sigma2_eta;

double sigma2_eta = 1.00/R::rgamma(alpha_sigma2_eta_update,
                                   (1.00/beta_sigma2_eta_update));

return(sigma2_eta);

}






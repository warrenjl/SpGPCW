#include "RcppArmadillo.h"
#include "SpGPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_theta_update(arma::vec theta,
                           arma::mat corr_inv_theta,
                           double alpha_sigma2_theta,
                           double beta_sigma2_theta){

int m = theta.size();
double alpha_sigma2_theta_update = 0.50*m + 
                                   alpha_sigma2_theta;

double beta_sigma2_theta_update = 0.50*dot(theta, ((corr_inv_theta)*theta)) + 
                                  beta_sigma2_theta;

double sigma2_theta = 1.00/R::rgamma(alpha_sigma2_theta_update,
                                     (1.00/beta_sigma2_theta_update));

return(sigma2_theta);

}






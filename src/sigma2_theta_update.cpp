#include "RcppArmadillo.h"
#include "SpGPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_theta_update(arma::mat MCAR,
                           arma::mat theta,
                           double rho_old,
                           arma::vec eta_old,
                           arma::mat corr_inv,
                           double alpha_sigma2_theta,
                           double beta_sigma2_theta){
  
int s = theta.n_rows;
int m = theta.n_cols;
arma::vec theta_full(s*m); theta_full.fill(0.00);
arma::vec eta_full(s*m); eta_full.fill(0.00);
for(int j = 0; j < s; ++j){
   theta_full.subvec(m*j, ((j + 1)*m - 1)) = trans(theta.row(j));
   eta_full.subvec(m*j, ((j + 1)*m - 1)) = eta_old;
   }              
  
double alpha_sigma2_theta_update = 0.50*(s*m) + 
                                   alpha_sigma2_theta;

double beta_sigma2_theta_update = 0.50*dot((theta_full - eta_full), kron((rho_old*MCAR + (1.00 - rho_old)*eye(s, s)), corr_inv)*(theta_full - eta_full)) + 
                                  beta_sigma2_theta;

double sigma2_theta = 1.00/R::rgamma(alpha_sigma2_theta_update,
                                     (1.00/beta_sigma2_theta_update));

return(sigma2_theta);

}
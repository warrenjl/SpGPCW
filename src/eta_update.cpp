#include "RcppArmadillo.h"
#include "SpGPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec eta_update(arma::mat MCAR,
                     arma::mat z_star,
                     arma::mat theta,
                     double sigma2_theta,
                     double rho,
                     double sigma2_eta_old,
                     double phi_old,
                     arma::mat corr_inv){

arma::mat z_star_trans = trans(z_star);
int s = theta.n_rows;
int m = theta.n_cols;
arma::vec theta_full(s*m); theta_full.fill(0.00);
for(int j = 0; j < s; ++j){
  theta_full.subvec(m*j, ((j + 1)*m - 1)) = trans(theta.row(j));
  }

arma::mat cov_eta = inv_sympd((1.00/sigma2_eta_old)*corr_inv + 
                              (z_star_trans*(((kron((rho*MCAR + (1.00 - rho)*eye(s, s)), corr_inv))/sigma2_theta)*z_star)));

arma::vec mean_eta = cov_eta*(z_star_trans*(((kron((rho*MCAR + (1.00 - rho)*eye(s, s)), corr_inv))/sigma2_theta)*theta_full));

arma::mat ind_norms = arma::randn(1, m);
arma::vec eta = mean_eta + 
                trans(ind_norms*arma::chol(cov_eta));

return(eta);

}







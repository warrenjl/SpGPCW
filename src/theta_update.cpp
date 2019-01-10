#include "RcppArmadillo.h"
#include "SpGPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::vec theta_update(arma::mat x, 
                       arma::mat z,
                       arma::vec site_id,
                       arma::vec w,
                       arma::vec gamma,
                       arma::vec beta,
                       arma::mat eta_old,
                       double sigma2_theta_old,
                       arma::mat corr_inv){

int m = z.n_cols;
int n = w.size();
arma::mat w_mat(n, m);
for(int j = 0; j < m; ++j){
   w_mat.col(j) = w;
   }

arma::mat z_trans = trans(z);

arma::mat cov_theta = inv_sympd(z_trans*(w_mat%z) + 
                                (1.00/sigma2_theta_old)*corr_inv);

arma::vec mean_temp(n); mean_temp.fill(0.00);
int s = eta_old.n_rows;
for(int j = 0; j < s; ++j){
   arma::uvec ids = find(site_id == (j + 1));
   mean_temp.elem(ids) = z.rows(ids)*trans(eta_old.row(j));
   }

arma::vec mean_theta = cov_theta*(z_trans*(w%(gamma - x*beta - mean_temp)));

arma::mat ind_norms = arma::randn(1, m);
arma::vec theta = mean_theta + 
                  trans(ind_norms*arma::chol(cov_theta));

return(theta);

}







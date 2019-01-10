#include "RcppArmadillo.h"
#include "SpGPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

arma::mat eta_update(arma::mat eta_old,
                     arma::mat x, 
                     arma::mat z,
                     arma::vec site_id,
                     arma::mat neighbors,
                     arma::vec w,
                     arma::vec gamma,
                     arma::vec beta,
                     arma::vec theta,
                     double rho_old,
                     double sigma2_eta_old,
                     arma::mat corr_inv){
  
int n = w.size();
int m = z.n_cols;
arma::mat w_mat(n, m);
for(int j = 0; j < m; ++j){
   w_mat.col(j) = w;
   }

int s = eta_old.n_rows;
arma::mat z_trans = trans(z);
arma::mat eta = eta_old;
arma::vec mean_eta_temp(m); mean_eta_temp.fill(0.00);
for(int j = 0; j < s; ++j){
   mean_eta_temp.fill(0.00);
   for(int k = 0; k < s; ++k){
      mean_eta_temp = mean_eta_temp +
                      neighbors(j,k)*trans(eta.row(k));
      }
   arma::uvec ids = find(site_id == (j + 1));
   arma::mat cov_eta = inv_sympd(z_trans.cols(ids)*(w_mat.rows(ids)%z.rows(ids)) +
                                 ((rho_old*sum(neighbors.row(j)) + 1.00 - rho_old)/sigma2_eta_old)*corr_inv);
   arma::vec mean_eta = cov_eta*(z_trans.cols(ids)*(w.elem(ids)%(gamma.elem(ids) - x.rows(ids)*beta - z.rows(ids)*theta)) +
                                 ((rho_old/sigma2_eta_old)*corr_inv)*mean_eta_temp);
  
   arma::mat ind_norms = arma::randn(1, m);
   eta.row(j) = trans(mean_eta + 
                      trans(ind_norms*arma::chol(cov_eta)));

   }

//Centering for Stability (MCAR Style)
for(int j = 0; j < m; ++ j){
   eta.col(j) = eta.col(j) - 
                mean(eta.col(j));
   }

return(eta);

}







  


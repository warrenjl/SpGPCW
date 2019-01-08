#include "RcppArmadillo.h"
#include "SpGPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_eta_update(arma::mat neighbors,
                         arma::mat eta,
                         double rho,
                         arma::mat corr_inv_eta,
                         double alpha_sigma2_eta,
                         double beta_sigma2_eta){
  
int s = eta.n_rows;
int m = eta.n_cols;
arma::vec eta_full(s*m); eta_full.fill(0.00);
arma::vec diag_neighbors(s); diag_neighbors.fill(0.00);
for(int j = 0; j < s; ++j){
   eta_full.subvec(m*j, ((j + 1)*m - 1)) = trans(eta.row(j));
   diag_neighbors(j) = sum(neighbors.row(j));
   }              
  
double alpha_sigma2_eta_update = 0.50*(s*m) + 
                                 alpha_sigma2_eta;

arma::mat MCAR = diagmat(diag_neighbors) - 
                 neighbors;
double beta_sigma2_eta_update = 0.50*dot(eta_full, kron((rho*MCAR + (1.00 - rho)*eye(s, s)), corr_inv_eta)*eta_full) + 
                                beta_sigma2_eta;

double sigma2_eta = 1.00/R::rgamma(alpha_sigma2_eta_update,
                                   (1.00/beta_sigma2_eta_update));

return(sigma2_eta);

}
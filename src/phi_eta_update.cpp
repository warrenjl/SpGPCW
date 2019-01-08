#include "RcppArmadillo.h"
#include "SpGPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List phi_eta_update(double phi_eta_old,
                          arma::mat neighbors,
                          arma::mat eta,
                          double rho,
                          double sigma2_eta,
                          Rcpp::List temporal_corr_info_eta,
                          double a_phi_eta,
                          double b_phi_eta,
                          double metrop_var_phi_eta_trans,
                          int acctot_phi_eta_trans){
  
int s = eta.n_rows;
int m = eta.n_cols;
arma::vec eta_full(s*m); eta_full.fill(0.00);
arma::vec diag_neighbors(s); diag_neighbors.fill(0.00);
for(int j = 0; j < s; ++j){
   eta_full.subvec(m*j, ((j + 1)*m - 1)) = trans(eta.row(j));
   diag_neighbors(j) = sum(neighbors.row(j));
   }
arma::mat MCAR = diagmat(diag_neighbors) - 
                 neighbors;

/*Second*/
Rcpp::List temporal_corr_info_eta_old = temporal_corr_info_eta;
arma::mat corr_inv_eta_old = temporal_corr_info_eta_old[0];
double log_deter_eta_old = temporal_corr_info_eta_old[1];
double phi_eta_trans_old = log((phi_eta_old - a_phi_eta)/(b_phi_eta - phi_eta_old));

double second = -0.50*s*log_deter_eta_old - 
                (1.00/sigma2_eta)*0.50*dot(eta_full, kron((rho*MCAR + (1.00 - rho)*eye(s, s)), corr_inv_eta_old)*eta_full) + 
                phi_eta_trans_old -
                2.00*log(1.00 + exp(phi_eta_trans_old));

/*First*/
double phi_eta_trans = R::rnorm(phi_eta_trans_old, 
                                sqrt(metrop_var_phi_eta_trans));
double phi_eta = (b_phi_eta*exp(phi_eta_trans) + a_phi_eta)/(exp(phi_eta_trans) + 1.00);
temporal_corr_info_eta = temporal_corr_fun(m, phi_eta);
arma::mat corr_inv_eta = temporal_corr_info_eta[0];
double log_deter_eta = temporal_corr_info_eta[1];

double first = -0.50*s*log_deter_eta - 
               (1.00/sigma2_eta)*0.50*dot(eta_full, kron((rho*MCAR + (1.00 - rho)*eye(s, s)), corr_inv_eta)*eta_full) +
               phi_eta_trans -
               2.00*log(1.00 + exp(phi_eta_trans));

/*Decision*/
double ratio = exp(first - second);   
int acc = 1;
if(ratio < R::runif(0.00, 1.00)){
  phi_eta = phi_eta_old;
  temporal_corr_info_eta = temporal_corr_info_eta_old;
  acc = 0;
  }
acctot_phi_eta_trans = acctot_phi_eta_trans + 
                       acc;

return Rcpp::List::create(Rcpp::Named("phi_eta") = phi_eta,
                          Rcpp::Named("acctot_phi_eta_trans") = acctot_phi_eta_trans,
                          Rcpp::Named("temporal_corr_info_eta") = temporal_corr_info_eta);

}




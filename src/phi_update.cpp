#include "RcppArmadillo.h"
#include "SpGPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List phi_update(double phi_old,
                      arma::mat neighbors,
                      arma::vec theta,
                      double sigma2_theta,
                      arma::mat eta,
                      double rho,
                      double sigma2_eta,
                      Rcpp::List temporal_corr_info,
                      double a_phi,
                      double b_phi,
                      double metrop_var_phi_trans,
                      int acctot_phi_trans){
  
int s = eta.n_rows;
int m = theta.size();
arma::vec eta_full(s*m); eta_full.fill(0.00);
arma::vec diag_neighbors(s); diag_neighbors.fill(0.00);
for(int j = 0; j < s; ++j){
   eta_full.subvec(m*j, ((j + 1)*m - 1)) = trans(eta.row(j));
   diag_neighbors(j) = sum(neighbors.row(j));
   }
arma::mat MCAR = diagmat(diag_neighbors) - 
                 neighbors;

/*Second*/
Rcpp::List temporal_corr_info_old = temporal_corr_info;
arma::mat corr_inv_old = temporal_corr_info_old[0];
double log_deter_old = temporal_corr_info_old[1];
double phi_trans_old = log((phi_old - a_phi)/(b_phi - phi_old));

double second = -0.50*(1 + s)*log_deter_old - 
                (1.00/sigma2_theta)*0.50*dot(theta, (corr_inv_old*theta)) -
                (1.00/sigma2_eta)*0.50*dot(eta_full, kron((rho*MCAR + (1.00 - rho)*eye(s, s)), corr_inv_old)*eta_full) + 
                phi_trans_old -
                2.00*log(1.00 + exp(phi_trans_old));

/*First*/
double phi_trans = R::rnorm(phi_trans_old, 
                            sqrt(metrop_var_phi_trans));
double phi = (b_phi*exp(phi_trans) + a_phi)/(exp(phi_trans) + 1.00);
temporal_corr_info = temporal_corr_fun(m, phi);
arma::mat corr_inv = temporal_corr_info[0];
double log_deter = temporal_corr_info[1];

double first = -0.50*(1 + s)*log_deter - 
               (1.00/sigma2_theta)*0.50*dot(theta, (corr_inv*theta)) -
               (1.00/sigma2_eta)*0.50*dot(eta_full, kron((rho*MCAR + (1.00 - rho)*eye(s, s)), corr_inv)*eta_full) +
               phi_trans -
               2.00*log(1.00 + exp(phi_trans));

/*Decision*/
double ratio = exp(first - second);   
int acc = 1;
if(ratio < R::runif(0.00, 1.00)){
  phi = phi_old;
  temporal_corr_info = temporal_corr_info_old;
  acc = 0;
  }
acctot_phi_trans = acctot_phi_trans + 
                   acc;

return Rcpp::List::create(Rcpp::Named("phi") = phi,
                          Rcpp::Named("acctot_phi_trans") = acctot_phi_trans,
                          Rcpp::Named("temporal_corr_info") = temporal_corr_info);

}




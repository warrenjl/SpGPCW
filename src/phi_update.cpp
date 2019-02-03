#include "RcppArmadillo.h"
#include "SpGPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List phi_update(double phi_old,
                      arma::mat MCAR,
                      arma::mat theta,
                      double sigma2_theta,
                      double rho,
                      arma::vec eta,
                      double sigma2_eta,
                      Rcpp::List temporal_corr_info,
                      double a_phi,
                      double b_phi,
                      double metrop_var_phi_trans,
                      int acctot_phi_trans){
  
int s = theta.n_rows;
int m = eta.size();
arma::vec theta_full(s*m); theta_full.fill(0.00);
arma::vec eta_full(s*m); eta_full.fill(0.00);
for(int j = 0; j < s; ++j){
   theta_full.subvec(m*j, ((j + 1)*m - 1)) = trans(theta.row(j));
   eta_full.subvec(m*j, ((j + 1)*m - 1)) = eta;
   }

/*Second*/
Rcpp::List temporal_corr_info_old = temporal_corr_info;
arma::mat corr_inv_old = temporal_corr_info_old[0];
double log_deter_old = temporal_corr_info_old[1];
double phi_trans_old = log((phi_old - a_phi)/(b_phi - phi_old));

double second = -0.50*(1 + s)*log_deter_old - 
                (1.00/sigma2_eta)*0.50*dot(eta, (corr_inv_old*eta)) -
                (1.00/sigma2_theta)*0.50*dot((theta_full - eta_full), kron((rho*MCAR + (1.00 - rho)*eye(s, s)), corr_inv_old)*(theta_full - eta_full)) + 
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
               (1.00/sigma2_eta)*0.50*dot(eta, (corr_inv*eta)) -
               (1.00/sigma2_theta)*0.50*dot((theta_full - eta_full), kron((rho*MCAR + (1.00 - rho)*eye(s, s)), corr_inv)*(theta_full - eta_full)) +
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




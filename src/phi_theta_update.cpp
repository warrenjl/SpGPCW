#include "RcppArmadillo.h"
#include "SpGPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List phi_theta_update(double phi_theta_old,
                            arma::vec theta,
                            double sigma2_theta,
                            Rcpp::List temporal_corr_info_theta,
                            double a_phi_theta,
                            double b_phi_theta,
                            double metrop_var_phi_theta_trans,
                            int acctot_phi_theta_trans){
  
int m = theta.size();

/*Second*/
Rcpp::List temporal_corr_info_theta_old = temporal_corr_info_theta;
arma::mat corr_inv_theta_old = temporal_corr_info_theta_old[0];
double log_deter_theta_old = temporal_corr_info_theta_old[1];
double phi_theta_trans_old = log((phi_theta_old - a_phi_theta)/(b_phi_theta - phi_theta_old));

double second = -0.50*log_deter_theta_old - 
                (1.00/sigma2_theta)*0.50*dot(theta, (corr_inv_theta_old*theta)) -
                phi_theta_trans_old -
                2.00*log(1.00 + exp(phi_theta_trans_old));

/*First*/
double phi_theta_trans = R::rnorm(phi_theta_trans_old, 
                                  sqrt(metrop_var_phi_theta_trans));
double phi_theta = (b_phi_theta*exp(phi_theta_trans) + a_phi_theta)/(exp(phi_theta_trans) + 1.00);
temporal_corr_info_theta = temporal_corr_fun(m, phi_theta);
arma::mat corr_inv_theta = temporal_corr_info_theta[0];
double log_deter_theta = temporal_corr_info_theta[1];

double first = -0.50*log_deter_theta - 
               (1.00/sigma2_theta)*0.50*dot(theta, (corr_inv_theta*theta)) -
               phi_theta_trans -
               2.00*log(1.00 + exp(phi_theta_trans));

/*Decision*/
double ratio = exp(first - second);   
int acc = 1;
if(ratio < R::runif(0.00, 1.00)){
  phi_theta = phi_theta_old;
  temporal_corr_info_theta = temporal_corr_info_theta_old;
  acc = 0;
  }
acctot_phi_theta_trans = acctot_phi_theta_trans + 
                         acc;

return Rcpp::List::create(Rcpp::Named("phi_theta") = phi_theta,
                          Rcpp::Named("acctot_phi_theta_trans") = acctot_phi_theta_trans,
                          Rcpp::Named("temporal_corr_info_theta") = temporal_corr_info_theta);

}




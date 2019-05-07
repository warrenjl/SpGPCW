#include "RcppArmadillo.h"
#include "SpGPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List SpGPCW(int mcmc_samples,
                  arma::vec y,
                  arma::mat x,
                  arma::mat z,
                  arma::vec site_id,
                  arma::mat neighbors,
                  double metrop_var_rho_trans,
                  double metrop_var_phi_trans,
                  Rcpp::Nullable<double> sigma2_beta_prior = R_NilValue,
                  Rcpp::Nullable<double> alpha_sigma2_theta_prior = R_NilValue,
                  Rcpp::Nullable<double> beta_sigma2_theta_prior = R_NilValue,
                  Rcpp::Nullable<double> a_rho_prior = R_NilValue,
                  Rcpp::Nullable<double> b_rho_prior = R_NilValue,
                  Rcpp::Nullable<double> alpha_sigma2_eta_prior = R_NilValue,
                  Rcpp::Nullable<double> beta_sigma2_eta_prior = R_NilValue,
                  Rcpp::Nullable<double> a_phi_prior = R_NilValue,
                  Rcpp::Nullable<double> b_phi_prior = R_NilValue,
                  Rcpp::Nullable<Rcpp::NumericVector> beta_init = R_NilValue,
                  Rcpp::Nullable<Rcpp::NumericMatrix> theta_init = R_NilValue,
                  Rcpp::Nullable<double> sigma2_theta_init = R_NilValue,
                  Rcpp::Nullable<double> rho_init = R_NilValue,
                  Rcpp::Nullable<Rcpp::NumericVector> eta_init = R_NilValue,
                  Rcpp::Nullable<double> sigma2_eta_init = R_NilValue,
                  Rcpp::Nullable<double> phi_init = R_NilValue,
                  Rcpp::Nullable<int> rho_zero_indicator = R_NilValue){

//Defining Parameters and Quantities of Interest
int p_x = x.n_cols;
int m = z.n_cols;
int s = neighbors.n_cols;
arma::mat beta(p_x, mcmc_samples); beta.fill(0.00);
Rcpp::List theta(mcmc_samples);
for(int j = 0; j < mcmc_samples; ++j){
   arma::mat theta_temp(s, m); theta_temp.fill(0.00);
   theta[j] = theta_temp;
   }
arma::vec sigma2_theta(mcmc_samples); sigma2_theta.fill(0.00);
arma::vec rho(mcmc_samples); rho.fill(0.00);
arma::mat eta(m, mcmc_samples); eta.fill(0.00);
arma::vec sigma2_eta(mcmc_samples); sigma2_eta.fill(0.00);
arma::vec phi(mcmc_samples); phi.fill(0.00);
arma::vec neg_two_loglike(mcmc_samples); neg_two_loglike.fill(0.00);

//Miscellaneous Information
arma::vec diag_neighbors(s); diag_neighbors.fill(0.00);
arma::mat z_star((s*m), m);
for(int j = 0; j < s; ++j){
   diag_neighbors(j) = sum(neighbors.row(j));
   z_star.submat(m*j, 0, ((j + 1)*m - 1), (m-1)) = eye(m, m);
   }
arma::mat MCAR = diagmat(diag_neighbors) - 
                 neighbors;

//Prior Information
double sigma2_beta = 10000.00;
if(sigma2_beta_prior.isNotNull()){
  sigma2_beta = Rcpp::as<double>(sigma2_beta_prior);
  }

double alpha_sigma2_theta = 3.00;
if(alpha_sigma2_theta_prior.isNotNull()){
  alpha_sigma2_theta = Rcpp::as<double>(alpha_sigma2_theta_prior);
  }
  
double beta_sigma2_theta = 2.00;
if(beta_sigma2_theta_prior.isNotNull()){
  beta_sigma2_theta = Rcpp::as<double>(beta_sigma2_theta_prior);
  }

double a_rho = 0.00;
if(a_rho_prior.isNotNull()){
  a_rho = Rcpp::as<double>(a_rho_prior);
  }

double b_rho = 1.00;
if(b_rho_prior.isNotNull()){
  b_rho = Rcpp::as<double>(b_rho_prior);
  }

double alpha_sigma2_eta = 3.00;
if(alpha_sigma2_eta_prior.isNotNull()){
  alpha_sigma2_eta = Rcpp::as<double>(alpha_sigma2_eta_prior);
  }

double beta_sigma2_eta = 2.00;
if(beta_sigma2_eta_prior.isNotNull()){
  beta_sigma2_eta = Rcpp::as<double>(beta_sigma2_eta_prior);
  }

double a_phi = log(0.9999)/(-(m - 1.00));  
if(a_phi_prior.isNotNull()){
  a_phi = Rcpp::as<double>(a_phi_prior);
  }

double b_phi = log(0.0001)/(-1.00);
if(b_phi_prior.isNotNull()){
  b_phi = Rcpp::as<double>(b_phi_prior);
  }

//Initial Values
beta.col(0).fill(0.00);
if(beta_init.isNotNull()){
  beta.col(0) = Rcpp::as<arma::vec>(beta_init);
  }

arma::mat theta_temp(s, m); theta_temp.fill(0.00);
if(theta_init.isNotNull()){
  theta_temp = Rcpp::as<arma::mat>(theta_init);
  }
theta[0] = theta_temp;

sigma2_theta(0) = 1.00;
if(sigma2_theta_init.isNotNull()){
  sigma2_theta(0) = Rcpp::as<double>(sigma2_theta_init);
  }

rho(0) = (b_rho - a_rho)*0.50;
if(rho_init.isNotNull()){
  rho(0) = Rcpp::as<double>(rho_init);
  }

eta.col(0).fill(0.00);
if(eta_init.isNotNull()){
  eta.col(0) = Rcpp::as<arma::vec>(eta_init);
  }

sigma2_eta(0) = 1.00;
if(sigma2_eta_init.isNotNull()){
  sigma2_eta(0) = Rcpp::as<double>(sigma2_eta_init);
  }

phi(0) = (b_phi - a_phi)*0.01;
if(phi_init.isNotNull()){
  phi(0) = Rcpp::as<double>(phi_init);
  }

Rcpp::List temporal_corr_info = temporal_corr_fun(m, phi(0));
neg_two_loglike(0) = neg_two_loglike_update(y,
                                            x,
                                            z,
                                            site_id,
                                            beta.col(0),
                                            theta[0]);

//Non Spatial Option (\rho fixed at 0):
//rho_zero = 0; Spatial
//rho_zero = Any Other Integer (Preferably One); Non Spatial
int rho_zero = 0;
if(rho_zero_indicator.isNotNull()){
  rho_zero = Rcpp::as<int>(rho_zero_indicator);
  }

//Metropolis Settings
int acctot_rho_trans = 0;
int acctot_phi_trans = 0;

//Main Sampling Loop
for(int j = 1; j < mcmc_samples; ++j){

   //w Update
   Rcpp::List w_output = w_update(y,
                                  x,
                                  z,
                                  site_id,
                                  beta.col(j-1),
                                  theta[j-1]);
   arma::vec w = w_output[0];
   arma::vec gamma = w_output[1];
  
   //beta Update
   beta.col(j) = beta_update(x, 
                             z,
                             site_id,
                             w,
                             gamma,
                             theta[j-1],
                             sigma2_beta);
   
   //theta Update
   theta[j] = theta_update(theta[j-1],
                           x, 
                           z,
                           site_id,
                           neighbors,
                           w,
                           gamma,
                           beta.col(j),
                           sigma2_theta(j-1),
                           rho(j-1),
                           eta.col(j-1),
                           temporal_corr_info[0]);
   
   //sigma2_theta Update
   sigma2_theta(j) = sigma2_theta_update(MCAR,
                                         theta[j],
                                         rho(j-1),
                                         eta.col(j-1),
                                         temporal_corr_info[0],
                                         alpha_sigma2_theta,
                                         beta_sigma2_theta);
   
   //rho Update
   //Only if rho_zero = 0
   rho(j) = 0;
   if(rho_zero == 0){
     Rcpp::List rho_output = rho_update(rho(j-1),
                                        MCAR,
                                        theta[j],
                                        sigma2_theta(j),
                                        eta.col(j-1),
                                        temporal_corr_info[0],
                                        a_rho,
                                        b_rho,
                                        metrop_var_rho_trans,
                                        acctot_rho_trans);
   
     rho(j) = Rcpp::as<double>(rho_output[0]);
     acctot_rho_trans = rho_output[1];
     }
   
   //eta Update
   eta.col(j) = eta_update(MCAR,
                           z_star,
                           theta[j],
                           sigma2_theta(j),
                           rho(j),
                           sigma2_eta(j-1),
                           phi(j-1),
                           temporal_corr_info[0]);
   
   //sigma2_eta Update
   sigma2_eta(j) = sigma2_eta_update(eta.col(j),
                                     temporal_corr_info[0],
                                     alpha_sigma2_eta,
                                     beta_sigma2_eta);
   
   //phi Update
   Rcpp::List phi_output = phi_update(phi(j-1),
                                      MCAR,
                                      theta[j],
                                      sigma2_theta(j),
                                      rho(j),
                                      eta.col(j),
                                      sigma2_eta(j),
                                      temporal_corr_info,
                                      a_phi,
                                      b_phi,
                                      metrop_var_phi_trans,
                                      acctot_phi_trans);
   
   phi(j) = Rcpp::as<double>(phi_output[0]);
   acctot_phi_trans = phi_output[1];
   temporal_corr_info = phi_output[2];

   //neg_two_loglike Update
   neg_two_loglike(j) = neg_two_loglike_update(y,
                                               x,
                                               z,
                                               site_id,
                                               beta.col(j),
                                               theta[j]);
   
   //Progress
   if((j + 1) % 10 == 0){ 
     Rcpp::checkUserInterrupt();
     }
  
   if(((j + 1) % int(round(mcmc_samples*0.05)) == 0)){
     double completion = round(100*((j + 1)/(double)mcmc_samples));
     Rcpp::Rcout << "Progress: " << completion << "%" << std::endl;
     if(rho_zero == 0){
       double accrate_rho_trans = round(100*(acctot_rho_trans/(double)j));
       Rcpp::Rcout << "rho Acceptance: " << accrate_rho_trans << "%" << std::endl;
       }
     double accrate_phi_trans = round(100*(acctot_phi_trans/(double)j));
     Rcpp::Rcout << "phi Acceptance: " << accrate_phi_trans << "%" << std::endl;
     Rcpp::Rcout << "*******************" << std::endl;
     }
  
   }
                                  
return Rcpp::List::create(Rcpp::Named("beta") = beta,
                          Rcpp::Named("theta") = theta,
                          Rcpp::Named("sigma2_theta") = sigma2_theta,
                          Rcpp::Named("rho") = rho,
                          Rcpp::Named("eta") = eta,
                          Rcpp::Named("sigma2_eta") = sigma2_eta,
                          Rcpp::Named("phi") = phi,
                          Rcpp::Named("neg_two_loglike") = neg_two_loglike,
                          Rcpp::Named("acctot_rho_trans") = acctot_rho_trans,
                          Rcpp::Named("acctot_phi_trans") = acctot_phi_trans);

}


#include "RcppArmadillo.h"
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
                  double metrop_var_phi_theta_trans,
                  double metrop_var_rho_trans,
                  double metrop_var_phi_eta_trans,
                  Rcpp::Nullable<double> sigma2_beta_prior = R_NilValue,
                  Rcpp::Nullable<double> alpha_sigma2_theta_prior = R_NilValue,
                  Rcpp::Nullable<double> beta_sigma2_theta_prior = R_NilValue,
                  Rcpp::Nullable<double> a_phi_theta_prior = R_NilValue,
                  Rcpp::Nullable<double> b_phi_theta_prior = R_NilValue,
                  Rcpp::Nullable<double> a_rho_prior = R_NilValue,
                  Rcpp::Nullable<double> b_rho_prior = R_NilValue,
                  Rcpp::Nullable<double> alpha_sigma2_eta_prior = R_NilValue,
                  Rcpp::Nullable<double> beta_sigma2_eta_prior = R_NilValue,
                  Rcpp::Nullable<double> a_phi_eta_prior = R_NilValue,
                  Rcpp::Nullable<double> b_phi_eta_prior = R_NilValue,
                  Rcpp::Nullable<Rcpp::NumericVector> beta_init = R_NilValue,
                  Rcpp::Nullable<Rcpp::NumericVector> theta_init = R_NilValue,
                  Rcpp::Nullable<double> sigma2_theta_init = R_NilValue,
                  Rcpp::Nullable<double> phi_theta_init = R_NilValue,
                  Rcpp::Nullable<Rcpp::NumericMatrix> eta_init = R_NilValue,
                  Rcpp::Nullable<double> rho_init = R_NilValue,
                  Rcpp::Nullable<double> sigma2_eta_init = R_NilValue,
                  Rcpp::Nullable<double> phi_eta_init = R_NilValue,
                  Rcpp::Nullable<int> rho_zero_indicator = R_NilValue){

//Defining Parameters and Quantities of Interest
arma::mat beta(x.n_cols, mcmc_samples); beta.fill(0.00);
arma::mat theta(z.n_cols, mcmc_samples); theta.fill(0.00);
arma::vec sigma2_theta(mcmc_samples); sigma2_theta.fill(0.00);
arma::vec phi_theta(mcmc_samples); phi_theta.fill(0.00);
Rcpp::List eta(mcmc_samples);
for(int j = 0; j < mcmc_samples; ++j){
   arma::mat eta_temp(neighbors.n_cols, z.n_cols); eta_temp.fill(0.00);
   eta[j] = eta_temp;
   }
arma::vec rho(mcmc_samples); rho.fill(0.00);
arma::vec sigma2_eta(mcmc_samples); sigma2_eta.fill(0.00);
arma::vec phi_eta(mcmc_samples); phi_eta.fill(0.00);
arma::vec neg_two_loglike(mcmc_samples); neg_two_loglike.fill(0.00);

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

double a_phi_theta = log(0.9999)/(-(z.n_cols - 1));  
if(a_phi_theta_prior.isNotNull()){
  a_phi_theta = Rcpp::as<double>(a_phi_theta_prior);
  }

double b_phi_theta = log(0.0001)/(-1);
if(b_phi_theta_prior.isNotNull()){
  b_phi_theta = Rcpp::as<double>(b_phi_theta_prior);
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

double a_phi_eta = log(0.9999)/(-(z.n_cols - 1));  
if(a_phi_eta_prior.isNotNull()){
  a_phi_eta = Rcpp::as<double>(a_phi_eta_prior);
  }
  
double b_phi_eta = log(0.0001)/(-1);
if(b_phi_eta_prior.isNotNull()){
  b_phi_eta = Rcpp::as<double>(b_phi_eta_prior);
  }

//Initial Values
beta.col(0).fill(0.00);
if(beta_init.isNotNull()){
  beta.col(0) = Rcpp::as<arma::vec>(beta_init);
  }

theta.col(0).fill(0.00);
if(theta_init.isNotNull()){
  theta.col(0) = Rcpp::as<arma::vec>(theta_init);
  }

sigma2_theta(0) = 1.00;
if(sigma2_theta_init.isNotNull()){
  sigma2_theta(0) = Rcpp::as<double>(sigma2_theta_init);
  }

phi_theta(0) = (b_phi_theta - a_phi_theta)*0.01;
if(phi_theta_init.isNotNull()){
  phi_theta(0) = Rcpp::as<double>(phi_theta_init);
  }

arma::mat eta_temp(neighbors.n_cols, z.n_cols); eta_temp.fill(0.00);
if(eta_init.isNotNull()){
  eta_temp = Rcpp::as<arma::mat>(eta_init);
  }
eta[0] = eta_temp;

rho(0) = (b_rho - a_rho)*0.50;
if(rho_init.isNotNull()){
  rho(0) = Rcpp::as<double>(rho_init);
  }

sigma2_eta(0) = 1.00;
if(sigma2_eta_init.isNotNull()){
  sigma2_eta(0) = Rcpp::as<double>(sigma2_eta_init);
  }

phi_eta(0) = (b_phi_eta - a_phi_eta)*0.01;
if(phi_eta_init.isNotNull()){
  phi_eta(0) = Rcpp::as<double>(phi_eta_init);
  }

Rcpp::List temporal_corr_info_theta = temporal_corr_fun(z.n_cols, phi_theta(0));
Rcpp::List temporal_corr_info_eta = temporal_corr_fun(z.n_cols, phi_eta(0));
neg_two_loglike(0) = neg_two_loglike_update(y,
                                            x,
                                            z,
                                            site_id,
                                            beta.col(0),
                                            theta.col(0),
                                            eta[0]);

//Non Spatial Option (\rho fixed at 0):
//rho_zero = 0; Spatial
//rho_zero = Any Other Integer (Preferably One); Non Spatial
int rho_zero = 0;
if(rho_zero_indicator.isNotNull()){
  rho_zero = Rcpp::as<int>(rho_zero_indicator);
  }

//Metropolis Settings
int acctot_phi_theta_trans = 0;
int acctot_rho_trans = 0;
int acctot_phi_eta_trans = 0;

//Main Sampling Loop
for(int j = 1; j < mcmc_samples; ++j){
  
   //w Update
   Rcpp::List w_output = w_update(y,
                                  x,
                                  z,
                                  site_id,
                                  beta.col(j-1),
                                  theta.col(j-1),
                                  eta[j-1]);
   arma::vec w = w_output[0];
   arma::vec gamma = w_output[1];
  
   //beta Update
   beta.col(j) = beta_update(x, 
                             z,
                             site_id,
                             sigma2_beta,
                             w,
                             gamma,
                             theta.col(j-1),
                             eta[j-1]);
   
   //theta Update
   theta.col(j) = theta_update(x, 
                               z,
                               site_id,
                               w,
                               gamma,
                               beta.col(j),
                               eta[j-1],
                               sigma2_theta(j-1),
                               temporal_corr_info_theta(0));
   
   //sigma2_theta Update
   sigma2_theta(j) = sigma2_theta_update(theta.col(j),
                                         temporal_corr_info_theta(0),
                                         alpha_sigma2_theta,
                                         beta_sigma2_theta);
   
   //phi_theta Update
   Rcpp::List phi_theta_output = phi_theta_update(phi_theta(j-1),
                                                  theta.col(j),
                                                  sigma2_theta(j),
                                                  temporal_corr_info_theta,
                                                  a_phi_theta,
                                                  b_phi_theta,
                                                  metrop_var_phi_theta_trans,
                                                  acctot_phi_theta_trans);
   
   phi_theta(j) = phi_theta_output[0];
   acctot_phi_theta_trans = phi_theta_output[1];
   temporal_corr_info_theta = phi_theta_output[2];
   
   //eta Update
   eta[j] = eta_update(eta[j-1],
                       x, 
                       z,
                       site_id,
                       neighbors,
                       w,
                       gamma,
                       beta.col(j),
                       theta.col(j),
                       rho(j-1),
                       sigma2_eta(j-1),
                       temporal_corr_info_eta[0]);
  
   //rho Update
   //Only if rho_zero = 0
   rho(j) = 0;
   if(rho_zero == 0){
     Rcpp::List rho_output = rho_update(rho(j-1),
                                        neighbors,
                                        eta[j],
                                        sigma2_eta(j-1),
                                        temporal_corr_info_eta[0],
                                        a_rho,
                                        b_rho,
                                        metrop_var_rho_trans,
                                        acctot_rho_trans);
   
     rho(j) = rho_output[0];
     acctot_rho_trans = rho_output[1];
     }
   
   //sigma2_eta Update
   sigma2_eta(j) = sigma2_eta_update(neighbors,
                                     eta[j],
                                     rho(j),
                                     temporal_corr_info_eta(0),
                                     alpha_sigma2_eta,
                                     beta_sigma2_eta);
   
   //phi_eta Update
   Rcpp::List phi_eta_output = phi_eta_update(phi_eta(j-1),
                                              neighbors,
                                              eta[j],
                                              rho(j),
                                              sigma2_eta(j),
                                              temporal_corr_info_eta,
                                              a_phi_eta,
                                              b_phi_eta,
                                              metrop_var_phi_eta_trans,
                                              acctot_phi_eta_trans);
     
   phi_eta(j) = phi_eta_output[0];
   acctot_phi_eta_trans = phi_eta_output[1];
   temporal_corr_info_eta = phi_eta_output[2];

   //neg_two_loglike Update
   neg_two_loglike(j) = neg_two_loglike_update(y,
                                               x,
                                               z,
                                               site_id,
                                               beta.col(j),
                                               theta.col(j),
                                               eta[j]);
   
   //Progress
   if((j + 1) % 10 == 0){ 
     Rcpp::checkUserInterrupt();
     }
  
   if(((j + 1) % int(round(mcmc_samples*0.05)) == 0)){
     double completion = round(100*((j + 1)/(double)mcmc_samples));
     Rcpp::Rcout << "Progress: " << completion << "%" << std::endl;
     double accrate_phi_theta_trans = round(100*(acctot_phi_theta_trans/(double)j));
     Rcpp::Rcout << "phi_theta Acceptance: " << accrate_phi_theta_trans << "%" << std::endl;
     if(rho_zero == 0){
       double accrate_rho_trans = round(100*(acctot_rho_trans/(double)j));
       Rcpp::Rcout << "rho Acceptance: " << accrate_rho_trans << "%" << std::endl;
       }
     double accrate_phi_eta_trans = round(100*(acctot_phi_eta_trans/(double)j));
     Rcpp::Rcout << "phi_eta Acceptance: " << accrate_phi_eta_trans << "%" << std::endl;
     Rcpp::Rcout << "*************************" << std::endl;
     }
  
   }
                                  
return Rcpp::List::create(Rcpp::Named("beta") = beta,
                          Rcpp::Named("theta") = theta,
                          Rcpp::Named("sigma2_theta") = sigma2_theta,
                          Rcpp::Named("phi_theta") = phi_theta,
                          Rcpp::Named("eta") = eta,
                          Rcpp::Named("rho") = rho,
                          Rcpp::Named("sigma2_eta") = sigma2_eta,
                          Rcpp::Named("phi_eta") = phi_eta,
                          Rcpp::Named("neg_two_loglike") = neg_two_loglike,
                          Rcpp::Named("acctot_rho_trans") = acctot_rho_trans,
                          Rcpp::Named("acctot_phi_eta_trans") = acctot_phi_eta_trans);

}


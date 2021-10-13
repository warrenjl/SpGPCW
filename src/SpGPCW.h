#ifndef __SpGPCW__
#define __SpGPCW__

arma::vec rcpp_pgdraw(arma::vec b, 
                      arma::vec c);

Rcpp::List temporal_corr_fun(int m,
                             double phi);

Rcpp::List w_update(arma::vec y,
                    arma::mat x,
                    arma::mat z,
                    arma::vec site_id,
                    arma::vec beta_old,
                    arma::mat theta_old);

arma::vec beta_update(arma::mat x, 
                      arma::mat z,
                      arma::vec site_id,
                      arma::vec w,
                      arma::vec gamma,
                      arma::mat theta_old,
                      double sigma2_beta);

arma::mat theta_update(arma::mat theta_old,
                       arma::mat x, 
                       arma::mat z,
                       arma::vec site_id,
                       arma::mat neighbors,
                       arma::vec w,
                       arma::vec gamma,
                       arma::vec beta,
                       double sigma2_theta_old,
                       double rho_old,
                       arma::vec eta_old,
                       arma::mat corr_inv);

double sigma2_theta_update(arma::mat MCAR,
                           arma::mat theta,
                           double rho_old,
                           arma::vec eta_old,
                           arma::mat corr_inv,
                           double alpha_sigma2_theta,
                           double beta_sigma2_theta);

Rcpp::List rho_update(double rho_old,
                      arma::mat MCAR,
                      arma::mat theta,
                      double sigma2_theta,
                      arma::vec eta_old,
                      arma::mat corr_inv,
                      double a_rho,
                      double b_rho,
                      double metrop_var_rho_trans,
                      int acctot_rho_trans);

arma::vec eta_update(arma::mat MCAR,
                     arma::mat z_star,
                     arma::mat theta,
                     double sigma2_theta,
                     double rho,
                     double sigma2_eta_old,
                     double phi_old,
                     arma::mat corr_inv);

double sigma2_eta_update(arma::vec eta,
                         arma::mat corr_inv,
                         double alpha_sigma2_eta,
                         double beta_sigma2_eta);

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
                      int acctot_phi_trans);

double neg_two_loglike_update(arma::vec y,
                              arma::mat x,
                              arma::mat z,
                              arma::vec site_id,
                              arma::vec beta,
                              arma::mat theta);

Rcpp::List SpGPCW(int mcmc_samples,
                  arma::vec y,
                  arma::mat x,
                  arma::mat z,
                  arma::vec site_id,
                  arma::mat neighbors,
                  double metrop_var_rho_trans,
                  double metrop_var_phi_trans,
                  Rcpp::Nullable<double> sigma2_beta_prior,
                  Rcpp::Nullable<double> alpha_sigma2_theta_prior,
                  Rcpp::Nullable<double> beta_sigma2_theta_prior,
                  Rcpp::Nullable<double> a_rho_prior,
                  Rcpp::Nullable<double> b_rho_prior,
                  Rcpp::Nullable<double> alpha_sigma2_eta_prior,
                  Rcpp::Nullable<double> beta_sigma2_eta_prior,
                  Rcpp::Nullable<double> a_phi_prior,
                  Rcpp::Nullable<double> b_phi_prior,
                  Rcpp::Nullable<Rcpp::NumericVector> beta_init,
                  Rcpp::Nullable<Rcpp::NumericMatrix> theta_init,
                  Rcpp::Nullable<double> sigma2_theta_init,
                  Rcpp::Nullable<double> rho_init,
                  Rcpp::Nullable<Rcpp::NumericVector> eta_init,
                  Rcpp::Nullable<double> sigma2_eta_init,
                  Rcpp::Nullable<double> phi_init,
                  Rcpp::Nullable<int> rho_zero_indicator); 

#endif // __SpGPCW__

#ifndef __SpGPCW__
#define __SpGPCW__

arma::vec rcpp_pgdraw(double b, 
                      arma::vec c);

Rcpp::List temporal_corr_fun(int m,
                             double phi);

Rcpp::List w_update(arma::vec y,
                    arma::mat x,
                    arma::mat z,
                    arma::vec site_id,
                    arma::vec beta_old,
                    arma::vec theta_old,
                    arma::mat eta_old);

arma::vec beta_update(arma::mat x, 
                      arma::mat z,
                      arma::vec site_id,
                      double sigma2_beta,
                      arma::vec w,
                      arma::vec gamma,
                      arma::vec theta_old,
                      arma::mat eta_old);

arma::vec theta_update(arma::mat x, 
                       arma::mat z,
                       arma::vec site_id,
                       arma::vec w,
                       arma::vec gamma,
                       arma::vec beta,
                       arma::mat eta_old,
                       double sigma2_theta_old,
                       arma::mat corr_inv_theta);

double sigma2_theta_update(arma::vec theta,
                           arma::mat corr_inv_theta,
                           double alpha_sigma2_theta,
                           double beta_sigma2_theta);

Rcpp::List phi_theta_update(double phi_theta_old,
                            arma::vec theta,
                            double sigma2_theta,
                            Rcpp::List temporal_corr_info_theta,
                            double a_phi_theta,
                            double b_phi_theta,
                            double metrop_var_phi_theta_trans,
                            int acctot_phi_theta_trans);

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
                     arma::mat corr_inv_eta);

Rcpp::List rho_update(double rho_old,
                      arma::mat neighbors,
                      arma::mat eta,
                      double sigma2_eta_old,
                      arma::mat corr_inv_eta,
                      double a_rho,
                      double b_rho,
                      double metrop_var_rho_trans,
                      int acctot_rho_trans);

double sigma2_eta_update(arma::mat neighbors,
                         arma::mat eta,
                         double rho,
                         arma::mat corr_inv_eta,
                         double alpha_sigma2_eta,
                         double beta_sigma2_eta);

Rcpp::List phi_eta_update(double phi_eta_old,
                          arma::mat neighbors,
                          arma::mat eta,
                          double rho,
                          double sigma2_eta,
                          Rcpp::List temporal_corr_info_eta,
                          double a_phi_eta,
                          double b_phi_eta,
                          double metrop_var_phi_eta_trans,
                          int acctot_phi_eta_trans);

double neg_two_loglike_update(arma::vec y,
                              arma::mat x,
                              arma::mat z,
                              arma::vec site_id,
                              arma::vec beta,
                              arma::vec theta,
                              arma::mat eta);

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
                  Rcpp::Nullable<double> a_phi_theta_prior,
                  Rcpp::Nullable<double> b_phi_theta_prior,
                  Rcpp::Nullable<double> a_rho_prior,
                  Rcpp::Nullable<double> b_rho_prior,
                  Rcpp::Nullable<double> alpha_sigma2_eta_prior,
                  Rcpp::Nullable<double> beta_sigma2_eta_prior,
                  Rcpp::Nullable<double> a_phi_eta_prior,
                  Rcpp::Nullable<double> b_phi_eta_prior,
                  Rcpp::Nullable<Rcpp::NumericVector> beta_init,
                  Rcpp::Nullable<Rcpp::NumericVector> theta_init,
                  Rcpp::Nullable<double> sigma2_theta_init,
                  Rcpp::Nullable<double> phi_theta_init,
                  Rcpp::Nullable<Rcpp::NumericMatrix> eta_init,
                  Rcpp::Nullable<double> rho_init,
                  Rcpp::Nullable<double> sigma2_eta_init,
                  Rcpp::Nullable<double> phi_eta_init,
                  Rcpp::Nullable<int> rho_zero_indicator); 

#endif // __SpGPCW__

#include "RcppArmadillo.h"
#include "SpGPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]



Rcpp::List rho_update(double rho_old,
                      arma::mat neighbors,
                      arma::mat eta,
                      double sigma2_eta_old,
                      arma::mat corr_inv_eta,
                      double a_rho,
                      double b_rho,
                      double metrop_var_rho_trans,
                      int acctot_rho_trans){
  
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
double rho_trans_old = log((rho_old - a_rho)/(b_rho - rho_old));
double MCAR_info_old = 0.00; 
double sign_old = 0.00;     
log_det(MCAR_info_old, sign_old, (rho_old*MCAR + (1.00 - rho_old)*eye(s, s)));

double second = 0.50*m*MCAR_info_old - 
                (1.00/sigma2_eta_old)*0.50*dot(eta_full, kron((rho_old*MCAR + (1.00 - rho_old)*eye(s, s)), corr_inv_eta)*eta_full) + 
                rho_trans_old -
                2.00*log(1 + exp(rho_trans_old));

/*First*/
double rho_trans = R::rnorm(rho_trans_old, 
                            sqrt(metrop_var_rho_trans));
double rho = (b_rho*exp(rho_trans) + a_rho)/(exp(rho_trans) + 1.00);
double MCAR_info = 0.00; 
double sign = 0.00;     
log_det(MCAR_info, sign, (rho*MCAR + (1.00 - rho)*eye(s, s)));

double first = 0.50*m*MCAR_info - 
               (1.00/sigma2_eta_old)*0.50*dot(eta_full, kron((rho*MCAR + (1.00 - rho)*eye(s, s)), corr_inv_eta)*eta_full) + 
               rho_trans -
               2.00*log(1.00 + exp(rho_trans));

/*Decision*/
double ratio = exp(first - second);   
double acc = 1;
if(ratio < R::runif(0.00, 1.00)){
  rho = rho_old;
  acc = 0;
  }
acctot_rho_trans = acctot_rho_trans + 
                   acc;

return Rcpp::List::create(Rcpp::Named("rho") = rho,
                          Rcpp::Named("acctot_rho_trans") = acctot_rho_trans);

}




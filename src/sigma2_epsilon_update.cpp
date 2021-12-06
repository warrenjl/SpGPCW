#include "RcppArmadillo.h"
#include "SpGPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double sigma2_epsilon_update(arma::vec y,
                             arma::mat x,
                             arma::mat z,
                             arma::vec site_id,
                             arma::vec off_set,
                             double a_sigma2_epsilon,
                             double b_sigma2_epsilon,
                             arma::vec beta_old,
                             arma::mat theta_old){

int n = y.size();
int s = theta_old.n_rows;
arma::vec mu(n); mu.fill(0.00);
for(int j = 0; j < s; ++j){
    
   arma::uvec ids = find(site_id == (j + 1));
   mu.elem(ids) = off_set.elem(ids) +
                  x.rows(ids)*beta_old +
                  z.rows(ids)*trans(theta_old.row(j));
    
   }

double a_sigma2_epsilon_update = 0.50*n + 
                                 a_sigma2_epsilon;

double b_sigma2_epsilon_update = 0.50*dot((y - mu), (y - mu)) + 
                                 b_sigma2_epsilon;

double sigma2_epsilon = 1.00/R::rgamma(a_sigma2_epsilon_update,
                                       (1.00/b_sigma2_epsilon_update));

return(sigma2_epsilon);

}






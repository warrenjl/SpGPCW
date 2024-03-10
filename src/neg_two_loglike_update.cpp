#include "RcppArmadillo.h"
#include "SpGPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

double neg_two_loglike_update(arma::vec y,
                              arma::mat x,
                              arma::mat z,
                              arma::vec site_id,
                              arma::vec off_set,
                              arma::vec tri_als,
                              int likelihood_indicator,
                              int r,
                              double sigma2_epsilon,
                              arma::vec beta,
                              arma::mat theta){

int n = y.size();
int s = theta.n_rows;
arma::vec mu(n); mu.fill(0.00);
arma::vec dens(n); dens.fill(0.00);
for(int j = 0; j < s; ++j){
   
   arma::uvec ids = find(site_id == (j + 1));
   mu.elem(ids) = off_set.elem(ids) +
                  x.rows(ids)*beta +
                  z.rows(ids)*trans(theta.row(j));
   
   }

if(likelihood_indicator == 0){
   
  arma::vec probs = exp(mu)/(1.00 + exp(mu));
  for(int j = 0; j < n; ++j){
     dens(j) = R::dbinom(y(j),
                         tri_als(j),
                         probs(j),
                         TRUE);
     }
   
  }

if(likelihood_indicator == 1){
  for(int j = 0; j < n; ++j){
     dens(j) = R::dnorm(y(j),
                        mu(j),
                        sqrt(sigma2_epsilon),
                        TRUE);
     }
  }

if(likelihood_indicator == 2){
   
  arma::vec probs = exp(mu)/(1.00 + exp(mu));
  for(int j = 0; j < n; ++j){
     dens(j) = R::dnbinom(y(j), 
                          r, 
                          (1.00 - probs(j)),        
                          TRUE);
     }
   
  }

double neg_two_loglike = -2.00*sum(dens);

return neg_two_loglike;

}






















































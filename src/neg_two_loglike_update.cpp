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
                              arma::vec beta,
                              arma::mat theta){

int n = y.size();
int s = theta.n_rows;
arma::vec logit_probs(n); logit_probs.fill(0.00);
for(int j = 0; j < s; ++j){
   arma::uvec ids = find(site_id == (j + 1));
   logit_probs.elem(ids) = x.rows(ids)*beta +
                           z.rows(ids)*trans(theta.row(j));
   }
arma::vec probs = exp(logit_probs)/(1 + exp(logit_probs));

arma::vec dens(n); dens.fill(0.00);
for(int j = 0; j < n; ++j){
   dens(j) = R::dbinom(y(j),
                       1,
                       probs(j),
                       TRUE);
   }

double neg_two_loglike = -2.00*sum(dens);

return neg_two_loglike;

}






















































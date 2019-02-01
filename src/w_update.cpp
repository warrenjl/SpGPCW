#include "RcppArmadillo.h"
#include "SpGPCW.h"
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

Rcpp::List w_update(arma::vec y,
                    arma::mat x,
                    arma::mat z,
                    arma::vec site_id,
                    arma::vec beta_old,
                    arma::mat theta_old){

int n = y.size();
int s = theta_old.n_rows;
arma::vec mean_w(n); mean_w.fill(0.00);
for(int j = 0; j < s; ++j){
   arma::uvec ids = find(site_id == (j + 1));
   mean_w.elem(ids) = x.rows(ids)*beta_old + 
                      z.rows(ids)*trans(theta_old.row(j));
   }

arma::vec w = rcpp_pgdraw(1.00,
                          mean_w);

arma::vec gamma = (y - 0.50)/w;

return Rcpp::List::create(Rcpp::Named("w") = w,
                          Rcpp::Named("gamma") = gamma);

}

































































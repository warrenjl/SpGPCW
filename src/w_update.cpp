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
                    arma::vec off_set,
                    int likelihood_indicator,
                    int r,
                    arma::vec beta_old,
                    arma::mat theta_old){

int n = y.size();
int s = theta_old.n_rows;
arma::vec mean_w(n); mean_w.fill(0.00);
for(int j = 0; j < s; ++j){
  
   arma::uvec ids = find(site_id == (j + 1));
   mean_w.elem(ids) = off_set.elem(ids) +
                      x.rows(ids)*beta_old + 
                      z.rows(ids)*trans(theta_old.row(j));
   
   }

arma::vec input0(1); input0.fill(1.00);
arma::vec input2 = (r + y);

arma::vec w(n); w.fill(0.00);
arma::vec gamma(n); gamma.fill(0.00);

if(likelihood_indicator == 0){
  
  w = rcpp_pgdraw(input0,
                  mean_w);
  gamma = (y - 0.50)/w;
  
  } 

if(likelihood_indicator == 2){
  
  w = rcpp_pgdraw(input2,
                  mean_w);
  gamma = 0.50*(y - r)/w;
  
  }

return Rcpp::List::create(Rcpp::Named("w") = w,
                          Rcpp::Named("gamma") = gamma);

}

































































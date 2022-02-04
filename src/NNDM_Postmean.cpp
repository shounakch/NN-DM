#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]

arma::mat NNDM_pmean(int n, float p, int s, arma::mat mu_vec, Rcpp::List bw_inv, float df_n, float nun, arma::mat inputpt){
  
  arma::vec f_hat(s, fill::zeros);
  arma::mat middle_matrix(s, n, fill::zeros);
  
  for(int j=0; j<n; j++){
    
    arma::mat K_j = bw_inv(j);
    arma::mat R = chol(K_j);
    float detR = det(R);
    float detKj = pow(detR, 2.0);
    
    for(int i=0; i<s; i++){
      
      arma::vec locvec = (mu_vec.row(j)).t();
      arma::vec dummyvec = (inputpt.row(i)).t();
      arma::vec diffvec = dummyvec - locvec;
      arma::mat expr(1, 1, fill::zeros);
      
      expr = (diffvec.t())*(K_j*diffvec);
      
      float dens_ij = (pow(detKj, 0.5)) * (pow((1 + (expr(0,0) / df_n)), -0.5*(df_n + p)));
      
      middle_matrix(i,j) = dens_ij;
      
    }
    
  }
  
  return middle_matrix;
  
}






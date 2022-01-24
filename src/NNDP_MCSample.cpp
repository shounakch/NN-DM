#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
double random_gamma_univ(float a) {
  return R::rgamma(a, 1);
}

// [[Rccp::export]]
double random_normal(float s){
  return R::rnorm(0.0, s);
}

// [[Rcpp::export]]

arma::mat fdraws_univ(int n, int s, int MC, arma::vec wtvec, arma::vec mu_vec, arma::vec rate_par, float gamman, float nun, arma::vec inputpt){
  
  arma::mat f_draws(s, MC, fill::zeros);
  
  for(int m=0; m<MC; m++){
    
    arma::vec gam(n, fill::zeros);
    arma::mat evalmat(s, n, fill::zeros);
    double pi = 3.14159265;
    
    for(int i=0; i<n; i++){
      float rg = wtvec(i);
      gam(i) = random_gamma_univ(rg);
    }
    
    arma::vec pi_vec = gam/accu(gam);
    
    for(int i=0; i<n; i++){
      
      float K = random_gamma_univ(gamman/2)/rate_par(i);
      float sigma = 1/K;
      float meanrow = mu_vec(i);
      float stdnorm  = random_normal(1.0);
      float theta = (sqrt(sigma/nun)*stdnorm) + meanrow;
      
      for(int j=0; j<s; j++){
        
        float dummyvec = inputpt(j);
        float diffvec = dummyvec - theta;
        
        float expr = pow(diffvec, 2.0)*K;
        evalmat(j,i) = (1/sqrt(sigma*(2.0)*pi))*exp(-0.5*expr);
      }
      
    }
    
    arma::vec draws = evalmat*pi_vec;
    
    f_draws.col(m) = draws;
    
    if((m+1)%500==0){
      Rcpp::Rcout << "Monte Carlo Sample: " << m+1 << std::endl;
    }
    
  }
  
  return f_draws;
  
}





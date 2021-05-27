#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
double random_gamma(float a) {
  return R::rgamma(a, 1.0);
}

// [[Rccp::export]]
double random_normal(float s){
  return R::rnorm(0.0, s);
}

// [[Rcpp::export]]

mat fdraws(int n, int s, int MC, vec wtvec, vec mu_vec, vec rate_par, float gamman, float nun, vec inputpt){
  
  mat f_draws(s, MC, fill::zeros);
  
  for(int m=0; m<MC; m++){
    
    vec gam(n, fill::zeros);
    mat evalmat(s, n, fill::zeros);
    double pi = 3.14159265;
    
    for(int i=0; i<n; i++){
      float rg = wtvec(i);
      gam(i) = random_gamma(rg);
    }
    
    vec pi_vec = gam/accu(gam);
    
    for(int i=0; i<n; i++){
      
      float K = random_gamma(gamman/2)/rate_par(i);
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
    
    vec draws = evalmat*pi_vec;
    
    f_draws.col(m) = draws;
    
    if((m+1)%100==0){
      Rcpp::Rcout << "Monte Carlo Sample: " << m+1 << std::endl;
    }
    
  }
  
  return f_draws;
  
}





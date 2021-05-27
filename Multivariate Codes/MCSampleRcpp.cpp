#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace arma;

// [[Rcpp::export]]
double random_gamma(float a) {
  return R::rgamma(a, 1);
}

// [[Rcpp::export]]
vec random_gaussian(int p){
  vec rng(p, fill::zeros);
  for(int i=0; i<p; ++i){
    rng(i) = R::rnorm(0, 1);
  }
  return rng;
}

// [[Rcpp::export]]
void set_seed(double seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}

// [[Rcpp::export]]

mat fdraws(int n, float p, int s, int MC, vec wtvec, mat mu_vec, Rcpp::List psi_inv, float gamman, float nun, mat inputpt){
  
  mat f_draws(s, MC, fill::zeros);
  double pi = 3.14159265;
  double con = -0.5*p;
  
  for(int m=0; m<MC; m++){
    
    vec gam(n, fill::zeros);
    mat evalmat(s, n, fill::zeros);
    
    for(int i=0; i<n; ++i){
      float rg = wtvec(i);
      gam(i) = random_gamma(rg);
    }
    
    vec pi_vec = gam/accu(gam);
    
    for(int i=0; i<n; ++i){
      
      mat psimat = psi_inv(i);
      mat K = wishrnd(psimat, gamman);
      mat R = chol(K);
      mat R1 = inv(R);
      mat sig = R1*(R1.t());
      float detR = det(R);
      mat sig1 = sig/nun;
      
      vec meanrow = ((mu_vec).row(i)).t();
      mat theta  = mvnrnd(meanrow, sig1, 1);
      vec locvec = theta.col(0);
       
      for(int j=0; j<s; ++j){
        
        vec dummyvec = (inputpt.row(j)).t();
        vec diffvec = dummyvec - locvec;
        mat expr(1, 1, fill::zeros);
        
        expr = (diffvec.t())*(K*diffvec);
        
        evalmat(j,i) = detR*exp(-0.5*expr(0,0))*pow((2*pi), con);
      
      }
      
    }
    
    vec draws = evalmat*pi_vec;
    
    f_draws.col(m) = draws;
    
    if((m+1)%500==0){
      Rcpp::Rcout << "Monte Carlo Sample: " << m+1 << std::endl;
    }
    
  }
  
  return f_draws;
  
}





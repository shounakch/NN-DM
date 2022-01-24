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
arma::vec random_gaussian(int p){
  arma::vec rng(p, fill::zeros);
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

arma::mat fdraws(int n, float p, int s, int MC, arma::vec wtvec, arma::mat mu_vec, Rcpp::List psi_inv, float gamman, float nun, arma::mat inputpt){
  
  arma::mat f_draws(s, MC, fill::zeros);
  double pi = 3.14159265;
  double con = -0.5*p;
  
  for(int m=0; m<MC; m++){
    
    arma::vec gam(n, fill::zeros);
    arma::mat evalmat(s, n, fill::zeros);
    
    for(int i=0; i<n; ++i){
      float rg = wtvec(i);
      gam(i) = random_gamma(rg);
    }
    
    arma::vec pi_vec = gam/accu(gam);
    
    for(int i=0; i<n; ++i){
      
      arma::mat psimat = psi_inv(i);
      arma::mat K = wishrnd(psimat, gamman);
      arma::mat R = chol(K);
      arma::mat R1 = inv(R);
      arma::mat sig = R1*(R1.t());
      float detR = det(R);
      arma::mat sig1 = sig/nun;
      
      arma::vec meanrow = ((mu_vec).row(i)).t();
      arma::mat theta  = mvnrnd(meanrow, sig1, 1);
      arma::vec locvec = theta.col(0);
       
      for(int j=0; j<s; ++j){
        
        arma::vec dummyvec = (inputpt.row(j)).t();
        arma::vec diffvec = dummyvec - locvec;
        arma::mat expr(1, 1, fill::zeros);
        
        expr = (diffvec.t())*(K*diffvec);
        
        evalmat(j,i) = detR*exp(-0.5*expr(0,0))*pow((2*pi), con);
      
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





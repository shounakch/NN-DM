#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

// [[Rcpp::export]]
mat adj_mat(mat grand_nbd) {
  int nrow = grand_nbd.n_rows, ncol = grand_nbd.n_rows;
  mat out(nrow,ncol,fill::zeros);
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      rowvec V = grand_nbd.row(j);
      vec V1 = V.t();
      uvec indices = find(V1==(i+1));
      if(indices.n_elem>=1){
        out(i,j)=1;
      }
      else{
        out(i,j)=0;
      }
    }
  }
  return out;
}

// [[Rcpp::export]]

cube Tmats_univ(int n, int L, int k, mat all_nbd, vec x, mat adj_mat, vec phi0seq, vec grand_means, vec grand_vars, float mu0, float nu0, float gamma0){
  
  cube all_T(n, n, L, fill::zeros);
  float nun = nu0 + k;
  float pwt = nu0/nun;
  float dwt = k/nun;
  float gamman = gamma0 + k;
  float df = gamman;
  float phi0sq = 1;
  float psi0 = 1;
  float psi = 1;
  float lambda = 1;
  
  double pi = 3.14159265;
  
  for(int i=0; i<n; ++i){
    
    float query = x(i);
    
    for(int j=0; j<n; ++j){
      
      float meanvec = 0;
      float varmat = 0;
      
      float meanj = grand_means(j);
      float Sj = grand_vars(j);
      float d  = 0;
      float diffvec = 0;
      
      if(adj_mat(i,j)==1){
        
        meanvec = (((k+1)*meanj)-query)/k;
        diffvec = meanj - query;
        varmat = Sj - ((k+1)*pow(diffvec, 2)/k);
        
      }
      else{
        
        d = x(all_nbd(j,k)-1);
        
        meanvec = (((k+1)*meanj)-d)/k;
        diffvec = meanj - d;
        varmat = Sj - ((k+1)*pow(diffvec, 2)/k);
        
      }
      
      float locvec = (pwt*mu0)+(dwt*meanvec);
      float diffvec2 = meanvec - mu0;
      float rest = varmat + (k*nu0*pow(diffvec2, 2)/nun);
      
      for(int u=0; u<L; ++u){
        
        if(i==j){
          all_T(i,j,u) = 0;
        }
        else{
          
          phi0sq = phi0seq(u);
          psi0 = gamma0*phi0sq;
          
          // Be careful about choice of psi0!
          
          psi = psi0 + rest;
          lambda = sqrt(((nun+1)/(nun*gamman))*psi);
          
          float qf = pow((query-locvec), 2)/pow(lambda, 2);
          float coef1 = exp(lgamma(0.5*(df+1))-lgamma(0.5*df))/(pow(df, 1/2)*pow(pi, 1/2));
          
          float expr = 1/((pow((1+(qf/df)), (df+1)/2))*lambda);
          
          all_T(i,j,u) = expr*coef1;
          
        }
        
        // Rcpp::Rcout << "i: " << i << std::endl;
        
      }
      
    }
    
  }
  
  return all_T;
  
}
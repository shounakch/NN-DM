//Leave One Out Functions.cpp
#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

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

cube Tmats(int n, int p, int L, int k, mat all_nbd, mat x, mat adj_mat, vec phi0seq, mat grand_means, cube grand_vars, vec mu0, float nu0, float gamma0){
  
  cube all_T(n, n, L, fill::zeros);
  float nun = nu0 + k;
  float pwt = nu0/nun;
  float dwt = k/nun;
  float gamman = gamma0 + k;
  float df = gamman - p + 1;
  double pi = 3.14159265;
  
  for(int i=0; i<n; ++i){
    
    vec query = (x.row(i)).t();
    
    for(int j=0; j<n; ++j){
      
      vec meanvec(p, fill::zeros);
      mat varmat(p, p, fill::zeros);
      
      vec meanj = (grand_means.row(j)).t();
      mat Sj = (grand_vars.slice(j));
      rowvec d(p,fill::zeros);
      vec D(p,fill::zeros);
      
      if(adj_mat(i,j)==1){
        
        meanvec = (((k+1)*meanj)-query)/k;
        vec diffvec = meanj - query;
        varmat = Sj - ((k+1)*diffvec*diffvec.t()/k);
        
      }
      else{
        
        rowvec d = x.row((all_nbd(j,k)-1));
        vec D = d.t();
        
        meanvec = (((k+1)*meanj)-D)/k;
        vec diffvec = meanj - D;
        varmat = Sj - ((k+1)*diffvec*diffvec.t()/k);
        
      }
      
      vec locvec = (pwt*mu0)+(dwt*meanvec);
      vec diffvec2 = meanvec - mu0;
      mat rest = varmat + (k*nu0*diffvec2*diffvec2.t()/nun);
      
      for(int u=0; u<L; ++u){
        
        if(i==j){
          all_T(i,j,u) = 0;
        }
        else{
          
          float phi0sq = phi0seq(u);
          mat psi0 = eye(p, p)*((gamma0-p+1)*phi0sq);
          
          // Be careful about choice of psi0!
          
          mat psi = psi0 + rest;
          mat lambda = ((nun+1)/(nun*(gamman-p+1)))*psi;
          
          mat R = chol(lambda);
          mat R1 = inv(R);
          float detR = det(R);
          vec maindiff = query - locvec;
          mat lambinv = R1*R1.t();
          
          float qf = as_scalar(maindiff.t()*lambinv*maindiff);
          float coef1 = exp(lgamma(0.5*(df+p))-lgamma(0.5*df))/(pow(df, p/2)*pow(pi, p/2));
          
          float expr = 1/((pow((1+(qf/df)), (df+p)/2))*detR);
          
          all_T(i,j,u) = expr*coef1;
          
        }
        
      }
      
    }
    
  }
  
  return all_T;
  
}

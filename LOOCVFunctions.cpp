//Leave One Out Functions.cpp
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <Rcpp.h>
#include <mvt.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppDist)]]

using namespace arma;

// [[Rcpp::export]]
mat adj_mat(mat grand_nbd) {
  int nrow = grand_nbd.n_rows, ncol = grand_nbd.n_rows;
  mat out(nrow,ncol,fill::zeros);
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
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
vec f1(vec phi0seq, mat x, Rcpp::List newloc, Rcpp::List psi_rest, int k, float gamma0, float nu0){
  
  int n = x.n_rows, p = x.n_cols, L = phi0seq.n_elem;
  vec crit(L,fill::zeros);
  float gamman = k + gamma0;
  
  for(int l = 0; l < L; l++){
    
    float phi0sq = phi0seq(l);
    float dfn = gamman - p + 1;
    float coef = (nu0+k+1)/((nu0+k)*(gamman-p+1));
    
    mat psi0 = eye(p, p)*(gamma0*phi0sq);
    mat Tmat(n, n, fill::zeros);
    
    for(int i = 0; i<n; i++){
      for(int j = 0; j<n; j++){
        
        mat m1 = newloc[i];
        rowvec locn1 = m1.row(j);
        vec locn = locn1.t();
        cube m2 = psi_rest[i];
        mat newbw = (psi0 + (m2.slice(j)))*coef;
        Tmat(i,j) = dmvt(x.row(i), locn, newbw, dfn, false)(0);
        
      }
    }
    
    Tmat.diag().zeros();
    
    crit(l) = accu(log(sum(Tmat,1)/(n-1)))/n;
    
  }
  
  return crit;
  
}

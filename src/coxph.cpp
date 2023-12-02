#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector p1d(NumericVector ee, NumericVector ff, IntegerVector delta, NumericMatrix Z){
  int n = delta.size();
  int m = Z.ncol();
  
  NumericVector gg(m);
  NumericVector d(m);
  
  for(int i = n-1; i >= 0; i--){
    
    for(int j = 0; j < m; j++){
      gg(j) += Z(i, j) * ee(i);
    }
    
    if(delta(i) == 1){
      for(int j = 0; j < m; j++){
        d(j) += Z(i, j) - gg(j) / ff(i);
      }
    }
    
  }
  
  return d;
}

// [[Rcpp::export]]
NumericMatrix p2d(NumericVector ee, NumericVector ff, IntegerVector delta, NumericMatrix Z){
  int n = delta.size();
  int m = Z.ncol();
  
  NumericMatrix d(m, m);
  NumericVector gg(m);
  NumericMatrix hh(m, m);
  
  for(int i = n-1; i >= 0; i--){
    
    for(int j = 0; j < m; j++){
      gg(j) += Z(i, j) * ee(i);
      for(int k = 0; k <= j; k++){
        hh(j, k) += Z(i, j) * Z(i, k) * ee(i);
        
        if(delta(i) == 1){
          d(j, k) += gg(j) * gg(k) / ff(i) / ff(i) - hh(j, k) / ff(i);
        }
      }
    }
  }
  
  for(int j = 0; j < m; j++){
    for(int k = 0; k < j; k++){
      d(k, j) = d(j, k);
    }
  }
  
  return d;
}

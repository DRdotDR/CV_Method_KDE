#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector J_cpp(
    double hN,
    NumericMatrix DZ,
    NumericVector W1,
    NumericVector W2,
    NumericMatrix V,
    NumericVector A1,
    NumericVector A2,
    double sigma
) {
    int L = W1.size();

    NumericMatrix Kmat(L, L);

    double inv_const = 1.0 / (sqrt(2.0 * M_PI) * sigma);

    // Gaussian kernel matrix
    for(int i = 0; i < L; i++) {
        for(int j = 0; j < L; j++) {

            double x = DZ(i, j) / hN / sigma;

            Kmat(i, j) = inv_const * exp(-0.5 * x * x);
        }
    }

    double S1 = 0.0;
    double S2 = 0.0;

    double LL1 = (double)L * (L - 1);

    for(int l = 0; l < L; l++) {

        double VL11 =
            (V(0,0) * L - W1[l] * W1[l]) / (L - 1);

        double VL22 =
            (V(1,1) * L - W2[l] * W2[l]) / (L - 1);

        double VL12 =
            (V(0,1) * L - W1[l] * W2[l]) / (L - 1);

        // Exact inverse of 2x2 matrix
        double detV = VL11 * VL22 - VL12 * VL12;

        double i11 =  VL22 / detV;
        double i22 =  VL11 / detV;
        double i12 = -VL12 / detV;

        double s1 = 0.0;
        double s2 = 0.0;

        for(int j = 0; j < L; j++) {

            if(j == l)
                continue;

            double ALMinus1 =
                i11 * W1[j] + i12 * W2[j];

            double ALMinus2 =
                i12 * W1[j] + i22 * W2[j];

            double kval = Kmat(l, j);

            s1 += ALMinus1 * kval;
            s2 += ALMinus2 * kval;
        }

        double coef1 = A1[l] / (LL1 * hN);
        double coef2 = A2[l] / (LL1 * hN);

        S1 += coef1 * s1;
        S2 += coef2 * s2;
    }

    return NumericVector::create(S1, S2);
}

// [[Rcpp::export]]
double Integ1_cpp_opt(double hN, NumericMatrix DZ2, NumericVector A1, double sigma) {
  int L = A1.size();
  double denom_exp = 4.0 * hN * hN * sigma * sigma;

  double diagSum = 0.0;
  double offDiagSum = 0.0;
  
  for(int i = 0; i < L; i++) {
    
    // diagonal
    diagSum += A1[i] * A1[i];
    
    // upper triangle only
    for(int j = i + 1; j < L; j++) {
      
      double val =
        A1[i] * A1[j] *
        exp(-DZ2(i, j) / denom_exp);
      
      offDiagSum += val;
    }
  }
  
  double sumE = diagSum + 2.0 * offDiagSum;
  
  double denom_integ = 2.0 * sqrt(M_PI) * (double)L * (double)L * hN * sigma;
  return sumE / denom_integ;
}

// [[Rcpp::export]]
double Integ2_cpp_opt(double hN, NumericMatrix DZ2, NumericVector A2, double sigma) {
  int L = A2.size();
  double denom_exp = 4.0 * hN * hN * sigma * sigma;
  
  double diagSum = 0.0;
  double offDiagSum = 0.0;
  
  for(int i = 0; i < L; i++) {
    
    // diagonal
    diagSum += A2[i] * A2[i];
    
    // upper triangle only
    for(int j = i + 1; j < L; j++) {
      
      double val =
        A2[i] * A2[j] *
        exp(-DZ2(i, j) / denom_exp);
      
      offDiagSum += val;
    }
  }
  
  double sumE = diagSum + 2.0 * offDiagSum;
  
  
  double denom_integ = 2.0 * sqrt(M_PI) * (double)L * (double)L * hN * sigma;
  return sumE / denom_integ;
}

// [[Rcpp::export]]
double CV1_cpp(double hN, NumericMatrix DZ, NumericMatrix DZ2, NumericVector W1, 
               NumericVector W2, NumericMatrix V, NumericVector A1, 
               NumericVector A2, double sigma) {
    
    NumericVector Jval = J_cpp(hN, DZ, W1, W2, V, A1, A2, sigma);
    double integ1 = Integ1_cpp_opt(hN, DZ2, A1, sigma);
    
    return -2.0 * Jval[0] + integ1; 
}

// [[Rcpp::export]]
double CV2_cpp(double hN, NumericMatrix DZ, NumericMatrix DZ2, NumericVector W1, 
               NumericVector W2, NumericMatrix V, NumericVector A1, 
               NumericVector A2, double sigma) {
    
    NumericVector Jval = J_cpp(hN, DZ, W1, W2, V, A1, A2, sigma);
    double integ2 = Integ2_cpp_opt(hN, DZ2, A2, sigma);
    
    return -2.0 * Jval[1] + integ2; 
}

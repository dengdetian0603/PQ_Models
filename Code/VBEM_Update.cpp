#include <Rcpp.h>
#include <iostream>
#include <vector>

using namespace std;
using namespace Rcpp;




// [[Rcpp::export]]
NumericMatrix loopyUpdateQD(NumericMatrix& qL, NumericMatrix qD,
                            NumericMatrix qLD, NumericMatrix& qLL,
                            NumericMatrix qL2D2,
                            NumericVector& muTheta, NumericVector& tauTheta,
                            double muRho, double tauRho, double psiAB,
                            int nLoop) {
  int n = qL.nrow();
  int K = qL.ncol();
  NumericMatrix varRhoL(n, K);
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < K; ++k) {
      varRhoL(i, k) = (muRho * muRho + 1/tauRho) * qL(i, k) -
                      muRho * muRho * qL(i, k) * qL(i, k);
    }
  }
  double qLDi_k1k2;
  double qLDi_k2k1;
  double qL2D2i_k1k2; 
  double qL2D2i_k2k1;
  double varSi_k1k2;
  double varSi_k2k1;
  double e0i_k1k2;
  double e0i_k2k1;
  double e1i_k1k2;
  double e1i_k2k1;
  double T0i;
  double T1i;
  for (int loop = 0; loop < nLoop; ++loop) {
    for (int k1 = 0; k1 < K; ++k1) {
      for (int k2 = k1 + 1; k2 < K; ++k2) {
        double T0sum = 0;
        double T1sum = 0;
        for (int i = 0; i < n; ++i) {
          qLDi_k1k2 = qLD(i, k1) - qL(i, k2) * qD(k1, k2);
          qLDi_k2k1 = qLD(i, k2) - qL(i, k1) * qD(k1, k2);
          qL2D2i_k1k2 = qL2D2(i, k1) - qL(i, k2) * qL(i, k2) * qD(k1, k2) * qD(k1, k2);
          qL2D2i_k2k1 = qL2D2(i, k2) - qL(i, k1) * qD(k1, k2) * qL(i, k1) * qD(k1, k2);
          varSi_k1k2 = 1/tauTheta[k1] + (1/tauRho + muRho * muRho) * 
                       (qLDi_k1k2 - qL2D2i_k1k2) + pow(qLDi_k1k2, 2.0)/tauRho;
          varSi_k2k1 = 1/tauTheta[k2] + (1/tauRho + muRho * muRho) * 
                       (qLDi_k2k1 - qL2D2i_k2k1) + pow(qLDi_k2k1, 2.0)/tauRho;
          e0i_k1k2 = exp(muTheta[k1] + muRho * qLDi_k1k2);
          e0i_k2k1 = exp(muTheta[k2] + muRho * qLDi_k2k1);        
          e1i_k1k2 = exp(muTheta[k1] + muRho * qLDi_k1k2 + muRho * qL(i, k2));
          e1i_k2k1 = exp(muTheta[k2] + muRho * qLDi_k2k1 + muRho * qL(i, k1));
          T0i = log(1 + e0i_k1k2) + e0i_k1k2/(2 * pow(1 + e0i_k1k2, 2.0)) * varSi_k1k2 +
                log(1 + e0i_k2k1) + e0i_k2k1/(2 * pow(1 + e0i_k2k1, 2.0)) * varSi_k2k1;
          T1i = log(1 + e1i_k1k2) + e1i_k1k2/(2 * pow(1 + e1i_k1k2, 2.0)) * (varSi_k1k2 + varRhoL(i, k2)) +
                log(1 + e1i_k2k1) + e1i_k2k1/(2 * pow(1 + e1i_k2k1, 2.0)) * (varSi_k2k1 + varRhoL(i, k1));
          T0sum += T0i;
          T1sum += T1i;
          // cout << i << "," << k1 << "," << k2 << "," << qLDi_k1k2 << "/";
        }
        // cout << "," << k1 << "," << k2 << "," << T1sum << "," << T0sum << "/";
        //double D0 = exp(-T0sum);
        //double D1 = exp(2 * muRho * qLL(k1, k2) - T1sum + psiAB);
        double qD_k1k2 = 1/(1 + exp(-T0sum - (2 * muRho * qLL(k1, k2) - T1sum + psiAB)));
        for (int i = 0; i < n; ++i) {
          qLD(i, k1) = qLD(i, k1) + qL(i, k2) * (qD_k1k2 - qD(k1, k2));
          qLD(i, k2) = qLD(i, k2) + qL(i, k1) * (qD_k1k2 - qD(k1, k2));
          qL2D2(i, k1) = qL2D2(i, k1) + pow(qL(i, k2), 2) *
                         (qD_k1k2 * qD_k1k2 - qD(k1, k2) * qD(k1, k2));
          qL2D2(i, k2) = qL2D2(i, k2) + pow(qL(i, k1), 2) *
                         (qD_k1k2 * qD_k1k2 - qD(k1, k2) * qD(k1, k2));
        }
        qD(k1, k2) = qD_k1k2;
        qD(k2, k1) = qD_k1k2;
      }
    }
  }
  return qD;
}


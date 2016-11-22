#include <Rcpp.h>
#include <iostream>
#include <vector>

using namespace std;
using namespace Rcpp;

NumericVector log_Prob_MSSBSigivenL(const int& K, IntegerMatrix& Lall_matrix,
                                    NumericVector& MSSi, NumericVector& MBSi,
                                    NumericVector& ss_tpr,
                                    NumericVector& bs_tpr,
                                    NumericVector& bs_fpr) {
    int J = Lall_matrix.nrow();
    NumericVector result(J);
    for(int j = 0; j < J; ++j) {
        IntegerVector Lj = Lall_matrix(j,_);
        double log_probs = 0.0;
        for(int k = 0; k < K; ++k) {
            double logP1 = 0.0;
            if (Lj[k] + MSSi[k] > 0.5) {
                logP1 = MSSi[k] * (log(Lj[k]) + Lj[k] * log(ss_tpr[k])) +
                        Lj[k] * (1 - MSSi[k]) * log(1 - ss_tpr[k]);
            }
            double logP2 = MBSi[k] * (Lj[k] * log(bs_tpr[k]) +
                           (1 - Lj[k]) * log(bs_fpr[k])) +
                           (1 - MBSi[k]) * (Lj[k] * log(1 - bs_tpr[k]) +
                           (1 - Lj[k]) * log(1 - bs_fpr[k]));
            log_probs += logP1 + logP2;
        }
        result[j] = log_probs;
    }
    return result;
}

// [[Rcpp::export]]
NumericMatrix predictLogProbMat(const int& K, IntegerMatrix& Lall_matrix,
                                NumericVector& mssi, NumericVector& mbsi,
                                NumericMatrix& post_log_cell_prob,
                                NumericMatrix& post_ss_tpr,
                                NumericMatrix& post_bs_tpr,
                                NumericMatrix& post_bs_fpr) {
    int J = Lall_matrix.nrow();
    int N = post_log_cell_prob.nrow();

    NumericMatrix pred(N, J);
    for(int i = 0; i < N; ++i){
        NumericVector ss_tpr = post_ss_tpr(i,_);
        NumericVector bs_tpr = post_bs_tpr(i,_);
        NumericVector bs_fpr = post_bs_fpr(i,_);
        NumericVector logProb1 = log_Prob_MSSBSigivenL(K, Lall_matrix,
                                                       mssi, mbsi, 
                                                       ss_tpr,
                                                       bs_tpr,
                                                       bs_fpr);
        pred(i,_) =  logProb1 + post_log_cell_prob(i,_);
    }
    return pred;
}



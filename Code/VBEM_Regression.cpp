#include <iostream>
#include <RcppArmadillo.h>
#include <vector>
#include </Users/dengdetian0603/Downloads/dlib-18.18/dlib/optimization.h>
// #include </users/ddeng/dlib-18.18/dlib/optimization.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

typedef dlib::matrix<double,0,1> column_vector;

class QBeta_pk
{
public:
  typedef::column_vector column_vector;
  typedef dlib::matrix<double> general_matrix;
  QBeta_pk(const int& p, const int& k, const mat& X, const mat& XX,
           const mat& qL, const mat& betaPostMean, const mat& betaPostTau,
           const mat& rhoQLD, const mat& varRhoR,
           const mat& betaPriorMean, const mat& betaPriorTau) {
    this->p = p;
    this->k = k;
    this->X = X;
    this->XX = XX;
    this->qL = qL;
    this->betaPostMean = betaPostMean;
    this->betaPostTau = betaPostTau;
    this->betaPostVar = 1/betaPostTau;
    this->rhoQLD = rhoQLD;
    this->varRhoR = varRhoR;
    this->betaPriorMean = betaPriorMean;
    this->betaPriorTau = betaPriorTau;
    
    this->betaPostMean(p, k) = 0.0;
    this->betaPostVar(p, k) = 0.0;
  }
  
  double operator() (const column_vector& beta_pk) const {
    double x = beta_pk(0);
    vec eVec = exp(this->X.col(this->p) * x + 
                   this->X * this->betaPostMean.col(this->k) +
                   this->rhoQLD.col(this->k));
    double qVal_1 = sum(this->qL.col(this->k) % X.col(this->p)) * x;
    vec qvarVec = this->varRhoR.col(this->k) +
                  this->XX * this->betaPostVar.col(this->k); 
    vec qVec = -log(1 + eVec) - eVec/(2 * square(1 + eVec)) % qvarVec;
    double qVal_2 = sum(qVec) - 0.5 * pow(x - this->betaPriorMean(p, k), 2) * 
                    this->betaPriorTau(p, k);
    return qVal_1 + qVal_2;
  }
  
  void get_derivative_and_hessian(const column_vector& beta_pk,
                                  column_vector& grad,
                                  general_matrix& hess) const {
  }

private:
  int p;
  int k;
  mat X;
  mat XX;
  mat qL;
  mat betaPostMean;
  mat betaPostTau;
  mat betaPostVar;
  mat rhoQLD;
  mat varRhoR;
  mat betaPriorMean;
  mat betaPriorTau;
};

class QBeta_gradient
{
public:
  typedef::column_vector column_vector;
  typedef dlib::matrix<double> general_matrix;
  QBeta_gradient(const int& p, const int& k, const mat& X, const mat& XX,
                 const mat& qL, const mat& betaPostMean, const mat& betaPostTau,
                 const mat& rhoQLD, const mat& varRhoR,
                 const mat& betaPriorMean, const mat& betaPriorTau) {
      this->p = p;
      this->k = k;
      this->X = X;
      this->XX = XX;
      this->qL = qL;
      this->betaPostMean = betaPostMean;
      this->betaPostTau = betaPostTau;
      this->betaPostVar = 1/betaPostTau;
      this->rhoQLD = rhoQLD;
      this->varRhoR = varRhoR;
      this->betaPriorMean = betaPriorMean;
      this->betaPriorTau = betaPriorTau;
      
      this->betaPostMean(p, k) = 0.0;
      this->betaPostVar(p, k) = 0.0;
    }
  
  column_vector operator() (const column_vector& beta_pk) const {
    double x = beta_pk(0);
    vec eVec = exp(this->X.col(this->p) * x + 
               this->X * this->betaPostMean.col(this->k) +
               this->rhoQLD.col(this->k));
    double qGradVal_1 = sum(this->qL.col(this->k) % X.col(this->p));
    vec qvarVec = this->varRhoR.col(this->k) +
                  this->XX * this->betaPostVar.col(this->k); 

    vec qVec= - (eVec % this->X.col(this->p))/(1 + eVec) -
              (this->X.col(this->p) % (eVec - pow(eVec, 3)))/pow(1 + eVec, 4) %
              qvarVec * 0.5;

    double qGradVal_2 = sum(qVec) - (x - this->betaPriorMean(p, k)) * 
                        this->betaPriorTau(p, k);
    column_vector gradVal(1);
    gradVal(0) = qGradVal_1 + qGradVal_2;
    return gradVal;
  }

private:
  int p;
  int k;
  mat X;
  mat XX;
  mat qL;
  mat betaPostMean;
  mat betaPostTau;
  mat betaPostVar;
  mat rhoQLD;
  mat varRhoR;
  mat betaPriorMean;
  mat betaPriorTau;
};

class QRho
{
public:
  typedef::column_vector column_vector;
  QRho(const mat& qL, const mat& qLD, const mat& qLLD,
       const mat& thetaMat, const mat& thetaVarMat,
       const double& rhoPriorMean, const double& rhoPriorTau) {
    this->qL = qL;
    this->qLD = qLD;
    this->qLLD = qLLD;
    this->thetaMat = thetaMat;
    this->thetaVarMat = thetaVarMat;
    this->rhoPriorMean = rhoPriorMean;
    this->rhoPriorTau = rhoPriorTau;
    this->dPois = zeros<cube>(qL.n_rows, qL.n_cols, qL.n_cols);
    for (int j = 0; j < this->qL.n_cols; ++j) {
      for (int i = 0; i < this->qL.n_rows; ++i) {
        for (int k = 0; k < this->qL.n_cols; ++k) {
          this->dPois(i, k, j) = R::dpois(j, this->qLD(i, k), false);
        }
      }
    }
  }
  
  double operator() (const column_vector& rho) const {
    double x = rho(0);
    mat qMat(size(this->qL));
    qMat.fill(0);
    mat eMat(size(this->thetaMat));
    mat fMat(size(this->thetaMat));
    for (int j = 0; j < qL.n_cols; ++j) {
      eMat = exp(this->thetaMat + j * x);
      fMat = log(1 + eMat) + 0.5 * (eMat % this->thetaVarMat)/square(1 + eMat);
      qMat += this->dPois.slice(j) % fMat;
    }
    double qVal = accu(this->qLLD) * x - accu(qMat) - 
                  0.5 * this->rhoPriorTau * pow(x - this->rhoPriorMean, 2);
    return qVal;
  }
  
private:
  mat qL;
  mat qLD;
  mat qLLD;
  mat thetaMat;
  mat thetaVarMat;
  mat betaPostVar;
  double rhoPriorMean;
  double rhoPriorTau;
  cube dPois;
};

class QRho_gradient
{
public:
  typedef::column_vector column_vector;
  QRho_gradient(const mat& qL, const mat& qLD, const mat& qLLD,
                const mat& thetaMat, const mat& thetaVarMat,
                const double& rhoPriorMean, const double& rhoPriorTau) {
    this->qL = qL;
    this->qLD = qLD;
    this->qLLD = qLLD;
    this->thetaMat = thetaMat;
    this->thetaVarMat = thetaVarMat;
    this->rhoPriorMean = rhoPriorMean;
    this->rhoPriorTau = rhoPriorTau;
    this->dPois = zeros<cube>(qL.n_rows, qL.n_cols, qL.n_cols);
    for (int j = 0; j < this->qL.n_cols; ++j) {
      for (int i = 0; i < this->qL.n_rows; ++i) {
        for (int k = 0; k < this->qL.n_cols; ++k) {
          this->dPois(i, k, j) = R::dpois(j, this->qLD(i, k), false);
        }
      }
    }
  }
  
  column_vector operator() (const column_vector& rho) const {
    double x = rho(0);
    mat qMat(size(this->qL));
    qMat.fill(0);
    mat eMat(size(this->thetaMat));
    mat fMat(size(this->thetaMat));
    for (int j = 0; j < qL.n_cols; ++j) {
      eMat = exp(this->thetaMat + j * x);
      fMat = (eMat * j)/(1 + eMat) +
             0.5 * j * (eMat - pow(eMat, 3)) %
             this->thetaVarMat/pow(1 + eMat, 4);
      qMat += this->dPois.slice(j) % fMat;
    }
    double gradVal = accu(this->qLLD) - accu(qMat) - 
                     this->rhoPriorTau * (x - this->rhoPriorMean);
    column_vector gradVec(1);
    gradVec(0) = gradVal;
    return gradVec;
  }
  
private:
  mat qL;
  mat qLD;
  mat qLLD;
  mat thetaMat;
  mat thetaVarMat;
  mat betaPostVar;
  double rhoPriorMean;
  double rhoPriorTau;
  cube dPois;
};
/* ========================================================================== */

void loopyUpdateQL(mat& qL, mat& qD, mat& H, double rhoPostMean, int nLoop = 4) {
  int K = qL.n_cols;
  mat qLD = qL * qD;
  vec qL_k, diff_qL_k;
  for (int loop = 0; loop < nLoop; ++loop) {
    for (int k = 0; k < K; ++k) {
      qL_k = qL.col(k);
      qL.col(k) = 1/(1 + exp(-(H.col(k) + qLD.col(k) * rhoPostMean)));
      diff_qL_k = qL.col(k) - qL_k;
      qLD += diff_qL_k * qD.row(k);
    }
  }
}

void updateBeta(mat& betaPostMean, mat& betaPostTau,
                mat& qL, mat& qD, mat& qLD, mat& X, mat& XX,
                double rhoPostMean, double rhoPostTau, 
                mat& betaPriorMean, mat& betaPriorTau) {
  mat qL2D2 = square(qL) * square(qD);
  mat varRhoR = (rhoPostMean * rhoPostMean + 1/rhoPostTau) * (qLD - qL2D2) +
                1/rhoPostTau * square(qLD);
  mat rhoQLD = rhoPostMean * qLD;
  int K = qL.n_cols;
  int P = X.n_cols;
  for (int k = 0; k < K; ++k) {
    for (int p = 0; p < P; ++p) {
      column_vector starting_point(1);
      starting_point(0) = betaPostMean(p, k);
      //cout << starting_point(0) << endl;
      QBeta_pk qbeta(p, k, X, XX, qL, betaPostMean, betaPostTau,
                     rhoQLD, varRhoR, betaPriorMean, betaPriorTau);
      QBeta_gradient qgrad(p, k, X, XX, qL, betaPostMean, betaPostTau,
                           rhoQLD, varRhoR, betaPriorMean, betaPriorTau);
      // dlib::find_max_using_approximate_derivatives(
      //  dlib::bfgs_search_strategy(),
      //  dlib::objective_delta_stop_strategy(1e-7),
      //  qbeta, starting_point, 100, 1e-7);
      dlib::find_max(
        dlib::bfgs_search_strategy(),
        dlib::objective_delta_stop_strategy(1e-7, 100),
        qbeta, qgrad, starting_point, 100);

      column_vector betaTau = dlib::derivative(qgrad)(starting_point);
      betaPostMean(p, k) = starting_point(0);
      betaPostTau(p, k) = abs(betaTau(0));
    }
  }
}

void updateRho(double& rhoPostMean, double& rhoPostTau,
               mat& qL, mat& qLD, mat& thetaMat, mat& thetaVarMat,
               double rhoPriorMean, double rhoPriorTau) {
  mat qLLD = qL % qLD;
  QRho qrho(qL, qLD, qLLD, thetaMat, thetaVarMat,
            rhoPriorMean, rhoPriorTau);
  QRho_gradient qrho_grad(qL, qLD, qLLD, thetaMat, thetaVarMat,
                          rhoPriorMean, rhoPriorTau);

  column_vector starting_point(1);
  starting_point(0) = rhoPostMean;
  dlib::find_max_box_constrained(
    dlib::bfgs_search_strategy(),
    dlib::objective_delta_stop_strategy(1e-7, 100),
    qrho, qrho_grad, starting_point, -20, 10);
  
  // column_vector gradVal1(1); 
  // gradVal1 = qrho_grad(starting_point);
  // double gradVal2 = dlib::derivative(qrho)(starting_point);
  // cout << starting_point(0) << "," << gradVal1(0) << "," << gradVal2 << endl;
 
  column_vector rhoTau = dlib::derivative(qrho_grad)(starting_point);
  rhoPostMean = starting_point(0);
  rhoPostTau = abs(rhoTau(0));
}

void loopyUpdateQD(mat& qL, mat& qD, mat& qLD, mat& qLL, mat& qL2D2,
                   mat& thetaMat, mat& thetaVarMat,
                   double rhoPostMean, double rhoPostTau, double psiD,
                   int nLoop) {
  int n = qL.n_rows;
  int K = qL.n_cols;
  mat varRhoL(n, K);
  varRhoL = (pow(rhoPostMean, 2) + 1/rhoPostTau) * qL -
            pow(rhoPostMean, 2) * square(qL);
 
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
          qL2D2i_k1k2 = qL2D2(i, k1) - qL(i, k2) * qL(i, k2) *
                        qD(k1, k2) * qD(k1, k2);
          qL2D2i_k2k1 = qL2D2(i, k2) - qL(i, k1) * qD(k1, k2) *
                        qL(i, k1) * qD(k1, k2);
          varSi_k1k2 = thetaVarMat(i, k1) + (1/rhoPostTau +
                       rhoPostMean * rhoPostMean) * (qLDi_k1k2 - qL2D2i_k1k2) +
                       pow(qLDi_k1k2, 2.0)/rhoPostTau;
          varSi_k2k1 = thetaVarMat(i, k2) + (1/rhoPostTau +
                       rhoPostMean * rhoPostMean) * (qLDi_k2k1 - qL2D2i_k2k1) +
                       pow(qLDi_k2k1, 2.0)/rhoPostTau;
          e0i_k1k2 = exp(thetaMat(i, k1) + rhoPostMean * qLDi_k1k2);
          e0i_k2k1 = exp(thetaMat(i, k2) + rhoPostMean * qLDi_k2k1);        
          e1i_k1k2 = exp(thetaMat(i, k1) + rhoPostMean * qLDi_k1k2 +
                     rhoPostMean * qL(i, k2));
          e1i_k2k1 = exp(thetaMat(i, k2) + rhoPostMean * qLDi_k2k1 +
                     rhoPostMean * qL(i, k1));
          T0i = log(1 + e0i_k1k2) + e0i_k1k2/(2 * pow(1 + e0i_k1k2, 2.0)) *
                varSi_k1k2 + log(1 + e0i_k2k1) + e0i_k2k1/(2 *
                pow(1 + e0i_k2k1, 2.0)) * varSi_k2k1;
          T1i = log(1 + e1i_k1k2) + e1i_k1k2/(2 * pow(1 + e1i_k1k2, 2.0)) *
                (varSi_k1k2 + varRhoL(i, k2)) + log(1 + e1i_k2k1) +
                e1i_k2k1/(2 * pow(1 + e1i_k2k1, 2.0)) *
                (varSi_k2k1 + varRhoL(i, k1));
          T0sum += T0i;
          T1sum += T1i;
        }
        
        double qD_k1k2 = 1/(1 + exp(-T0sum -
                         (2 * rhoPostMean * qLL(k1, k2) - T1sum + psiD)));
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
}


vec digammaVec(vec A, uvec& avail) {
  vec digammaA(size(avail));
  for (int i = 0; i < avail.n_elem; ++i) {
    digammaA(i) = R::digamma(A(avail(i)));
  }
  return digammaA;
}

// [[Rcpp::export]]
List regVBEM(List inputObj, List hyperPar, List initVal,
             int maxIter, double tol, int nLoop) {
  // prepare data ==============================================================
  uvec ssAvail = inputObj["ss.available"];
  ssAvail -= 1;
  uvec bsAvail = inputObj["bs.available"];
  bsAvail -= 1;
  mat mssCase = inputObj["MSS.case"];
  mat mbsCase = inputObj["MBS.case"];
  mat mbsCtrl = inputObj["MBS.ctrl"];
  mat X = inputObj["X"];
  mat XX = square(X);
  int nCase = mbsCase.n_rows;
  int nCtrl = mbsCtrl.n_rows;
  int K = mbsCase.n_cols;
  mat mssCaseAvail;
  mat mssCaseInfty;
  if (!ssAvail.is_empty()) {
    mssCaseAvail = mssCase.cols(ssAvail);
    mssCaseInfty = mssCaseAvail * 1E+16;
  }
  mat mbsCaseAvail = mbsCase.cols(bsAvail);
  mat mbsCtrlAvail = mbsCtrl.cols(bsAvail);
  mat notMbsCase = ones<mat>(nCase, bsAvail.n_elem) - mbsCaseAvail;
  vec mbsCaseSum = sum(mbsCaseAvail, 0).t();
  vec mbsCtrlSum = sum(mbsCtrlAvail, 0).t();
  // prepare parameters --------------------------------------------------------
  vec A_st = initVal["A_st"];
  vec B_st = initVal["B_st"];
  vec A_bt = initVal["A_bt"];
  vec B_bt = initVal["B_bt"];
  vec A_bf = initVal["A_bf"];
  vec B_bf = initVal["B_bf"];
  
  mat qL = initVal["qL"];
  mat betaPostMean = initVal["mu_theta"];
  mat betaPostTau = initVal["tau_theta"];
  double rhoPostMean = initVal["mu_rho"];
  double rhoPostTau = initVal["tau_rho"];
  mat qD = initVal["qD"];
  // prepare hyper parameters --------------------------------------------------
  vec aa = hyperPar["aa"];
  vec bb = hyperPar["bb"];
  double cc = hyperPar["cc"];
  double dd = hyperPar["dd"];
  double ee = hyperPar["ee"];
  double ff = hyperPar["ff"];
  mat betaPriorMean = hyperPar["theta_mu"];
  mat betaPriorTau(size(betaPriorMean));
  betaPriorTau.fill(hyperPar["theta_tau"]);
  double rhoPriorMean = hyperPar["rho_mu"];
  double rhoPriorTau = hyperPar["rho_tau"];
  NumericVector pind_a = hyperPar["pind_a"];
  NumericVector pind_b = hyperPar["pind_b"];
  NumericVector psiD = digamma(pind_a) - digamma(pind_b);
  // prepare iteration ---------------------------------------------------------
  double parDiff = 0.9;
  int iter = 0;
  mat H(nCase, K);
  mat Hbs(nCase, bsAvail.n_elem);
  mat thetaMat = X * betaPostMean;
  mat thetaVarMat;
  vec psiAA, psiBB, psiAB, psiBAB;
  vec qSum, qmsSum, qmbSum;
  vec A_st_new, B_st_new, A_bt_new, B_bt_new, A_bf_new, B_bf_new;
  mat qLD, qLL, qL2D2, qD_prev;
  mat betaMean_prev, betaTau_prev;
  double rhoMean_prev, rhoTau_prev;
  // VB-EM iteration ===========================================================
  while (parDiff > tol && iter < maxIter) {
    // update qL (E-step) ------------------------------------------------------
    psiAA = digammaVec(A_bt, bsAvail) - digammaVec(A_bf, bsAvail);
    psiBB = digammaVec(B_bt, bsAvail) - digammaVec(B_bf, bsAvail);
    psiAB = digammaVec(A_bf + B_bf, bsAvail) -
            digammaVec(A_bt + B_bt, bsAvail);
    // psiAA = vec(digamma(as<NumericVector>(wrap(A_bt))))(bsAvail) -
    //         vec(digamma(as<NumericVector>(wrap(A_bf))))(bsAvail);
    // psiBB = vec(digamma(as<NumericVector>(wrap(B_bt))))(bsAvail) -
    //         vec(digamma(as<NumericVector>(wrap(B_bf))))(bsAvail);
    // psiAB = vec(digamma(as<NumericVector>(wrap(A_bf + B_bf))))(bsAvail) -
    //         vec(digamma(as<NumericVector>(wrap(A_bt + B_bt))))(bsAvail);
    H = thetaMat;
    Hbs = mbsCaseAvail.each_row() % psiAA.t() +
          notMbsCase.each_row() % psiBB.t();
    Hbs.each_row() += psiAB.t();
    H.cols(bsAvail) += Hbs;
    if (!ssAvail.is_empty()) {
      psiBAB = digammaVec(B_st, ssAvail) - digammaVec(A_st  + B_st, ssAvail);
      // psiBAB = vec(digamma(as<NumericVector>(wrap(B_st))))(ssAvail) -
      //          vec(digamma(as<NumericVector>(wrap(A_st + B_st))))(ssAvail);
      mat increMat = mssCaseInfty.each_row() + psiBAB.t();
      H.cols(ssAvail) += increMat;
    }
    //qL.print("before:");
    loopyUpdateQL(qL, qD, H, rhoPostMean, nLoop);
    //qL.print("after:");

    // update A, B -------------------------------------------------------------
    parDiff = 0.0;
    qSum = sum(qL, 0).t();
    // Silver standard data
    if (!ssAvail.is_empty()) {
      qmsSum = sum(qL.cols(ssAvail) % mssCaseAvail, 0).t();
      A_st_new = qmsSum + aa(ssAvail);
      B_st_new = qSum(ssAvail) - qmsSum + bb(ssAvail);
      parDiff += sum(square(A_st_new - A_st(ssAvail)) +
                     square(B_st_new - B_st(ssAvail)));
      A_st(ssAvail) = A_st_new;
      B_st(ssAvail) = B_st_new;
    }
    // Bronze standard data
    qmbSum = sum(qL.cols(bsAvail) % mbsCaseAvail, 0).t();
    A_bt_new = qmbSum + cc;
    B_bt_new = qSum(bsAvail) - qmbSum + dd;
    parDiff += sum(square(A_bt_new - A_bt(bsAvail)) +
                   square(B_bt_new - B_bt(bsAvail)));
    A_bt(bsAvail) = A_bt_new;
    B_bt(bsAvail) = B_bt_new;
    
    A_bf_new = mbsCaseSum + mbsCtrlSum - qmbSum + ee;
    B_bf_new = nCase + nCtrl + qmbSum - qSum(bsAvail) -
               mbsCaseSum - mbsCtrlSum + ff;
    parDiff += sum(square(A_bf_new - A_bf(bsAvail)) +
                   square(B_bf_new - B_bf(bsAvail)));
    A_bf(bsAvail) = A_bf_new;
    B_bf(bsAvail) = B_bf_new;
    // update beta -------------------------------------------------------------
    qLD = qL * qD;
    betaMean_prev = betaPostMean;
    betaTau_prev = betaPostTau;
    updateBeta(betaPostMean, betaPostTau,
               qL, qD, qLD, X, XX, rhoPostMean, rhoPostTau, 
               betaPriorMean, betaPriorTau);
    parDiff += accu(square(betaPostMean - betaMean_prev)) +
               accu(square(betaPostTau - betaTau_prev));

    // update rho --------------------------------------------------------------
    thetaMat = X * betaPostMean;
    thetaVarMat = XX * (1/betaPostTau);
    rhoMean_prev = rhoPostMean;
    rhoTau_prev = rhoPostTau;
    updateRho(rhoPostMean, rhoPostTau, qL, qLD, thetaMat, thetaVarMat,
              rhoPriorMean, rhoPriorTau);
    parDiff += pow(rhoPostMean - rhoMean_prev, 2) + 
               pow(rhoPostTau - rhoTau_prev, 2);
    
    // update qD ---------------------------------------------------------------
    qLL = qL.t() * qL;
    qL2D2 = square(qL) * square(qD);
    qD_prev = qD;
    loopyUpdateQD(qL, qD, qLD, qLL, qL2D2, thetaMat, thetaVarMat,
                  rhoPostMean, rhoPostTau, psiD(0), nLoop);
    parDiff += accu(square(qD - qD_prev));
    // stop/next signal --------------------------------------------------------
    parDiff = sqrt(parDiff);
    iter++;
    if (iter % 50 == 0) {
      cout << iter << ": " << parDiff << endl;
      betaPostMean.print("beta:");
    }
  }

  List parList;
  parList["A_st"] = A_st;
  parList["B_st"] = B_st;
  parList["A_bt"] = A_bt;
  parList["B_bt"] = B_bt;
  parList["A_bf"] = A_bf;
  parList["B_bf"] = B_bf;
  parList["Beta.mean"] = betaPostMean;
  parList["Beta.tau"] = betaPostTau;
  parList["Rho.mean"] = rhoPostMean;
  parList["Rho.tau"] = rhoPostTau;
  parList["qD"] = qD;
  parList["qL"] = qL;
  return parList;
}

// [[Rcpp::export]]
mat a1(mat m0, vec v0) {
  mat m1 = m0.t() * v0;
  return m1;
}


// [[Rcpp::export]]
double LpseudoExpQ(mat& parVec, mat& betaTauMat, double rhoPostTau,
                   List inputObj, List hyperPar) {
  // prepare data ==============================================================
  uvec ssAvail = inputObj["ss.available"];
  ssAvail -= 1;
  uvec bsAvail = inputObj["bs.available"];
  bsAvail -= 1;
  mat mssCase = inputObj["MSS.case"];
  mat mbsCase = inputObj["MBS.case"];
  mat mbsCtrl = inputObj["MBS.ctrl"];
  mat X = inputObj["X"];
  mat XX = square(X);
  int nCase = mbsCase.n_rows;
  int nCtrl = mbsCtrl.n_rows;
  int K = mbsCase.n_cols;
  mat mssCaseAvail;
  if (!ssAvail.is_empty()) {
    mssCaseAvail = mssCase.cols(ssAvail);
  }
  mat mbsCaseAvail = mbsCase.cols(bsAvail);
  mat mbsCtrlAvail = mbsCtrl.cols(bsAvail);
  vec mbsCaseSum = sum(mbsCaseAvail, 0).t();
  vec mbsCtrlSum = sum(mbsCtrlAvail, 0).t();
  // re-arrange parameters =====================================================
  int thisRow = nCase * K;
  mat qL = parVec.head_rows(thisRow);
  qL.reshape(nCase, K);
  mat logSStpr, log1_SStpr;
  if (!ssAvail.is_empty()) {
    logSStpr = parVec.rows(thisRow, thisRow + ssAvail.n_elem - 1);
    thisRow += ssAvail.n_elem;
    log1_SStpr = parVec.rows(thisRow, thisRow + ssAvail.n_elem - 1);
    thisRow += ssAvail.n_elem;
  }
  mat logBStpr = parVec.rows(thisRow, thisRow + bsAvail.n_elem - 1);
  thisRow += bsAvail.n_elem;
  mat log1_BStpr = parVec.rows(thisRow, thisRow + bsAvail.n_elem - 1);
  thisRow += bsAvail.n_elem;
  mat logBSfpr = parVec.rows(thisRow, thisRow + bsAvail.n_elem - 1);
  thisRow += bsAvail.n_elem;
  mat log1_BSfpr = parVec.rows(thisRow, thisRow + bsAvail.n_elem - 1);
  thisRow += bsAvail.n_elem;
  mat qD(K, K);
  qD.zeros();
  for (int r = 0; r < K; ++r) {
    for (int c = r + 1; c < K; ++c) {
      qD(r, c) = parVec(thisRow, 0);
      qD(c, r) = parVec(thisRow, 0);
      thisRow++;
    }
  }
  mat betaMat = parVec.rows(thisRow, thisRow + K * X.n_cols - 1);
  betaMat.reshape(X.n_cols, K);
  thisRow += K * X.n_cols;
  if (thisRow != parVec.n_rows - 1) {
    cout << "parVec index does not match." << endl;
  }
  double rhoPostMean = parVec(thisRow, 0);
  // prepare hyper parameters ==================================================
  vec aa = hyperPar["aa"];
  vec bb = hyperPar["bb"];
  double cc = hyperPar["cc"];
  double dd = hyperPar["dd"];
  double ee = hyperPar["ee"];
  double ff = hyperPar["ff"];
  mat betaPriorMean = hyperPar["theta_mu"];
  mat betaPriorTau(size(betaPriorMean));
  betaPriorTau.fill(hyperPar["theta_tau"]);
  double rhoPriorMean = hyperPar["rho_mu"];
  double rhoPriorTau = hyperPar["rho_tau"];
  NumericVector pind_a = hyperPar["pind_a"];
  NumericVector pind_b = hyperPar["pind_b"];
  NumericVector psiD = digamma(pind_a) - digamma(pind_b);
  // Expected log density ======================================================
  // measurement density ---------------
  vec qSum, qmsSum, qmbSum;
  vec A_st, B_st, A_bt, B_bt, A_bf, B_bf;
  mat logLSS, logLBS;
  double result = 0.0;
  qSum = sum(qL, 0).t();
  if (!ssAvail.is_empty()) {
    qmsSum = sum(qL.cols(ssAvail) % mssCaseAvail, 0).t();
    A_st = qmsSum + aa(ssAvail);
    B_st = qSum(ssAvail) - qmsSum + bb(ssAvail);
    logLSS = logSStpr.t() * A_st + log1_SStpr.t() * B_st;
    result = logLSS(0, 0);
  }
  qmbSum = sum(qL.cols(bsAvail) % mbsCaseAvail, 0).t();
  A_bt = qmbSum + cc;
  B_bt = qSum(bsAvail) - qmbSum + dd;
  A_bf = mbsCaseSum + mbsCtrlSum - qmbSum + ee;
  B_bf = nCase + nCtrl + qmbSum - qSum(bsAvail) -
         mbsCaseSum - mbsCtrlSum + ff;
  logLBS = logBStpr.t() * A_bt + log1_BStpr.t() * B_bt +
           logBSfpr.t() * A_bf + log1_BSfpr.t() * B_bf;
  result += logLBS(0, 0);
  // latent density -------------------
  mat qLD = qL * qD;
  mat thetaMat = X * betaMat;
  double qLtheta = accu(qL % thetaMat);
  double rhoLLD = rhoPostMean * accu(qL % qLD);
  cube dPois = zeros<cube>(qL.n_rows, qL.n_cols, qL.n_cols);
  for (int j = 0; j < qL.n_cols; ++j) {
    for (int i = 0; i < qL.n_rows; ++i) {
      for (int k = 0; k < qL.n_cols; ++k) {
        dPois(i, k, j) = R::dpois(j, qLD(i, k), false);
      }
    }
  }
  mat thetaVarMat = XX * (1/betaTauMat);
  mat qMat(size(qL));
  qMat.fill(0);
  mat eMat(size(thetaMat));
  mat fMat(size(thetaMat));
  for (int j = 0; j < qL.n_cols; ++j) {
    eMat = exp(thetaMat + j * rhoPostMean);
    fMat = log(1 + eMat) + 0.5 * (eMat % (thetaVarMat + j * j/rhoPostTau))/square(1 + eMat);
    qMat += dPois.slice(j) % fMat;
  }
  double logLatent = qLtheta + rhoLLD - accu(qMat);
  double logPrior = -0.5 * rhoPriorTau *
                    (1/rhoPostTau + rhoPostMean * rhoPostMean -
                    2 * rhoPostMean * rhoPriorMean) -
                    0.5 * accu(betaPriorTau %
                    (1/betaTauMat + square(betaMat) -
                    2 * betaMat % betaPriorMean)) +
                    accu(qD) * psiD(0)/2;
  result += logLatent + logPrior;
  return result;
}



vec factorial(int K);
cube gradPois(mat& qLD, vec& factJ); 
cube hessPois(mat& qLD, vec& factJ);
List approxIntegral(mat& thetaMat, mat& thetaVarMat, double rhoPostMean, double rhoPostTau);

// [[Rcpp::export]]
List hessBlock(List vbResult, List inputObj, List hyperPar) {
  // prepare data ==============================================================
  uvec ssAvail = inputObj["ss.available"];
  ssAvail -= 1;
  uvec bsAvail = inputObj["bs.available"];
  bsAvail -= 1;
  mat mssCase = inputObj["MSS.case"];
  mat mbsCase = inputObj["MBS.case"];
  mat mbsCtrl = inputObj["MBS.ctrl"];
  mat X = inputObj["X"];
  mat XX = square(X);
  int nCase = mbsCase.n_rows;
  // int nCtrl = mbsCtrl.n_rows;
  int K = mbsCase.n_cols;
  mat mssCaseAvail;
  if (!ssAvail.is_empty()) {
    mssCaseAvail = mssCase.cols(ssAvail);
  }
  mat mbsCaseAvail = mbsCase.cols(bsAvail);
  // mat mbsCtrlAvail = mbsCtrl.cols(bsAvail);
  // vec mbsCaseSum = sum(mbsCaseAvail, 0).t();
  // vec mbsCtrlSum = sum(mbsCtrlAvail, 0).t();
  // prepare parameters =====================================================
  // vec A_st = vbResult["A_st"];
  // vec B_st = vbResult["B_st"];
  // vec A_bt = vbResult["A_bt"]; 
  // vec B_bt = vbResult["B_bt"];
  // vec A_bf = vbResult["A_bf"];
  // vec B_bf = vbResult["B_bf"];
  mat qL = vbResult["qL"];
  // mat logSStpr, log1_SStpr;
  // if (!ssAvail.is_empty()) {
  //   logSStpr = digammaVec(A_st, ssAvail) - digammaVec(A_st + B_st, ssAvail);
  //   log1_SStpr = digammaVec(B_st, ssAvail) - digammaVec(A_st + B_st, ssAvail);
  // }
  // mat logBStpr = digammaVec(A_bt, ssAvail) - digammaVec(A_bt + B_bt, ssAvail);
  // mat log1_BStpr = digammaVec(B_bt, ssAvail) - digammaVec(A_bt + B_bt, ssAvail);
  // mat logBSfpr = digammaVec(A_bf, ssAvail) - digammaVec(A_bf + B_bf, ssAvail);
  // mat log1_BSfpr = digammaVec(B_bf, ssAvail) - digammaVec(A_bf + B_bf, ssAvail);
  mat qD = vbResult["qD"];
  mat betaMat = vbResult["Beta.mean"];
  mat betaTauMat = vbResult["Beta.tau"];
  double rhoPostMean = vbResult["Rho.mean"];
  double rhoPostTau = vbResult["Rho.tau"];
  // prepare hyper parameters ==================================================
  // vec aa = hyperPar["aa"];
  // vec bb = hyperPar["bb"];
  // double cc = hyperPar["cc"];
  // double dd = hyperPar["dd"];
  // double ee = hyperPar["ee"];
  // double ff = hyperPar["ff"];
  // mat betaPriorMean = hyperPar["theta_mu"];
  // mat betaPriorTau(size(betaPriorMean));
  // betaPriorTau.fill(hyperPar["theta_tau"]);
  // double rhoPriorMean = hyperPar["rho_mu"];
  // double rhoPriorTau = hyperPar["rho_tau"];
  // NumericVector pind_a = hyperPar["pind_a"];
  // NumericVector pind_b = hyperPar["pind_b"];
  // NumericVector psiD = digamma(pind_a) - digamma(pind_b);
  // hessin matrix by block ====================================================
  mat qLD = qL * qD;
  mat thetaMat = X * betaMat;
  mat thetaVarMat = XX * (1/betaTauMat);
  vec factorialJ = factorial(K);
  cube gPois = gradPois(qLD, factorialJ);
  cube hPois = hessPois(qLD, factorialJ);
  List approx = approxIntegral(thetaMat, thetaVarMat, rhoPostMean, rhoPostTau);
  cube intgeral = approx["res1"];
  cube gradIntegral = approx["res2"];
  // Hzz ------------------------------
  int lenZ = qL.n_cols * qL.n_rows + 2 * ssAvail.n_elem + 4 * bsAvail.n_elem;
  mat Hzz(lenZ, lenZ);
  Hzz.fill(0.0);
  mat hQL(size(qL));
  hQL.fill(0.0);
  for (int j = 0; j < K; ++j) {
    hQL += hPois.slice(j) % intgeral.slice(j);
  }
  int r = 0; int c = 0;
  for (int k1 = 0; k1 < K; ++k1) {
    for (int i = 0; i < nCase; ++i) {
      r = k1 * nCase + i;
      /* qL, qL */
      for (int k2 = k1; k2 < K; ++k2) {
        c = k2 * nCase + i;
        Hzz(r, c) = 2 * rhoPostMean * qD(k1, k2);
        for (int k = 0; k < K; ++k) {
          Hzz(r, c) -= hQL(i, k) * qD(k1, k) * qD(k2, k); 
        }
        Hzz(c, r) = Hzz(r, c);
      }
      /* qL, log(sstpr) */
      c = K * nCase;
      if (!ssAvail.is_empty()) {
        for (int kss = 0; kss < ssAvail.n_elem; ++kss) {
          Hzz(r, c) = mssCaseAvail(i, kss);
          Hzz(c, r) = Hzz(r, c);
          Hzz(r, c + ssAvail.n_elem) = 1 - mssCaseAvail(i, kss);
          Hzz(c + ssAvail.n_elem, r) = Hzz(r, c + ssAvail.n_elem);
          c++;
        }
      }
      /* qL, log(bstpr) & log(bsfpr) */
      c = K * nCase + 2 * ssAvail.n_elem;
      for (int kbs = 0; kbs < bsAvail.n_elem; ++kbs) {
        Hzz(r, c) = mbsCaseAvail(i, kbs);
        Hzz(c, r) = Hzz(r, c);
        Hzz(r, c + bsAvail.n_elem) = 1 - mbsCaseAvail(i, kbs);
        Hzz(c + bsAvail.n_elem, r) = Hzz(r, c + bsAvail.n_elem);
        Hzz(r, c + 2 * bsAvail.n_elem) = - mbsCaseAvail(i, kbs);
        Hzz(c + 2 * bsAvail.n_elem, r) = Hzz(r, c + 2 * bsAvail.n_elem);
        Hzz(r, c + 3 * bsAvail.n_elem) = mbsCaseAvail(i, kbs) - 1;
        Hzz(c + 3 * bsAvail.n_elem, r) = Hzz(r, c + 3 * bsAvail.n_elem);
        c++;
      }
    }
  }
  // Hza -----------------------------
  int lenA = (K * K - K)/2 + betaMat.n_elem + 1;
  mat Hza(lenZ, lenA);
  Hza.fill(0.0);
  mat gQL(size(qL));
  gQL.fill(0.0);
  for (int j = 0; j < K; ++j) {
    gQL += gPois.slice(j) % intgeral.slice(j);
  }
  mat hQLtheta(size(qL));
  hQLtheta.fill(0.0);
  for (int j = 0; j < K; ++j) {
    hQLtheta += gPois.slice(j) % gradIntegral.slice(j);
  }
  mat hQLrho(size(qL));
  hQLrho.fill(0.0);
  for (int j = 1; j < K; ++j) {
    hQLrho += gPois.slice(j) % gradIntegral.slice(j) * j;
  }
  hQLrho = hQLrho * qD;
  for (int k = 0; k < K; ++k) {
    for (int i = 0; i < nCase; ++i) {
      r = k * nCase + i;
      c = 0;
      /* qL, qD */
      for (int k2 = 0; k2 < K; ++k2) {
        for (int k1 = k2 + 1; k1 < K; ++k1) {
          if (k == k1) {
            Hza(r, c) = 2 * rhoPostMean * qL(i, k2) -
                        hQL(i, k2) * qL(i, k1) * qL(i, k1) + gQL(i, k2);
          } else if (k == k2) {
            Hza(r, c) = 2 * rhoPostMean * qL(i, k1) -
                        hQL(i, k1) * qL(i, k2) * qL(i, k2) + gQL(i, k1);
          } else {
            Hza(r, c) = - hQL(i, k2) * qL(i, k) * qL(i, k1) -
                        hQL(i, k1) * qL(i, k) * qL(i, k2);
          }
          c++; 
        }
      }
      /* qL, beta */
      c =  (K * K - K)/2;
      for (int kx = 0; kx < K; ++kx) {
        double hqlbeta = hQLtheta(i, k) * qD(k, kx);
        for (int px = 0; px < betaMat.n_rows; ++px) {
          if (kx == k) {
            Hza(r, c) = X(i, px);
          } else {
            Hza(r, c) = hqlbeta * X(i, px);
          }
          c++;
        }
      }
      /* qL, rho */
      c = lenA - 1;
      Hza(r, c) = 2 * qLD(i, k) - hQLrho(i, k);
    }
  }
  List hessList;
  hessList["Hzz"] = Hzz;
  hessList["Hza"] = Hza;
  return hessList;
}

// ----------------------------------------------------------------------------
vec factorial(int K) {
  vec facVec(K + 1);
  facVec(0) = 1;
  for (int j = 1; j <= K; ++j) {
    facVec(j) = j * facVec(j - 1);
  }
  return facVec;
}

cube gradPois(mat& qLD, vec& factJ) {
  cube gPois = zeros<cube>(qLD.n_rows, qLD.n_cols, qLD.n_cols);
  for (int j = 0; j < qLD.n_cols; ++j) {
    for (int i = 0; i < qLD.n_rows; ++i) {
      for (int k = 0; k < qLD.n_cols; ++k) {
        double lambdaik = qLD(i, k);
        gPois(i, k, j) = pow(lambdaik, j - 1) * 
                         exp(-lambdaik) * (j - lambdaik) / factJ(j);
      }
    }
  }
  return gPois;
}

cube hessPois(mat& qLD, vec& factJ) {
  cube hPois = zeros<cube>(qLD.n_rows, qLD.n_cols, qLD.n_cols);
  for (int j = 0; j < qLD.n_cols; ++j) {
    for (int i = 0; i < qLD.n_rows; ++i) {
      for (int k = 0; k < qLD.n_cols; ++k) {
        double lambdaik = qLD(i, k);
        hPois(i, k, j) = pow(lambdaik, j - 2) * 
                         exp(-lambdaik) / factJ(j) *
                         (j * (j - 1) - 2 * j * lambdaik + lambdaik * lambdaik);
      }
    }
  }
  return hPois;
}

List approxIntegral(mat& thetaMat, mat& thetaVarMat,
                    double rhoPostMean, double rhoPostTau) {
  cube res1 = zeros<cube>(thetaMat.n_rows, thetaMat.n_cols, thetaMat.n_cols);
  cube res2 = zeros<cube>(thetaMat.n_rows, thetaMat.n_cols, thetaMat.n_cols);
  mat eMat(size(thetaMat));
  mat eMatPlus1(size(thetaMat));
  for (int j = 0; j < thetaMat.n_cols; ++j) {
    eMat = exp(thetaMat + j * rhoPostMean);
    eMatPlus1 = 1 + eMat;
    res1.slice(j) = log(eMatPlus1) +
      0.5 * (eMat % (thetaVarMat + j * j/rhoPostTau))/square(eMatPlus1);
    res2.slice(j) = eMat/(eMatPlus1) + 0.5 * (eMat - pow(eMat, 3)) %
                    (thetaVarMat + j * j/rhoPostTau) / pow(eMatPlus1, 4);
  }
  List result;
  result["res1"] = res1;
  result["res2"] = res2;
  return result;
}
  
  
  
  
  
  
  
  
  
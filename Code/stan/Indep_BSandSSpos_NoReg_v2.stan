data {
// control data
  int<lower=0> N_ctrl;
  int<lower=0> K;
  int mbs_ctrl[N_ctrl, K];
// case data
  int<lower=0> N_case;
  int mbs_case[N_case, K];
  int mss_case[N_case, K];
// hyper parameters
  real<lower=0> ma;
  real<lower=0> mb;
  real<lower=0> aa;
  real<lower=0> bb;
  real<lower=0> cc;
  real<lower=0> dd;
  real<lower=0> ee;
  real<lower=0> ff;
}
parameters {
  vector<lower=0, upper=1>[K] mu;
  vector<lower=0, upper=1>[K] ss_tpr;
  vector<lower=0, upper=1>[K] bs_tpr;
  vector<lower=0, upper=1>[K] bs_fpr;
}
// transformed parameters {
//   real cont_mgs[N_case, K]; 
// }
model {
// ctrl likelihood
  for (i in 1:N_ctrl) {
    for (j in 1:K) {
      mbs_ctrl[i, j] ~ bernoulli(bs_fpr[j]);
    }
  }
// case likelihood
  for (i in 1:N_case) {
    for (j in 1:K) {
      real pbs_a;
      real pbs_b;
      pbs_a = bs_tpr[j] * ((1 - ss_tpr[j]) * mu[j] / (1 - ss_tpr[j] * mu[j])) ^
              (1 - mss_case[i, j]);
      pbs_b = (1 - mu[j]) * bs_fpr[j] / (1 - ss_tpr[j] * mu[j]) *
              (1 - mss_case[i, j]);        
      mbs_case[i, j] ~ bernoulli(pbs_a + pbs_b);
      mss_case[i, j] ~ bernoulli(ss_tpr[j] * mu[j]);
    }
  }
// prior
  for (j in 1:K) {
    mu[j] ~ beta(ma, mb);
    ss_tpr[j] ~ beta(aa, bb);
    bs_tpr[j] ~ beta(cc, dd);
    bs_fpr[j] ~ beta(ee, ff);
  }
}



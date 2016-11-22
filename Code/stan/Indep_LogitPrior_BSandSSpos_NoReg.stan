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
  vector[K] mu_logit;
  cov_matrix[K] V;
  real<lower=0> cc;
  real<lower=0> dd;
  real<lower=0> ee;
  real<lower=0> ff;
}
parameters {
  vector<lower=0, upper=1>[K] bs_fpr;
  vector<lower=0, upper=1>[K] bs_tpr;
  vector[K] logit_mu;
}
transformed parameters {
  vector<lower=0, upper=1>[K] mu;
  for (j in 1:K) {
    mu[j] = inv_logit(logit_mu[j]);
  }
}
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
      pbs_a = bs_tpr[j] * mu[j]^(1 - mss_case[i, j]);
      pbs_b = (1 - mu[j]) * bs_fpr[j] * (1 - mss_case[i, j]); 
      target += (mbs_case[i, j] * log(pbs_a + pbs_b) +
                (1 - mbs_case[i, j]) * log(1 - pbs_a - pbs_b));
    }
  }
// prior
  logit_mu ~ multi_normal(mu_logit, V);
  for (j in 1:K) {
    bs_tpr[j] ~ beta(cc, dd);
    bs_fpr[j] ~ beta(ee, ff);
  }
}



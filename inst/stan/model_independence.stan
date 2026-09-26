data {
  int<lower=1> N;
  int<lower=1> S;
  array[N] int<lower=1,upper=S> study;
  vector[N] X;
  vector[N] Y;
  real<lower=1> nu;

  // Hyperpriors on the generative means
  real prior_mu_beta;      real<lower=0> prior_scale_beta;
  real prior_mu_alpha;  real<lower=0> prior_scale_alpha;
  real prior_mu_sig;  real<lower=0> prior_scale_sig;

  // fixed large heterogeneities = the tau -> inf limit
  // (99.9th percentile of each prior, passed as data)
  real<lower=0> prior_tau_beta_fixed;
  real<lower=0> prior_tau_alpha_fixed;
  real<lower=0> prior_tau_sig_fixed;
}
parameters {
  real          mu_beta;
  real<lower=0> mu_alpha;
  real<lower=0> mu_sig;
  vector[S]          z_beta;
  vector<lower=0>[S] alpha;
  vector<lower=0>[S] sigma;
}
transformed parameters {
  // study slopes with tau fixed at the large constant (independence limit)
  vector[S] beta = mu_beta + prior_tau_beta_fixed * z_beta;
}
model {
  mu_beta     ~ student_t(nu, prior_mu_beta,     prior_scale_beta);
  mu_alpha ~ student_t(nu, prior_mu_alpha, prior_scale_alpha);   // auto-truncated >0
  mu_sig ~ student_t(nu, prior_mu_sig, prior_scale_sig);   // auto-truncated >0

  z_beta     ~ student_t(nu, 0, 1);
  alpha   ~ student_t(nu, mu_alpha, prior_tau_alpha_fixed);
  target += -S * student_t_lccdf(0 | nu, mu_alpha, prior_tau_alpha_fixed);
  sigma ~ student_t(nu, mu_sig, prior_tau_sig_fixed);
  target += -S * student_t_lccdf(0 | nu, mu_sig, prior_tau_sig_fixed);

  for (n in 1:N)
    Y[n] ~ normal(alpha[study[n]] + beta[study[n]] * X[n], sigma[study[n]]);
}
generated quantities {
  vector[N] y_rep;
  vector[N] log_lik_Y;
  for (n in 1:N) {
    real mu_y = alpha[study[n]] + beta[study[n]] * X[n];
    y_rep[n]     = normal_rng(mu_y, sigma[study[n]]);
    log_lik_Y[n] = normal_lpdf(Y[n] | mu_y, sigma[study[n]]);
  }
}
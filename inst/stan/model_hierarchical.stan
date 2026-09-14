data {
  int<lower=1> N;
  int<lower=1> S;
  array[N] int<lower=1,upper=S> study;
  vector[N] X;
  vector[N] M;
  real<lower=1> nu;                        // Student-t degrees of freedom

  // Hyperpriors: generative means / scales
  real prior_mu_a;         real<lower=0> prior_scale_a;         // slope beta (untruncated)
  real prior_mu_alpha;     real<lower=0> prior_scale_alpha;     // intercept alpha (>0)
  real prior_mu_sig_m;     real<lower=0> prior_scale_sig_m;     // residual sigma (>0)

  // Hyperpriors: between-study heterogeneities (tau, >0)
  real prior_mu_tau_a;     real<lower=0> prior_scale_tau_a;
  real prior_mu_tau_int1;  real<lower=0> prior_scale_tau_int1;
  real prior_mu_tau_sig_m; real<lower=0> prior_scale_tau_sig_m;
}
parameters {
  real          mu_a;                      // generative slope mean
  real<lower=0> mu_alpha;                  // generative intercept mean (>0)
  real<lower=0> mu_sig_m;                  // generative residual mean (>0)

  real<lower=0> tau_a;                     // between-study scales (>0)
  real<lower=0> tau_intercept1;
  real<lower=0> tau_sig_m;

  vector[S]          z_a;                  // non-centered slopes
  vector<lower=0>[S] alpha;                // centered positive intercepts
  vector<lower=0>[S] sigma_M;              // centered positive residual scales
}
transformed parameters {
  vector[S] a = mu_a + tau_a * z_a;        // study-specific slopes
}
model {
  // generative means
  mu_a     ~ student_t(nu, prior_mu_a,     prior_scale_a);
  mu_alpha ~ student_t(nu, prior_mu_alpha, prior_scale_alpha);   // auto-truncated >0
  mu_sig_m ~ student_t(nu, prior_mu_sig_m, prior_scale_sig_m);   // auto-truncated >0

  // between-study heterogeneities (truncated-t >0)
  tau_a          ~ student_t(nu, prior_mu_tau_a,     prior_scale_tau_a);
  tau_intercept1 ~ student_t(nu, prior_mu_tau_int1,  prior_scale_tau_int1);
  tau_sig_m      ~ student_t(nu, prior_mu_tau_sig_m, prior_scale_tau_sig_m);

  // study-level draws (truncation normalized because the scale is a parameter)
  z_a     ~ student_t(nu, 0, 1);
  alpha   ~ student_t(nu, mu_alpha, tau_intercept1);
  target += -S * student_t_lccdf(0 | nu, mu_alpha, tau_intercept1);
  sigma_M ~ student_t(nu, mu_sig_m, tau_sig_m);
  target += -S * student_t_lccdf(0 | nu, mu_sig_m, tau_sig_m);

  // likelihood
  for (n in 1:N)
    M[n] ~ normal(alpha[study[n]] + a[study[n]] * X[n], sigma_M[study[n]]);
}
generated quantities {
  vector[N] m_rep;
  vector[N] log_lik_M;
  for (n in 1:N) {
    real mu_m = alpha[study[n]] + a[study[n]] * X[n];
    m_rep[n]     = normal_rng(mu_m, sigma_M[study[n]]);
    log_lik_M[n] = normal_lpdf(M[n] | mu_m, sigma_M[study[n]]);
  }
}
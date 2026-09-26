data {
  int<lower=1> N;
  int<lower=1> S;
  array[N] int<lower=1,upper=S> study;
  vector[N] X;
  vector[N] Y;
  real<lower=1> nu;                        // Student-t degrees of freedom

  // Hyperpriors: generative means / scales
  real prior_mu_beta;         real<lower=0> prior_scale_beta;         // slope beta (untruncated)
  real prior_mu_alpha;     real<lower=0> prior_scale_alpha;     // intercept alpha (>0)
  real prior_mu_sig;     real<lower=0> prior_scale_sig;     // residual sigma (>0)

  // Hyperpriors: between-study heterogeneities (tau, >0)
  real prior_mu_tau_beta;     real<lower=0> prior_scale_tau_beta;
  real prior_mu_tau_alpha;  real<lower=0> prior_scale_tau_alpha;
  real prior_mu_tau_sig; real<lower=0> prior_scale_tau_sig;
}
parameters {
  real          mu_beta;                      // generative slope mean
  real<lower=0> mu_alpha;                  // generative intercept mean (>0)
  real<lower=0> mu_sig;                  // generative residual mean (>0)

  real<lower=0> tau_beta;                     // between-study scales (>0)
  real<lower=0> tau_alpha;
  real<lower=0> tau_sig;

  vector[S]          z_beta;                  // non-centered slopes
  vector<lower=0>[S] alpha;                // centered positive intercepts
  vector<lower=0>[S] sigma;              // centered positive residual scales
}
transformed parameters {
  vector[S] beta = mu_beta + tau_beta * z_beta;        // study-specific slopes
}
model {
  // generative means
  mu_beta     ~ student_t(nu, prior_mu_beta,     prior_scale_beta);
  mu_alpha ~ student_t(nu, prior_mu_alpha, prior_scale_alpha);   // auto-truncated >0
  mu_sig ~ student_t(nu, prior_mu_sig, prior_scale_sig);   // auto-truncated >0

  // between-study heterogeneities (truncated-t >0)
  tau_beta          ~ student_t(nu, prior_mu_tau_beta,     prior_scale_tau_beta);
  tau_alpha ~ student_t(nu, prior_mu_tau_alpha,  prior_scale_tau_alpha);
  tau_sig      ~ student_t(nu, prior_mu_tau_sig, prior_scale_tau_sig);

  // study-level draws (truncation normalized because the scale is a parameter)
  z_beta     ~ student_t(nu, 0, 1);
  alpha   ~ student_t(nu, mu_alpha, tau_alpha);
  target += -S * student_t_lccdf(0 | nu, mu_alpha, tau_alpha);
  sigma ~ student_t(nu, mu_sig, tau_sig);
  target += -S * student_t_lccdf(0 | nu, mu_sig, tau_sig);

  // likelihood
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
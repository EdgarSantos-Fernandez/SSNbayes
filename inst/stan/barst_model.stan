  data {
  int<lower = 1> T;
  int<lower = 1> niter;
  int<lower = 1> iter_num;
    int<lower = 0> N_y_obs; // number of observed values
    int<lower = 0> N_y_mis; // number of missing values
    int<lower = 1> i_y_obs[N_y_obs];
    int<lower = 1> i_y_mis[N_y_mis];
    vector[N_y_obs] y_obs;
    int<lower = 1> n; // number of observations
    int<lower = 1> K; // number of covariates
    matrix[n,K] X[T] ; // real X[N,K,T]; //
    int<lower = 1, upper = n> m; // number of knots
    matrix[n, n] D;
    int<lower = 0, upper = n>  replic[m, niter + 2];  // replicates of the random selection of spatial locations

  vector[K]beta_mean; // informative prior betas
  vector[K]beta_sd; // informative prior betas

  real<lower=0> eta_mean; // informative prior sill
  real<lower=0> eta_sd; // informative prior sill

  real<lower=0> alpha_mean; // informative prior range
  real<lower=0> alpha_sd; // informative prior range

  real<lower=-1> phi_mean; // informative prior range
  real<lower=0> phi_sd; // informative prior range

  real<lower=0> sigma_nug_mean;
  real<lower=0> sigma_nug_sd;

  real<lower=0> sigma_mean;
  real<lower=0> sigma_sd;

  vector[m]  w_z_mean;
  vector[m]  w_z_sd;

  vector[n]  e_z_mean;
  vector[n]  e_z_sd;

  //real<lower=0> e_z_mean;



  }



  parameters {
    vector[K] beta;
    real<lower = 0> eta_;
    real <lower=-1, upper = 1> phi;

    real<lower = 0> sigma;
    real<lower = 0> sigma_nug;

    real<lower=0> alpha;

    vector[m] w_z;
    vector[n] e_z;


    vector[N_y_mis] y_mis; // declaring the missing y
    }

  transformed parameters {
    //vector[n] y;
    vector[n] w;
    vector[n] sigma_e_tilde;
    matrix[m, m] Cstar; // cov matrix knots
    vector[m] w_star;
    matrix[m, m] inv_Cstar;
    matrix[n, m] C_site_star; // cov between knots and sites.
    matrix[n, m] C_ss_inv_Cstar;
    real eta_sq;
    real sig_sq;
    vector[n * T] mu2;

    vector[n * T] y;
    vector[n] Y[T];

    vector[n] epsilon[T]; // error term
    vector[n] mu[T]; // mean


    matrix[m, m] D_star;  // distance between knots
    matrix[n, m] D_site_star; // distance between points and knots

    real<lower=0> var_nug; // nugget

     var_nug = sigma_nug ^ 2; // variance nugget

      D_star = D[replic[,iter_num], replic[,iter_num]];
      D_site_star = D[1:n, replic[,iter_num]];


    y[i_y_obs] = y_obs;
    y[i_y_mis] = y_mis;

    for (t in 1:T){
        Y[t] = y[((t - 1) * n + 1):(t * n)];
    }


    eta_sq = pow(eta_, 2);
    sig_sq = pow(sigma, 2);

    for (i in 1:(m-1)) { // cov matrix knots
      for (j in (i + 1):m) {
        Cstar[i, j] = eta_sq * exp(-D_star[i, j] / alpha);
        Cstar[j, i] = Cstar[i, j];
      }
    }

    for (k in 1:m) Cstar[k, k] = eta_sq + sig_sq;
    inv_Cstar = inverse(Cstar); // inverse of Cov matrix of knots
    w_star = cholesky_decompose(Cstar) * w_z;

    // latent gp at sample locations
    C_site_star = eta_sq * exp(-D_site_star / alpha); // covariance mat between sites and knots
    C_ss_inv_Cstar = C_site_star * inv_Cstar;
    w = C_site_star * inv_Cstar * w_star;

    // bias adjustment from Finley et al. 2009
    sigma_e_tilde = eta_sq + sig_sq - rows_dot_product(C_ss_inv_Cstar, C_site_star);
    for (i in 1:n) {
      w[i] = w[i] + e_z[i] * sqrt(sigma_e_tilde[i]);
    }

      mu[1] = X[1] * beta;
      epsilon[1] = Y[1] - mu[1];

    //mu = X * beta + w;
      for (t in 2:T){
          mu[t] = X[t] * beta;
          epsilon[t] = Y[t] - mu[t];
          mu[t] = mu[t] + phi * epsilon[t-1] + w; //
      }

      for (t in 1:T){
        mu2[((t - 1) * n + 1):(t * n)] = mu[t];
      }

  }

  model {
    eta_ ~ normal(eta_mean, eta_sd); // sill
    alpha ~ normal(alpha_mean, alpha_sd); // spat range

  //  beta ~ normal(0, 3);
    sigma ~ normal(sigma_mean, sigma_sd);

    sigma_nug ~ normal(sigma_nug_mean, sigma_nug_sd);

    w_z ~ normal(w_z_mean, w_z_sd);
    e_z ~ normal(e_z_mean, e_z_sd);

    target += normal_lpdf(y | mu2, sigma_nug);

    for (k in 1:K){
      beta[k] ~ normal(beta_mean[k], beta_sd[k]);
    }

    phi ~ normal(phi_mean, phi_sd);

  }

  generated quantities {
    vector[n*T] ypred;
    for (i in 1:(n*T)){
      ypred[i] = normal_rng(mu2[i], sigma_nug);
    }
  }

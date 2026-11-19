functions {
  /**
  * Return the log probability of a proper conditional autoregressive (CAR) prior 
  * with a sparse representation for the adjacency matrix
  *
  * @param phi Vector containing the parameters with a CAR prior
  * @param tau Precision parameter for the CAR prior (real)
  * @param alpha Dependence (usually spatial) parameter for the CAR prior (real)
  * @param W_sparse Sparse representation of adjacency matrix (int array)
  * @param n Length of phi (int)
  * @param W_n Number of adjacent pairs (int)
  * @param D_sparse Number of neighbors for each location (vector)
  * @param lambda Eigenvalues of D^{-1/2}*W*D^{-1/2} (vector)
  *
  * @return Log probability density of CAR prior up to additive constant
  */
  real sparse_pcar_lpdf(vector phi, real tau, real alpha,
  array[,] int W_sparse, vector D_sparse, vector lambda,
  int n, int W_n) {
    
    row_vector[n] phit_D; // phi' * D
    row_vector[n] phit_W; // phi' * W
    vector[n] ldet_terms;
    
    phit_D = (phi .* D_sparse)';
    phit_W = rep_row_vector(0, n);
    for (i in 1 : W_n) {
      phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
      phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
    }
    
    for (i in 1 : n) {
      ldet_terms[i] = log1m(alpha * lambda[i]);// log(1 - alpha * eigen_d2_w); Determinant term
    }
    return 0.5
    * (n * log(tau) + sum(ldet_terms)
    - tau * (phit_D * phi - alpha * (phit_W * phi)));
  }
}

data {
  
  int <lower=0> N; // Number of spatial areas
  int <lower=0> T; // Number of timepoints
  vector [N] pop;//population at risk
  vector [N] log_pop;
  int <lower=0> death_natural[N,T]; // natural deaths
  
  int <lower=0> N_edges; // number of edges
  array[N_edges, 2] int nodes; // Adjacency pairs
  vector [N] num;
  vector [N] eigen_d2_w;
  
  matrix [N,T] risk;
  matrix [N,T] lag1_risk;
  //matrix [N,T] lag2_risk;

  vector [N] lpopd;
  vector [N] mintry;
  
  vector [N] age_dist1;
  vector [N] age_dist2;
  vector [N] mcc;
  vector [N] cardio;
  
  vector [N] svi1;
  vector [N] svi2;
  vector [N] svi3;
  vector [N] svi4;//Percentile ranking for Housing Type/ Transportation theme
  //vector [N] svi_tot;
  // 
}


transformed data {
  // Pre-compute constant terms to avoid repeated computation
  matrix[N,T] lpopd_mat = rep_matrix(lpopd, T);
  matrix[N,T] age_dist1_mat = rep_matrix(age_dist1, T);
  matrix[N,T] age_dist2_mat = rep_matrix(age_dist2, T);
  matrix[N,T] cardio_mat = rep_matrix(cardio, T);

  matrix[N,T] mcc_mat = rep_matrix(mcc, T);

  matrix[N,T] svi1_mat = rep_matrix(svi1, T);
  matrix[N,T] svi2_mat = rep_matrix(svi2, T);
  matrix[N,T] svi3_mat = rep_matrix(svi3, T);
  matrix[N,T] svi4_mat = rep_matrix(svi4, T);
  
  //sensitivity check on f_infect and f_vacc
  //real f_infect=0.5;
  //real f_vacc=0.5;
}

parameters {
  
  real alpha_s; // spatial parameter
  
  vector <lower=0>[T] sigma_s; // Scale parameter spatial

  real b0;//poisson-gamma intercept
  vector [T] b_infect;

  real <lower=0, upper=1>f_infect;
  
  real b_mcc;//start from fixed
  real b_svi3;
  real b_svi4;
  real b1;
  real b2;
  real b3;
  real b4;
  real b5;
  real b6;
  
  real rho;

  matrix[N, T] ar_wn_s;

}

transformed parameters {
  
  vector [T] tau_s = 1/(sigma_s .* sigma_s);

  matrix[N, T] glmm_re2;

  matrix[N, T] log_eta;
  matrix[N, T] log_eta_nocovid;

  matrix[N,T] tot_risk_infect = f_infect * risk + (1-f_infect) * lag1_risk;

  glmm_re2[,1] = ar_wn_s[,1];
  for (t in 2:T) {
    glmm_re2[,t] = rho * glmm_re2[,t-1] + ar_wn_s[,t];//AR(1) pCAR for spatial RE
  }
  
  log_eta = b0 + tot_risk_infect .* rep_matrix(b_infect',N)
    + mcc_mat * b_mcc
    + svi3_mat * b_svi3
    + svi4_mat * b_svi4
    + lpopd_mat * b1 + age_dist1_mat * b2 + age_dist2_mat * b3
    + cardio_mat* b4 + svi1_mat * b5 + svi2_mat * b6
    + glmm_re2;

  //no covid circulation
  log_eta_nocovid = b0 
    + mcc_mat * b_mcc
    + svi3_mat * b_svi3
    + svi4_mat * b_svi4
    + lpopd_mat * b1 + + age_dist1_mat * b2 + age_dist2_mat * b3
    + cardio_mat* b4 + svi1_mat * b5 + svi2_mat * b6
    + glmm_re2;

}

model {
  
  alpha_s ~ beta(1,1);
  
  sigma_s ~ normal(0, 5);
  //sigma_s ~ lognormal(0, 10);
  //sigma_s ~ uniform(0, 10);
  
  // vague prior on b0 and b's;
  b0 ~ normal(0,10);
  b_infect ~ normal(0,10);
  
  b_mcc ~ normal(0,10);
  b_svi3 ~ normal(0,10);
  b_svi4 ~ normal(0,10);
  b1 ~ normal(0,10);
  b2 ~ normal(0,10);
  b3 ~ normal(0,10);
  b4 ~ normal(0,10);
  b5 ~ normal(0,10);
  b6 ~ normal(0,10);

  f_infect ~ beta(1,1);
  rho ~ beta(1,1);

  for (t in 1:T) {
    ar_wn_s[, t] ~ sparse_pcar(tau_s[t], alpha_s, nodes, num, eigen_d2_w, N, N_edges);

    death_natural[, t] ~ poisson_log(log_eta[, t] + log_pop);
  }
  
  
}

generated quantities {

  matrix[N,T] death_w_covid  = exp(log_eta        + rep_matrix(log_pop, T));
  matrix[N,T] death_wo_covid = exp(log_eta_nocovid + rep_matrix(log_pop, T));


  real total_excess_covid = sum(exp(log_eta + rep_matrix(log_pop, T))
                              - exp(log_eta_nocovid + rep_matrix(log_pop, T)));
  vector[T] excess_covid;

  for (t in 1:T) {

    excess_covid[t] = sum(exp(log_eta[, t] + log_pop)
                            - exp(log_eta_nocovid[, t] + log_pop));
  }

}



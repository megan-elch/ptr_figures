data{
  int<lower=0> n_g;
  int<lower=0> total_obs_mrna;
  int<lower=0> total_obs_protein; 

  int<lower=0> mrna_sum[total_obs_mrna];
  vector[total_obs_protein] protein_av;

  vector<lower=0>[total_obs_mrna] n_mrna;
  vector<lower=0>[total_obs_protein] n_protein;

  int<lower=0> mrna_gene_rec[total_obs_mrna];
  int<lower=0> protein_gene_rec[total_obs_protein];

  int<lower=0> mrna_celltype_rec[total_obs_mrna];
  int<lower=0> protein_celltype_rec[total_obs_protein];

  int<lower=0> mrna_pop_rec[total_obs_mrna];
  int<lower=0> protein_pop_rec[total_obs_protein];
}

transformed data{
  int<lower=0> n_ct;
  int<lower=0> n_pop_mrna;
  int<lower=0> n_pop_protein;

  n_ct = max(protein_celltype_rec);
  n_pop_mrna = max(mrna_pop_rec);
  n_pop_protein = max(protein_pop_rec);
}

parameters{
  matrix[n_g, n_ct] r_raw;
  matrix[n_g, n_ct] mu_raw;
  vector<lower=0>[n_g] scale_param;
  vector<lower=0>[n_g] gamma[n_pop_mrna];
  vector<lower=0>[n_ct] a_mrna[n_pop_mrna];
  matrix[n_pop_protein, n_g] kappa;
  matrix[n_pop_protein, n_ct] b_protein;
  vector<lower=0>[n_g] sigma_r;
  vector<lower=0>[n_g] sigma_mu;
  vector<lower=0>[n_pop_protein] sigma_p;
  vector<lower=0>[n_pop_protein] tau_p;
  vector<lower=0>[n_pop_mrna] phi_m;
}

 transformed parameters{
  matrix[n_g, n_ct] r;
  matrix[n_g, n_ct] mu;
  vector<lower=0>[total_obs_mrna] mrna_mean;
  vector<lower=0>[total_obs_mrna] phi_transform;
  vector[total_obs_protein] protein_mean;
  vector<lower=0>[total_obs_protein] protein_sd;
  
  mu = mu_raw .* rep_matrix(sqrt(sigma_mu), cols(mu_raw));
  r = r_raw .* rep_matrix(sqrt(sigma_r), cols(r_raw));
  
  for(i in 1:total_obs_mrna){
    phi_transform[i] = phi_m[mrna_pop_rec[i]];
    mrna_mean[i] = 1.0*n_mrna[i]*exp2(mu[mrna_gene_rec[i], mrna_celltype_rec[i]])*gamma[mrna_pop_rec[i], mrna_gene_rec[i]]*a_mrna[mrna_pop_rec[i], mrna_celltype_rec[i]];
  }
  
  for(i in 1:total_obs_protein){
    protein_mean[i] = inv(scale_param[protein_gene_rec[i]])*(mu[protein_gene_rec[i], protein_celltype_rec[i]] + r[protein_gene_rec[i], protein_celltype_rec[i]]) + kappa[protein_pop_rec[i], protein_gene_rec[i]] + b_protein[protein_pop_rec[i], protein_celltype_rec[i]];
    protein_sd[i] = sqrt(sigma_p[protein_pop_rec[i]]^2/n_protein[i] + tau_p[protein_pop_rec[i]]^2);
  }
}

model{
  to_vector(r_raw) ~ normal(0, 1);
  to_vector(mu_raw) ~ normal(0, 1);
  scale_param ~ normal(1, 0.5);

  for(k in 1:n_pop_mrna){
    gamma[k] ~ gamma(1, 0.25);
    a_mrna[k] ~ gamma(1, 0.25);
    pow(phi_m[k], -1) ~ normal(0, 1); // prior on 1/()
    target+= log(pow(phi_m[k], -2)); //jacobian adj if prior is on 1/()
   }
    
  sigma_mu ~ normal(0, 0.5);
  sigma_r ~ normal(0, 0.25);
  to_vector(sigma_p) ~ lognormal(0, 1);
  tau_p ~ lognormal(0, 1);
  mrna_sum ~ neg_binomial_2(mrna_mean, phi_transform);
  protein_av ~ normal(protein_mean, protein_sd);
  to_vector(kappa) ~ normal(0, 1);
  to_vector(b_protein) ~ normal(0, 1);
}

generated quantities{
  int<lower=0> mrna_sum_rep[total_obs_mrna];
  real protein_bar_rep[total_obs_protein];

  mrna_sum_rep = neg_binomial_2_rng(mrna_mean, phi_transform);
  protein_bar_rep = normal_rng(protein_mean, protein_sd);
} 
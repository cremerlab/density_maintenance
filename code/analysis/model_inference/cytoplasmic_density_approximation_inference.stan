data {
    int<lower=1>N_size; // Total number of volume measurements
    vector<lower=0>[N_size] surface_areas; 
    vector<lower=0>[N_size] volume;
    vector<lower=0>[N_size] size_lam;

    int<lower=1>N_dna; // Total number of measurements of DNA/Protein ratio
    vector<lower=0>[N_dna] dna_protein_ratio;

    int<lower=1>N_prot;
    vector<lower=0>[N_prot] prot_per_cell;
    vector<lower=0>[N_prot] prot_per_cell_lam;

    int<lower=1>N_ms; // Number of mass spec measurements
    vector<lower=0>[N_ms] ms_lam;
    vector<lower=0>[N_ms] phi_rib;
    vector<lower=0>[N_ms] phi_cyt;

}

parameters {
    real log_V_0;
    real<lower=0> k_V;
    real<lower=0> V_sigma;

    real<lower=0> SA_0;
    real<lower=0> k_SA;
    real<lower=0> SA_sigma;

    real<lower=0> theta_dna_mu;
    real<lower=0> theta_dna_sigma;

    real log_prot_intercept;
    real<lower=0> log_prot_slope;
    real<lower=0> log_prot_sigma;
}


model { 
    // Priors and likelihood for volume trend 
    log_V_0 ~ std_normal();
    k_V ~ std_normal();
    V_sigma ~ std_normal();
    log(volume) ~ normal(log_V_0 + k_V .* size_lam, V_sigma);

    // Literature surface area inference
    SA_0 ~ normal(0, 3);
    k_SA ~ normal(0, 2);
    SA_sigma ~ std_normal();
    surface_areas ~ normal(SA_0 + k_SA .* size_lam, SA_sigma);

    // Total protein inference
    log_prot_intercept ~ std_normal();
    log_prot_slope ~ normal(3, 2);
    log_prot_sigma ~ normal(0, 0.1);
    log(prot_per_cell) ~ normal(log_prot_intercept + log_prot_slope .* prot_per_cell_lam, log_prot_sigma);

    // DNA to protein ratio 
    theta_dna_mu ~ normal(0, 0.1);
    theta_dna_sigma ~ normal(0, 1);
    dna_protein_ratio ~ normal(theta_dna_mu, theta_dna_sigma); 
}

generated quantities {
    vector[N_ms] M_prot_tot = exp(log_prot_intercept + log_prot_slope .* ms_lam);
    vector[N_ms] V_cyt = exp(log_V_0 + k_V .* ms_lam) - 0.0246 .* (SA_0 + k_SA .* ms_lam);
    vector[N_ms] M_RNA = phi_rib .* M_prot_tot / 0.4558;
    vector[N_ms] M_DNA = theta_dna_mu .* M_prot_tot;
    vector[N_ms] M_prot_cyt = phi_cyt .* M_prot_tot;
    vector[N_ms] rho_cyt = (M_RNA + M_DNA + M_prot_cyt) ./ V_cyt;
}
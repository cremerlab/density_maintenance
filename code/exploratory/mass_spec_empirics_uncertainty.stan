data { 
    // Size data
    int<lower=1> N_size;
    vector<lower=0>[N_size] size_lam;
    vector<lower=0>[N_size] widths;
    vector<lower=0>[N_size] lengths;


    // Protein data 
    int<lower=0> N_prot;
    vector<lower=0>[N_prot] prot_lam;
    vector<lower=0>[N_prot] prot_per_cell;

    // Mass spec data
    int<lower=0> N_ms;
    vector<lower=0>[N_ms] ms_lam;
    vector<lower=0>[N_ms] phi_mem;
    vector<lower=0>[N_ms] phi_peri;

    // PhiRb data
    int<lower=0> N_phi;
    vector<lower=0>[N_phi] phi_lam;
}

parameters {
    real<lower=0> width_slope; 
    real width_intercept; 
    real<lower=0> width_sigma; 
    real<lower=0> length_slope; 
    real length_intercept; 
    real<lower=0> length_sigma; 
    real log_prot_slope;
    real log_prot_intercept; 
    real log_prot_sigma;
}

transformed parameters {
    vector<lower=0>[N_ms] ms_prot = exp(log_prot_intercept + log_prot_slope .* ms_lam);
    vector<lower=0>[N_ms] ms_width = width_intercept + width_slope .* ms_lam;
    vector<lower=0>[N_ms] ms_length = length_intercept + length_slope .* ms_lam;
    vector<lower=0>[N_ms] ms_surface_area = pi() .* ms_length .* ms_width;
    vector<lower=0>[N_ms] m_mem = phi_mem .* ms_prot; 
    vector<lower=0>[N_ms] m_peri = phi_peri .* ms_prot;
    vector<lower=0>[N_ms] rho_mem = m_mem ./ (2 .* ms_surface_area);
    vector<lower=0>[N_ms] rho_peri = m_peri ./ (ms_surface_area .* 0.0246);
    vector<lower=0>[N_phi] phi_width = width_intercept + width_slope .* phi_lam;
}

model  {
    width_slope ~ std_normal();
    width_intercept ~ std_normal();
    width_sigma ~ std_normal();
    widths ~ normal(width_intercept + width_slope .* size_lam, width_sigma);

    length_slope ~ std_normal();
    length_intercept ~ std_normal();
    length_sigma ~ std_normal();
    lengths ~ normal(length_intercept + length_slope .* size_lam, length_sigma);

    log_prot_intercept ~ std_normal();
    log_prot_slope ~ std_normal();
    log_prot_sigma ~ std_normal(); 
    log(prot_per_cell) ~ normal(log_prot_intercept + log_prot_slope .* prot_lam, log_prot_sigma);
}

generated quantities { 
    vector<lower=0>[N_size] width_ppc;
    vector<lower=0>[N_size] length_ppc;
    vector<lower=0>[N_prot] prot_per_cell_ppc;

    for (i in 1:N_size) {
        width_ppc[i] = normal_rng(width_intercept + width_slope * size_lam[i], width_sigma);
        length_ppc[i] = normal_rng(length_intercept + length_slope * size_lam[i], length_sigma);

    }
    for (i in 1:N_prot) {
        prot_per_cell_ppc[i] = exp(normal_rng(log_prot_intercept + log_prot_slope * prot_lam[i], log_prot_sigma));
    }
}
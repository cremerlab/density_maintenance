data {
    // Literature size data.
    int<lower=1> N_size;
    vector<lower=0>[N_size] size_lam;
    vector<lower=0>[N_size] surface_areas;

    // Protein-per-cell data 
    int<lower=0> N_prot;
    vector<lower=0>[N_prot] prot_lam;
    vector<lower=0>[N_prot] prot_per_cell;

    // Mass spec data
    int<lower=0> N_ms;
    vector<lower=0>[N_ms] ms_lam;
    vector<lower=0>[N_ms] phi_mem_outer;
    vector<lower=0>[N_ms] phi_mem_inner;
}

parameters {
    // Shape parameters for calculating rho_mem from the mass spec data.
    real<lower=0> surface_area_slope;
    real<lower=0> surface_area_intercept;
    real<lower=0> surface_area_sigma;
    real log_prot_intercept;
    real<lower=0> log_prot_slope;
    real<lower=0> log_prot_sigma;

}



model { 
    // Literature surface area inference
    surface_area_intercept ~ normal(0, 3);
    surface_area_slope ~ normal(0, 2);
    surface_area_sigma ~ std_normal();
    surface_areas ~ normal(surface_area_intercept + surface_area_slope .* size_lam, surface_area_sigma);

    // Protein per cell
    log_prot_intercept ~ std_normal();
    log_prot_slope ~ normal(3, 2);
    log_prot_sigma ~ normal(0, 0.1);
    log(prot_per_cell) ~ normal(log_prot_intercept + log_prot_slope .* prot_lam, log_prot_sigma);

}

generated quantities{
    vector<lower=0>[N_ms] ms_total_prot = exp(log_prot_intercept + log_prot_slope .* ms_lam);
    vector<lower=0>[N_ms] ms_surface_area = surface_area_intercept + surface_area_slope .* ms_lam; 
    vector<lower=0>[N_ms] ms_rho_mem_outer = phi_mem_outer .* ms_total_prot ./ (ms_surface_area);
    vector<lower=0>[N_ms] ms_rho_mem_inner = phi_mem_inner .* ms_total_prot ./ (ms_surface_area);
}
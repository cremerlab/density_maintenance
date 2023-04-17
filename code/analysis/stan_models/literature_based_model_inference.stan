data {
    // Dimensional information
    int<lower=1> N_size;
    int<lower=1> N_prot;
    int<lower=1> N_mass_spec;
    real<lower=0> delta;

    // Observed data
    vector<lower=0>[N_size] widths;
    vector<lower=0>[N_size] lengths;
    vector<lower=0>[N_size] volumes;
    vector<lower=0>[N_size] size_lam;
    vector<lower=0>[N_prot] prot_per_cell;
    vector<lower=0>[N_prot] prot_lam;
    vector<lower=0>[N_mass_spec] phi_mem;
    vector<lower=0>[N_mass_spec] phi_peri;
    vector<lower=0>[N_mass_spec] ms_lam;
}

parameters { 
    real<lower=0> w_min;
    real<lower=0> w_slope;
    real<lower=0> w_sigma;
    real<lower=0> ell_sigma;
    real<lower=0> vol_sigma;
    real<lower=0> prot_sigma;
    real<lower=1> alpha;
    real<lower=0> rho_prot;
    real<lower=0> m_peri;
    real<lower=0> phi_mem_mu;
    real<lower=0> phi_mem_sigma;
    real<lower=0> phi_peri_sigma;
    
}

model {

   w_min ~ normal(0, 0.1);
   w_slope ~ std_normal();
   w_sigma ~ std_normal();
   alpha ~ normal(1, 3);
   rho_prot ~ normal(200, 100);
   m_peri ~ normal(0, 100);
   phi_mem_mu ~ beta(2, 10);

   // Likelihoods for size measurements
   widths ~ normal(w_min + w_slope .* size_lam, w_sigma);
   lengths ~ normal(alpha .* (w_min + w_slope .* size_lam), ell_sigma);
   volumes ~ normal((pi()/12) * (w_min + w_slope .* size_lam).^3 * (3 * alpha - 1), vol_sigma);

   // Likelihoods for protein measurements
   prot_per_cell ~ normal(rho_prot * (pi()/12) .* (w_min + w_slope .* prot_lam).^3 * (3 * alpha - 1), prot_sigma);

   // Likelihoods based on protein measurements
   phi_mem ~ normal(phi_mem_mu, phi_mem_sigma);
   phi_peri ~ normal(m_peri / (rho_prot * (pi()/12) .* (w_min + w_slope .* ms_lam).^3 * (3 * alpha - 1)), phi_peri_sigma);
}

generated quantities {
    // PPC for fit quantities
    vector[N_size] w_rep;
    vector[N_size] ell_rep;
    vector[N_size] vol_rep;
    vector[N_prot] prot_per_cell_rep;
    vector[N_mass_spec] phi_mem_rep;
    vector[N_mass_spec] phi_peri_rep;
    vector[N_mass_spec] rho_peri;
    vector[N_mass_spec] rho_mem;

    for (i in 1:N_size) {
        w_rep[i] = normal_rng(w_min + w_slope * size_lam[i], w_sigma);
        ell_rep[i] = normal_rng(alpha * w_rep[i], ell_sigma);
        vol_rep[i] = normal_rng((pi()/12) * (w_min + w_slope * size_lam[i])^3 * (3 * alpha - 1), vol_sigma);
    }

    for (i in 1:N_prot) {
        prot_per_cell_rep[i] = normal_rng(rho_prot * (pi()/12) * (w_min + w_slope * prot_lam[i])^3 * (3 * alpha - 1), prot_sigma);
    }

    for (i in 1:N_mass_spec) {
        phi_mem_rep[i] = normal_rng(phi_mem_mu, phi_mem_sigma);
        phi_peri_rep[i] = normal_rng(m_peri / (rho_prot * (pi()/12) * (w_min + w_slope * ms_lam[i])^3 * (3 * alpha - 1)), phi_peri_sigma);
        rho_mem[i] = phi_mem_mu * rho_prot * (pi()/12) * (w_min + w_slope * ms_lam[i])^3 * (3 * alpha - 1) / (pi() * alpha * (w_min + w_slope * ms_lam[i])^2);
        rho_peri[i] = m_peri / (pi() * alpha * delta * (w_min + w_slope * ms_lam[i])^2);
    }
}


data {
    // Dimensional information
    int<lower=1> N_size;
    int<lower=1> N_prot;
    int<lower=1> N_mass_spec;
    int<lower=1> N_sim;
    real<lower=0> delta;
    vector<lower=0>[N_sim] lam_sim;

    // Observed data
    vector<lower=0>[N_size] widths;
    vector<lower=0>[N_size] lengths;
    vector<lower=0>[N_size] volumes;
    vector<lower=0>[N_size] size_lam;
    vector<lower=0>[N_prot] prot_per_cell;
    vector<lower=0>[N_prot] prot_lam;
    vector<lower=0>[N_mass_spec] phi_mem;
    vector<lower=0>[N_mass_spec] phi_peri;
    vector<lower=0>[N_mass_spec] m_peri_meas;
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
    real<lower=0> m_peri_mu;
    real<lower=0> m_peri_sigma;
    real<lower=0> phi_mem_mu;
    real<lower=0> phi_mem_sigma;
    real<lower=0> phi_peri_sigma;
    
}

model {

   w_min ~ normal(0, 0.1);
   w_slope ~ std_normal();
   w_sigma ~ normal(0, 0.1);
   vol_sigma ~ normal(0, 0.1);
   ell_sigma ~ normal(0, 0.1);
   phi_peri_sigma ~ normal(0, 0.1);
   phi_mem_sigma ~ normal(0, 0.1);
   alpha ~ normal(1, 3);
   rho_prot ~ normal(200, 100);
   m_peri_mu ~ normal(0, 100);
   m_peri_sigma ~ std_normal();
   phi_mem_mu ~ beta(2, 10);

   // Likelihoods for size measurements
   widths ~ normal(w_min + w_slope .* size_lam, w_sigma);
   lengths ~ normal(alpha .* (w_min + w_slope .* size_lam), ell_sigma);
   log(volumes) ~ normal(log((pi()/12) * (w_min + w_slope .* size_lam).^3 * (3 * alpha - 1)), vol_sigma);

   // Likelihoods for protein measurements
   log(prot_per_cell) ~ normal(log(rho_prot * (pi()/12) .* (w_min + w_slope .* prot_lam).^3 * (3 * alpha - 1)), prot_sigma);

   // Likelihoods based on protein measurements
   phi_mem ~ normal(phi_mem_mu, phi_mem_sigma);
   phi_peri ~ normal(m_peri_mu / (rho_prot * (pi()/12) .* (w_min + w_slope .* ms_lam).^3 * (3 * alpha - 1)), phi_peri_sigma);
   m_peri_meas ~ normal(m_peri_mu, m_peri_sigma);
}

generated quantities {
    // PPC for fit quantities
    vector[N_sim] w_rep;
    vector[N_sim] ell_rep;
    vector[N_sim] vol_rep;
    vector[N_sim] prot_per_cell_rep;
    vector[N_sim] phi_mem_rep;
    vector[N_sim] phi_peri_rep;
    vector[N_sim] rho_peri;
    vector[N_sim] rho_mem;
    vector[N_sim] m_peri;

    for (i in 1:N_sim) {
        m_peri[i] = normal_rng(m_peri_mu, m_peri_sigma);
        w_rep[i] = normal_rng(w_min + w_slope * lam_sim[i], w_sigma);
        ell_rep[i] = normal_rng(alpha * w_rep[i], ell_sigma);
        vol_rep[i] = exp(normal_rng(log((pi()/12) * (w_min + w_slope * lam_sim[i])^3 * (3 * alpha - 1)), vol_sigma));
        prot_per_cell_rep[i] = exp(normal_rng(log(rho_prot * (pi()/12) * w_rep[i]^3 * (3 * alpha - 1)), prot_sigma));    
        phi_mem_rep[i] = normal_rng(phi_mem_mu, phi_mem_sigma);
        phi_peri_rep[i] = normal_rng(m_peri[i] / (rho_prot * (pi()/12) * w_rep[i]^3 * (3 * alpha - 1)), phi_peri_sigma);
        rho_mem[i] = phi_mem_rep[i] * rho_prot * (pi()/12) * w_rep[i]^3 * (3 * alpha - 1) / (pi() * alpha * w_rep[i]^2);
        rho_peri[i] = m_peri_mu / (pi() * alpha * delta * w_rep[i]^2);
        
    }
}

